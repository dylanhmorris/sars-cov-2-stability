#!/usr/bin/env Rscript

###################################
## plot fitted exponential decay
## rates for viruses
####################################

script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'magrittr',   # for pipe operator %>%
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'ggplot2',    # for plotting
    'modelr',     # for data_grid()
    'tidyr',      # for crossing()
    'cowplot',    # publication ready ggplot
    'scales',     # for trans_breaks(), etc.
    'forcats'    # for fct_relevel()
)   

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                 character.only=TRUE))
}


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly=TRUE)
data_path <- args[1]
results_path <- args[2]
plotting_params_path <- args[3]
logged <- as.logical(args[4])
material <- as.character(args[5])
outpath <- args[6]


## read data / style files
source(plotting_params_path) # shared project style

chains <- readRDS(results_path)

dat <- read_data_for_plotting(data_path) # function in plotting_style.R

#################################
## overall plot styling
#################################
set.seed(989327) # reproducible! (since we use random draws)


##################################################
## calculate posterior draws for regression lines
##################################################
experiment_dat_virus <- dat[dat$material == material, ]
print(experiment_dat_virus)

experiment_dat_virus <- experiment_dat_virus %>%
    summarise(detection_limit = 10^first(detection_limit_log10_titer))
print(experiment_dat_virus)

mat_dat <- dat[dat$material == material, ]

if(material != "Aerosols") {
    scaling <- 1
    ylab_expression <- expression("titer (TCID"[50] * "/mL media)")
    max_x <- mat_dat %>%
        group_by(trial_unique_id, replicate) %>%
        filter(log10_titer == detection_limit_log10_titer) %>%
        select(time, trial_unique_id, replicate) %>%
        summarise(min_time = min(time)) %>%
        ungroup() %>%
        select(min_time) %>%
        max()
} else {
    scaling <- 10 / 3 ## convert to tcid50/L air
    ylab_expression <- expression("titer (TCID"[50] * "/L air)")
    max_x <- max(mat_dat$time)
}
    

print(max_x)

plot_times <- dat %>%
    data_grid(time = seq_range(c(0, max_x),
                               n = fineness))

print(material)

## get needed draws and add human readable names
tidy_draws <- chains %>%
    spread_draws(decay_rate[trial_unique_id],
                 intercept[trial_unique_id][replicate]) 
tidy_draws <- dat %>%
    select(trial_unique_id,
           grouping_name,
           replicate,
           virus,
           material) %>% 
    inner_join(tidy_draws,
               by = c("trial_unique_id", "replicate")) %>%
    ungroup()

tidy_draws <- tidy_draws[tidy_draws$material == material,]

print(tidy_draws)

## draw n_lines random regression lines
func_samples <- tidy_draws %>%
    group_by(trial_unique_id,
             replicate) %>%
    sample_n(n_lines) %>%
    ungroup()

print(func_samples)

## annotate lines so that each
## has a unique id for ggplot overplotting
## (else two lines from the same draw but
## different replicates can get confused
## with each other)
func_samples <- func_samples %>%
    mutate(line_id = as.numeric(rownames(func_samples)))

## cross product decay_rates with x (time) values
## and calculate y (titer) values
to_plot <- func_samples %>% crossing(plot_times)

to_plot <- to_plot %>%
    mutate(predicted_titer = scaling * 10^(intercept - decay_rate * time))

dat <- dat[dat$material == material, ]

max_titer <- scaling * max(dat$titer)

## plot either logged or not
if (logged) {
    ytrans <- 'log10'
    ylim <- c(10^(lowest_log_titer - 0.5),
              max_titer)
    ybreaks <- trans_breaks('log10', function(x) 10^x)
    yformat <- trans_format(ytrans,
                            math_format(10^.x))
    max_nonzero_time <- to_plot %>%
        filter(log10(predicted_titer) > lowest_log_titer) %>%
        select(time) %>% max()
    xlim <- c(0, max_x)


} else {
    ytrans <- 'identity'
    ylim <- c(0, max_titer)
    ybreaks <- pretty_breaks(n = n_breaks)(ylim)
    yformat <- trans_format(ytrans,
                            sci10_formatter)
    xlim <- c(0, max_x)
}


fit_fig <- to_plot %>%
    ggplot(aes(
        x = time,
        y = predicted_titer,
        color = virus,
        group = line_id)) +
    geom_line(alpha = line_alpha,
              size = line_size) +
    scale_colour_manual(values = unlist(virus_colors)) +
    geom_point(aes(x = time,
                   y = scaling * 10^log10_titer,
                   group = replicate),
               data = dat,
               color = pointborder,
               fill = pointfill,
               alpha = pointalpha,
               size = pointsize,
               stroke = pointstroke) +
    geom_hline(
        data = experiment_dat_virus,
        aes(yintercept = scaling * detection_limit),
        linetype = detection_linestyle) + 
    scale_y_continuous(trans = ytrans,
                       breaks = ybreaks,
                       labels = yformat) +
    coord_cartesian(ylim = ylim,
                    xlim = xlim) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_project(base_size = 50) +
    theme(legend.position = "none") +
    facet_grid(vars(virus),
               vars(replicate))


## style the figure a bit

fit_fig <- fit_fig +
    xlab("time (hrs)") +
    ylab(ylab_expression) +
    ggtitle(material)
                                       
save_plot(outpath,
          fit_fig,
          base_height = 15,  # set in plotting_style.R
          base_asp = 2.15)     # set in plotting_style.R
warnings()
