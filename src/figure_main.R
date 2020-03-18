#!/usr/bin/env Rscript

########################################
## plot main figure, which contains raw
## data, fitted exponential decay rates
## for viruses, and violin plots of
## half-life posterior distributions
#######################################

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
outpath <- args[4]

## read data / style files
source(plotting_params_path) # shared project style

chains <- readRDS(results_path)

dat <- read_data_for_plotting(data_path) # function in plotting_style.R

#################################
## overall plot styling
#################################
set.seed(989327) # reproducible! (since we use random draws)
text_size = 30
detection_linesize = 0.75
titer_ylab <- expression("titer (TCID"[50] * "/mL media)")
aer_titer_ylab <- expression("titer (TCID"[50] * "/L air)")
bar_width <- 0.75 # width for bars in barplot
bar_cap_width <- 0.5 # cap width for error bar caps in barplot
bar_err_size <- 0.75 # line thickness for error bars in barplot
scatter_err_size <- 1 # line thickness for error bars in scatterplot
dodge_width <- 0.5 # how much to dodge dots
convert_mL_media_to_L_air <- 10 / 3
##################################################
## calculate posterior draws for regression lines
##################################################


surface_dat <- dat %>% filter(material != "Aerosols")
aer_dat <- dat %>% filter(material == "Aerosols")

cat('calculating max x values...\n')
surface_max_x <- surface_dat %>%
    group_by(trial_unique_id, replicate) %>%
    filter(log10_titer == detection_limit_log10_titer) %>%
    select(time, trial_unique_id, replicate) %>%
    summarise(min_time = min(time)) %>%
    ungroup() %>%
    select(min_time) %>%
    max()

aer_max_x <- max(aer_dat$time)

surface_plot_times <- surface_dat %>%
    data_grid(time = seq_range(c(0, surface_max_x),
                               n = fineness))
aer_plot_times <- aer_dat %>%
    data_grid(time = seq_range(c(0, aer_max_x),
                               n = fineness))


## get needed draws and add human readable names
tidy_draws <- chains %>%
    spread_draws(decay_rate[trial_unique_id],
                 intercept[trial_unique_id][replicate])

tidy_draws <- dat %>%
    distinct(trial_unique_id,
             .keep_all = TRUE) %>%
    select(trial_unique_id,
           grouping_name,
           virus,
           material) %>% 
    inner_join(tidy_draws,
               by = "trial_unique_id")

## sort facets by the posterior median half life
## of the first virus in our virus order
## (HCoV-19, SARS-CoV-1, MERS-CoV)
## in the plotted data
sorting_virus <- virus_order[virus_order %in% tidy_draws$virus][1]
cat("sorting panels by ", sorting_virus, "half-life...\n")

ordering <- tidy_draws %>%
    filter(virus == sorting_virus,
           material != "Aerosols") %>%
    group_by(material) %>%
    summarise(med = median(-decay_rate)) %>%
    arrange(med)

biggest_first_below_limit_time <- dat %>%
    filter(log10_titer <= detection_limit_log10_titer) %>%
    group_by(material, virus, replicate) %>%
    summarise(first_below = min(time)) %>%
    ungroup() %>%
    select(first_below) %>%
    max()

tidy_draws_aer <- tidy_draws %>%
    filter(material == "Aerosols")
tidy_draws_aer$material <- factor(
    tidy_draws_aer$material)

tidy_draws_surface <- tidy_draws %>%
    filter(material != "Aerosols")
tidy_draws_surface$material <- factor(
    tidy_draws_surface$material,
    levels = ordering$material)



surface_dat$material <- factor(
    surface_dat$material,
    levels = ordering$material)
aer_dat$material <- factor(
    aer_dat$material)

surface_dat$fct_time <- factor(
    surface_dat$time,
    levels = sort(unique(surface_dat$time)))

aer_dat$fct_time <- factor(
    aer_dat$time,
    levels = sort(unique(aer_dat$time)))

## calculate y limits
max_titer <- max(dat$titer)
max_aer_titer <- convert_mL_media_to_L_air * max(aer_dat$titer)
ytrans <- 'log10'

surface_ylim <- c(10^(lowest_log_titer - 0.5),
                  max_titer)

aer_ylim <- c(10^(lowest_log_titer - 0.5),
              max_aer_titer)

ybreaks <- trans_breaks('log10', function(x) 10^x)
yformat <- trans_format(ytrans,
                        math_format(10^.x))

## calculate detection limits
experiment_dat_surface <- surface_dat %>%
    group_by(material) %>%
    summarise(detection_limit = 10^first(detection_limit_log10_titer))
experiment_dat_aer <- aer_dat %>%
    group_by(material) %>%
    summarise(detection_limit = 10^first(detection_limit_log10_titer))

experiment_dat_virus_surface <- surface_dat %>%
    group_by(material, virus) %>%
    summarise(detection_limit = 10^first(detection_limit_log10_titer))

experiment_dat_virus_aer <- aer_dat %>%
    group_by(material, virus) %>%
    summarise(detection_limit = 10^first(detection_limit_log10_titer))


###################################
## plot panel showing raw surface
## data
###################################
cat('plotting raw data...\n')


print(biggest_first_below_limit_time)

raw_panel_surface <- surface_dat %>%
    filter(time <= biggest_first_below_limit_time) %>%
    ggplot(aes(
        x = fct_time,
        y = 10^log10_titer,
        group = virus,
        fill = virus,
        color = virus)) +
    geom_point(
        stat = "summary",
        fun.y = "mean",
        size = pointsize * 0.8,
        position = position_dodge(width = dodge_width)) +
    stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        size = bar_err_size,
        width = bar_cap_width,
        position = position_dodge(width = dodge_width)) +
    geom_hline(
        data = experiment_dat_surface,
        aes(yintercept = detection_limit),
        linetype = detection_linestyle,
        size = detection_linesize) + 
    scale_fill_manual(values = unlist(virus_colors)) +
    scale_color_manual(values = unlist(virus_colors)) +
    scale_y_continuous(trans = ytrans,
                       breaks = ybreaks,
                       labels = yformat) +
    coord_cartesian(ylim = surface_ylim) +
    xlab("time (hrs)") +
    ylab(titer_ylab) +
    facet_grid(cols = vars(material),
               drop = TRUE) +
    theme_project(base_size = text_size) +
    theme(legend.position = "none")


raw_panel_aer <- aer_dat %>%
    filter(time <= biggest_first_below_limit_time) %>%
    ggplot(aes(
        x = fct_time,
        y = convert_mL_media_to_L_air * 10^(log10_titer),
        group = virus,
        fill = virus,
        color = virus)) +
    geom_point(
        stat = "summary",
        fun.y = "mean",
        size = pointsize * 0.8,
        position = position_dodge(width = dodge_width)) +
    stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        size = bar_err_size,
        width = bar_cap_width,
        position = position_dodge(width = dodge_width)) +
    geom_hline(
        data = experiment_dat_aer,
        aes(yintercept = detection_limit * convert_mL_media_to_L_air),
        linetype = detection_linestyle,
        size = detection_linesize) + 
    scale_fill_manual(values = unlist(virus_colors)) +
    scale_color_manual(values = unlist(virus_colors)) +
    scale_y_continuous(trans = ytrans,
                       breaks = ybreaks,
                       labels = yformat) +
    coord_cartesian(ylim = aer_ylim) +
    xlab("time (hrs)") +
    ylab(aer_titer_ylab) +
    facet_grid(cols = vars(material),
               drop = TRUE) +
    theme_project(base_size = text_size) +
    theme(legend.position = c(0.65, 0.2),
          legend.margin = margin(t = -10,
                                 l = 0,
                                 r = 0,
                                 b = 0),
          legend.title = element_blank(),
          legend.box.background = element_rect(color = "black",
                                               size = 2))


###################################
## plot panel showing fit of
## regression lines to real data
###################################
cat('plotting regression lines...\n')
## draw n_lines random regression lines
func_samples_surface <- tidy_draws_surface %>%
    group_by(trial_unique_id,
             replicate) %>%
    sample_n(n_lines) %>%
    ungroup()
func_samples_aer <- tidy_draws_aer %>%
    group_by(trial_unique_id,
             replicate) %>%
    sample_n(n_lines) %>%
    ungroup()

## annotate lines so that each
## has a unique id for ggplot overplotting
## (else two lines from the same draw but
## different replicates can get confused
## with each other)
func_samples_surface <- func_samples_surface %>%
    mutate(line_id = as.numeric(rownames(func_samples_surface)))
func_samples_aer <- func_samples_aer %>%
    mutate(line_id = as.numeric(rownames(func_samples_aer)))

## cross product decay_rates with x (time) values
## and calculate y (titer) values
cat('setting up x values...\n')

to_plot_surface <- func_samples_surface %>%
    crossing(surface_plot_times)

to_plot_aer <- func_samples_aer %>%
    crossing(aer_plot_times)


to_plot_surface <- to_plot_surface %>%
    mutate(predicted_titer = 10^(intercept - decay_rate * time))
to_plot_aer <- to_plot_aer %>%
    mutate(predicted_titer = convert_mL_media_to_L_air * 10^(intercept - decay_rate * time))


max_nonzero_time <- to_plot_surface %>%
    filter(log10(predicted_titer) > lowest_log_titer) %>%
    select(time) %>% max()
surface_xlim <- c(0, max_nonzero_time)
aer_xlim <- c(0, aer_max_x)
print(aer_xlim)
aer_jitwid <- 3/100

fit_panel_surface <- to_plot_surface %>%
    ggplot(aes(
        x = time,
        y = predicted_titer,
        color = virus,
        group = line_id)) +
    geom_line(alpha = line_alpha,
              size = line_size) +
    scale_colour_manual(values = unlist(virus_colors)) +
    geom_point(aes(x = time,
                   y = 10^(log10_titer),
                   group = trial_unique_id),
               data = surface_dat,
               color = pointborder,
               fill = pointfill,
               alpha = pointalpha,
               size = pointsize,
               stroke = pointstroke,
               position = position_jitter(
                   width = jitwid,
                   height = jith,
                   seed = 5)) +
    geom_hline(
        data = experiment_dat_virus_surface,
        aes(yintercept = detection_limit),
        linetype = detection_linestyle,
        size = detection_linesize) + 
    scale_y_continuous(trans = ytrans,
                       breaks = ybreaks,
                       labels = yformat) +
    coord_cartesian(ylim = surface_ylim,
                    xlim = surface_xlim) +
    facet_grid(vars(virus),
               vars(material),
               drop = TRUE)

## styling: no facet labels because is background plot
fit_panel_surface <- fit_panel_surface +
    theme_project(base_size = text_size) +
    xlab("time (hrs)") +
    ylab(titer_ylab)


fit_panel_aer <- to_plot_aer %>%
    ggplot(aes(
        x = time,
        y = predicted_titer,
        color = virus,
        group = line_id)) +
    geom_line(alpha = line_alpha,
              size = line_size) +
    scale_colour_manual(values = unlist(virus_colors)) +
    geom_point(aes(x = time,
                   y = convert_mL_media_to_L_air * 10^(log10_titer),
                   group = trial_unique_id),
               data = aer_dat,
               color = pointborder,
               fill = pointfill,
               alpha = pointalpha,
               size = pointsize,
               stroke = pointstroke,
               position = position_jitter(
                   width = aer_jitwid,
                   height = jith,
                   seed = 5)) +
    geom_hline(
        data = experiment_dat_virus_aer,
        aes(yintercept = detection_limit * convert_mL_media_to_L_air),
        linetype = detection_linestyle,
        size = detection_linesize) + 
    scale_y_continuous(trans = ytrans,
                       breaks = ybreaks,
                       labels = yformat) +
    coord_cartesian(ylim = aer_ylim,
                    xlim = aer_xlim) +
    facet_grid(rows = vars(virus),
               drop = TRUE)


## styling: no facet labels because is background plot
fit_panel_aer <- fit_panel_aer +
    theme_project(base_size = text_size) +
    xlab("time (hrs)") +
    ylab(aer_titer_ylab) +




##################################################
## calculate posterior draws for half lives
##################################################

cat('extracting posterior draws...\n')

half_life <- tidy_draws %>%
    mutate(half_life = log10(2) / decay_rate)

cat('Plotting half-lives...\n')

quants <- half_life %>%
    group_by(grouping_name) %>%
    summarise(
        q025 = quantile(half_life, 0.025),
        q50 = quantile(half_life, 0.5),
        q975 = quantile(half_life, 0.975),
        q995 = quantile(half_life, 0.995))

## sort facets by the posterior median half life
## of the first virus in our virus order
## (HCoV-19, SARS-CoV-1, MERS-CoV)
## in the plotted data
sorting_virus <- virus_order[virus_order %in% half_life$virus][1]
cat("sorting panels by ", sorting_virus, "half-life...\n")

ordering <- half_life %>%
    filter(virus == sorting_virus) %>%
    group_by(material) %>%
    summarise(med = median(-decay_rate)) %>%
    arrange(med)

half_life$material <- factor(
    half_life$material,
    levels = ordering$material)

dat$material <- factor(
    dat$material,
    levels = ordering$material)

## establish axis limits and ticks
hl_max_y <- max(max(quants$q995), 10)
hl_ylim <- c(0, hl_max_y)
hl_yticks <- seq(0, 10, 2)

half_life_surface <- half_life %>%
    filter(material != "Aerosols") %>%
    ggplot(aes(
        x = virus,
        y = half_life,
        fill = virus)) +
    geom_eye(size = text_size / 2) +
    scale_fill_manual(values = virus_colors) +
    facet_grid(cols = vars(material),
               switch = "y",
               drop = TRUE) +
    coord_cartesian(ylim = hl_ylim) +
    scale_y_continuous(breaks = hl_yticks,
                       expand = c(0, 0))

## styling -- 
half_life_surface <- half_life_surface +
    theme_project(base_size = text_size) +
    ylab("half-life (hrs)") +
    xlab("virus")

half_life_aer <- half_life %>%
    filter(material == "Aerosols") %>%
    ggplot(aes(
        x = virus,
        y = half_life,
        fill = virus)) +
    geom_eye(size = text_size / 2) +
    scale_fill_manual(values = virus_colors) +
    facet_grid(cols = vars(material),
               switch = "y",
               drop = TRUE) +
    coord_cartesian(ylim = hl_ylim) +
    scale_y_continuous(breaks = hl_yticks,
                       expand = c(0, 0))

## styling -- 
half_life_aer <- half_life_aer +
    theme_project(base_size = text_size) +
    ylab("half-life (hrs)") +
    xlab("virus")


####################################
## compose full figure from panels
####################################

upper_panel_theme <- theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = text_size),
    strip.placement = "outside",
    strip.switch.pad.grid = unit("0.5", "in"))

## no facet labels for lower panel because 
## can just read the upper ones
lower_panel_theme <- theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = text_size))

raw_panel_aer <- raw_panel_aer + upper_panel_theme
raw_panel_surface <- raw_panel_surface + upper_panel_theme

fit_panel_aer <- fit_panel_aer + lower_panel_theme +
    theme(strip.text.y = element_blank()) 
fit_panel_surface <- fit_panel_surface + lower_panel_theme

half_life_aer <- half_life_aer + lower_panel_theme
half_life_surface <- half_life_surface + lower_panel_theme

rel_widths <- c(1.25, 4)
rel_heights <- c(1, 1.25, 1)

cat('making aerosol panels...\n')
aer_panel <- plot_grid(
    raw_panel_aer,
    fit_panel_aer,
    half_life_aer,
    align = "v",
    axis = "lr",
    labels = "AUTO",
    rel_heights = rel_heights,
    nrow = 3)

cat('making surface panels...\n')
surface_panel <- plot_grid(
    raw_panel_surface,
    fit_panel_surface,
    half_life_surface,
    align = "v",
    axis = "lr",
    labels = "",
    rel_heights = rel_heights,
    nrow = 3)

cat('making full figure...\n')
full_fig <- plot_grid(
    aer_panel,
    surface_panel,
    label_size = text_size,
    align = "h",
    axis = "lr",
    rel_heights = rel_heights,
    rel_widths = rel_widths,
    ncol = 2)

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          full_fig,
          base_height = panel_fig_height,
          base_asp = 1.3)
warnings()
