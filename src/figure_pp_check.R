#!/usr/bin/env Rscript

#####################################
## name: figure_pp_check.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## plot posterior predictive checks 
## for all models
##
####################################


script_packages <- c(
    'rstan',     # for stan model processing
    'bayesplot', # for pp check functions
    'cowplot',   # for publication ready ggplot
    'readr',     # for read_csv()
    'magrittr',
    'tidybayes',
    'dplyr')

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
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
mcmc_fit_path <- args[2]
style_path <- args[3]
outpath <- args[4]

## get original data and mcmc output
dat <- read_csv(data_path,
                col_types = cols())

fit <- readRDS(mcmc_fit_path)

## read style files
source(style_path) # shared project style


#####################################
## styling parameters
#####################################
#set.seed(543522) # reproducibility
n_to_plot = 50 # only plot a random subsample, for efficiency


## define plot styling
## pp_color for predictive draw lines, set in plotting_style.R

cat("\n\nPlotting posterior predictive checks...\n")
cat("Ignore coordinate system and colour warnings;",
    "these are expected behavior\n")


chains <- rstan::extract(fit)
n_checks = length(chains$titer_rep_censored[, 1])
sample_checks = sample(1:n_checks, n_to_plot)
print(sample_checks)

print(names(chains))

cat(sprintf("\nPlotting pp checks...\n\n"))
real_titers <- dat$log10_titer
rep_titers <- chains$titer_rep_censored[sample_checks, ]

plot_upper <- pp_check(10^(real_titers),
                 yrep = 10^(rep_titers), ## convert to RML titer
                 fun = ppc_dens_overlay,
                 alpha = 0.2) +
    ggtitle("Predictive checks") +
    scale_color_manual(name="", 
                       labels = c("real data",
                                  "posterior predictive draws"),
                       values = c("black",
                                  pp_color)) +
    scale_x_continuous(
        trans = 'log10',
        labels = trans_format('log10',
                              math_format(10^.x))) +
    coord_cartesian(xlim = c(10^0.5, 10^6.5)) +
    expand_limits(x = 0.5, y = 0) +
    theme_project(base_size = 30)



## save plot to outpath
cat("\nSaving plot...\n\n")
save_plot(outpath,
          plot_upper,
          base_height = 10,
          base_asp = 2)

