#!/usr/bin/env Rscript

#####################################
## name: chain_diagnostics.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## read in mcmc chains
## and output diagnostics
##
####################################

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(bayesplot)) # for rhat, neff_ratio

## read command line args
args <- commandArgs(trailingOnly=TRUE)
n_models <- length(args) - 1
mcmc_fit_paths <- args[1:n_models]
outpath <- args[n_models + 1]

rhat_max <- function(chains){
    rhat_vals <- rhat(fit)

    if (any(is.nan(rhat_vals))) {
        warning(paste0(
            "Warning: removing NaN rhat values. This may be benign ",
            "if, for example, there is a parameter that is not sampled ",
            "because it is fully constrained"))
        rhat_vals <- rhat_vals[!is.nan(rhat_vals)]
    }

    rhat_max_parm <- names(which.max(rhat_vals))
    rhat_max_val <- max(rhat_vals)

    return (list(rhat_max_parm, rhat_max_val))

}

neff_min <- function(chains){

    neff_ratios <- neff_ratio(fit)
    
    if (any(is.nan(neff_ratios))) {
        warning(paste0(
            "Warning: removing NaN neff_ratio ",
            "values. This may be benign ",
            "if, for example, there is a parameter ",
            "that is not sampled ",
            "because it is fully constrained"))
        neff_ratios <- neff_ratios[!is.nan(neff_ratios)]
    }

    neff_min_parm <- names(which.min(neff_ratios))
    neff_min_val <- min(neff_ratios)
    return (list(neff_min_parm, neff_min_val))

}

models <- rep('', n_models)
rhat_max_parms <- rep('', n_models)
rhat_max_vals <- rep(NA, n_models)
neff_min_parms <- rep('', n_models)
neff_min_vals <- rep(NA, n_models)

## get mcmc chains
for (k in 1:length(mcmc_fit_paths)) {

    path <- mcmc_fit_paths[k]
    model_name <- path
    
    cat(sprintf("\nDiagnostics for model %s...\n",
                model_name))
    fit <- readRDS(path)
    chains <- extract(fit)
    rhat_maxes <- rhat_max(chains)
    neff_mins <- neff_min(chains)

    models[k] <- model_name
    rhat_max_parms[k] <- rhat_maxes[[1]]
    rhat_max_vals[k] <- rhat_maxes[[2]]
    neff_min_parms[k] <- neff_mins[[1]]
    neff_min_vals[k] <- neff_mins[[2]]
}

diagnostic_data <- data.frame(
    models = models,
    rhat_max_parm = rhat_max_parms,
    rhat_max_val = rhat_max_vals,
    neff_min_parm = neff_min_parms,
    neff_min_val = neff_min_vals)

write.csv(diagnostic_data, outpath)
