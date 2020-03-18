#!/usr/bin/env Rscript

#####################################
## name: clean_data.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## process raw data and save cleaned
## data for use in model fitting
##
####################################


script_packages <- c(
    'readr',     # csv read-in
    'dplyr'      # grouping and combining
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}

## functions

round_nearest <- function(x, stepsize){
    return( round(x / stepsize) * stepsize )
}

## specify needed cleaning
clean_titer_data <- function(dat,
                             use_mers = TRUE){

    ## use current HCoV-19 nomenclature
    ## rather than SARS-CoV-2
    dat <- dat %>%
        mutate(virus = replace(virus, virus == "SARS-CoV-2",
                               "HCoV-19"))

    dat$n_wells <- ifelse((grepl("MERS", dat$virus) &
                             dat$material == "Copper"),
                            3,
                          4)

    dat$detection_limit_log10_titer <- ifelse(
        dat$material == "Copper",
        1.5,
        0.5)
    
    log_titer_stepsize = 1 / dat$n_wells
    
    if (!"log10_titer" %in% names(dat)){

        ## undetectable titers are always
        ## reported as 0.5 log10 units
        ## (3.16... unlogged), but others
        ## are in intervals of either 1/3s
        ## or 1/4s
        dat$log10_titer = ifelse(
            dat$titer > 3.17,
            0.5 + round_nearest(
                      log10(dat$titer) - 0.5,
                      log_titer_stepsize),
            0.5)
    }

    if (!"strain" %in% names(dat)) {
        dat$strain <- dat$virus

    }

    if (!use_mers) {
        dat <- dat %>%
            filter(virus != "MERS-CoV")
    }

    ## for multi-virus fitting,
    ## include both the current
    ## MERS-CoV EMC/12 strain experiments
    ## and the prior ones (but as different trials)
    if ("xp" %in% names(dat)) {
        dat <- dat %>%
            filter(xp == "new" |
                   strain == "MERS-CoV_EMC/12") %>%
            mutate(trial_unique_name = paste0(
                       material,
                       " ",
                       temperature,
                       "°C ",
                       humidity,
                       "% RH ",
                       sub("MERS_", "", xp)),
                   grouping_name = paste0(
                       material,
                       " ",
                       temperature,
                       "°C ",
                       humidity,
                       "% RH"))
        
        dat <- dat %>%
            arrange(virus,
                    strain,
                    temperature,
                    humidity,
                    material,
                    xp) %>%
            group_by(virus,
                     strain,
                     temperature,
                     humidity,
                     material,
                     xp) %>%
            mutate(trial_unique_id = group_indices())
        
        dat <- dat %>%
            arrange(virus,
                    material,
                    temperature,
                    humidity,
                    replicate,
                    time,
                    xp) %>%
            group_by(virus,
                     material,
                     temperature,
                     humidity,
                     replicate,
                     time,
                     xp) %>%
            mutate(titer_id = group_indices())


    } else {
        dat <- dat %>%
            mutate(trial_unique_name = paste0(
                       virus,
                       " ",
                       material,
                       " ",
                       temperature,
                       "°C ",
                       humidity,
                       "% RH"),
                   grouping_name = paste0(
                       material,
                       " ",
                       temperature,
                       "°C ",
                       humidity,
                       "% RH"))

        dat <- dat %>%
            arrange(virus,
                    strain,
                    temperature,
                    humidity,
                    material) %>%
            group_by(virus,
                     strain,
                     temperature,
                     humidity,
                     material) %>%
            mutate(trial_unique_id =
                       group_indices())

        dat <- dat %>%
            arrange(virus,
                    material,
                    temperature,
                    humidity,
                    replicate,
                    time) %>%
            group_by(virus,
                     material,
                     temperature,
                     humidity,
                     replicate,
                     time) %>%
            mutate(titer_id = group_indices())
    }

    dat <- dat %>% arrange(titer_id)
    return(dat)
}


clean_well_data <- function(dat, use_mers = FALSE) {
    ## use current HCoV-19 nomenclature
    ## rather than SARS-CoV-2
    dat <- dat %>%
        mutate(virus = replace(virus, virus == "SARS-CoV-2",
                               "HCoV-19"))

    ## filter out MERS-CoV if not to be used
    if (!use_mers) {
        dat <- dat %>%
            filter(virus != "MERS-CoV")
    }

    ## for now, do not estimate dose deposited
    ## loss rate -- maybe add if everything
    ## else works
    dat <- dat %>% filter(time >= 0)
    

    ## add naming and trial ids
    dat <- dat %>%
        mutate(trial_unique_name = paste0(
                   virus,
                   " ",
                   material,
                   " ",
                   temperature,
                   "°C ",
                   humidity,
                   "% RH"),
               grouping_name = paste0(
                   material,
                   " ",
                   temperature,
                   "°C ",
                   humidity,
                   "% RH"),
               )

    dat <- dat %>%
        group_by(virus,
                 temperature,
                 humidity,
                 material) %>%
        mutate(trial_unique_id = group_indices())

    dat <- dat %>%
        group_by(virus,
                 material,
                 temperature,
                 humidity,
                 replicate,
                 time) %>%
        mutate(titer_id = group_indices())

    dat <- dat %>% arrange(titer_id)
    return (dat)
}


## load and clean data
args <- commandArgs(trailingOnly=TRUE)
raw_data_path <- args[1]
use_mers <- as.logical(args[2])

delim <- ";"

dat <- read_delim(raw_data_path,
                  delim = delim,
                  col_types = cols())

if ( any(grepl("well", names(dat))) ){
    cleaned <- clean_well_data(
        dat,
        use_mers = use_mers)
} else {

    cleaned <- clean_titer_data(dat,
                                use_mers = use_mers)
}

print(cleaned)

write_csv(cleaned,
          args[3])

warnings()
