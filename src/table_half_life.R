#!/usr/bin/env Rscript

###################################
## analyze the results of a fit
## stan model of titer decay
####################################

script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'magrittr',   # for pipe operator %>%
    'dplyr',      # for filter()
    'tidyr',      # for spread()
    'tidybayes',  # for spread_draws()
    'huxtable',   # for programmatic results output
    'officer',    # for outputing tables to docx
    'flextable'   # for table output
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


## read data
chains <- readRDS(results_path)

source(plotting_params_path) # shared project style

dat <- read_data_for_plotting(data_path)
sorting_virus <- "HCoV-19"

half_life <- chains %>%
    spread_draws(decay_rate[trial_unique_id]) %>%
    mutate(half_life = log10(2) / decay_rate) %>%
    ungroup()

half_life <- dat %>%
    distinct(trial_unique_id, .keep_all = TRUE) %>%
    select(trial_unique_id,
           virus,
           material) %>%
    inner_join(half_life, by = "trial_unique_id") %>%
    group_by(trial_unique_id)

ordering <- half_life %>%
    filter(virus == sorting_virus,
           material != "Aerosols") %>%
    group_by(material) %>%
    summarise(med = median(-decay_rate)) %>%
    arrange(med)

full_summary <- half_life %>%
    ungroup() %>%
    pivot_wider(names_from = virus,
                values_from = half_life,
                id_cols = c(".draw", "material")) %>%
    rename(HCoV19 = `HCoV-19`, SARSCoV1 = `SARS-CoV-1`) %>%
    select(.draw, material, HCoV19, SARSCoV1) %>%
    mutate(hl_diff = HCoV19 - SARSCoV1) %>%
    group_by(material) %>%
    summarise(
        hcov_median = quantile(HCoV19, 0.5),
        hcov_q025 = quantile(HCoV19, 0.025),
        hcov_q975 = quantile(HCoV19, 0.975),
        sars_median = quantile(SARSCoV1, 0.5),
        sars_q025 = quantile(SARSCoV1, 0.025),
        sars_q975 = quantile(SARSCoV1, 0.975),
        diff_median = quantile(hl_diff, 0.5),
        diff_q025 = quantile(hl_diff, 0.025),
        diff_q975 = quantile(hl_diff, 0.975)) %>%
    ungroup() %>%
    mutate(material = factor(material,
                             levels = c("Aerosols", ordering$material))) %>%
    arrange(material)


#################################
## Construct and format table
#################################

output_table <- as_hux(full_summary)
header_row <- c("Material", rep(c("median", "2.5%", "97.5%"), 3))
unit_row <- c("",
              "half-life (hrs)", "", "",
              "half-life (hrs)", "", "",
              "difference (hrs)", "", "")
virus_row <- c("",
               hcov_name, "", "",
               sars_name, "", "",
               "HCoV-19 - SARS-CoV-1", "", "")



output_table <- output_table %>%
    insert_row(header_row,
               after = 0) %>%
    insert_row(unit_row,
               after = 0) %>%
    insert_row(virus_row,
               after = 0) %>%
    merge_cells(1, 2:4) %>%
    merge_cells(1, 5:7) %>%
    merge_cells(1, 8:10) %>%
    merge_cells(2, 2:4) %>%
    merge_cells(2, 5:7) %>%
    merge_cells(2, 8:10)

    


## style the table a bit

bottom_border(output_table)[2, ] <- 1



cat('Saving table to ', outpath, '...\n')
quick_docx(output_table,
           file = outpath)
warnings()
