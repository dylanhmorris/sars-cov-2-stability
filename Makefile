#####################################
# name: Makefile
# author: Dylan Morris <dhmorris@princeton.edu>
#
# Makefile to generate analyses
# for van Doremalen et al study
# of SARS-CoV-2 decay
# in environment
####################################

#####################################
# Directory structure
####################################

default: all

SRC := src
OUTPUT := out
DATA := dat
RAW := $(DATA)/raw
CLEANED := $(DATA)/cleaned

MCMC_CHAINS := $(OUTPUT)/mcmc_chains
FIGURE_DIR = $(OUTPUT)/figures
TABLE_DIR = $(OUTPUT)/tables

#####################################
# File extensions and the like
####################################
CHAINS_SUFFIX = _chains.Rds


#####################################
# Expected bash settings
#
# Check these vs your local
# machine setup if you are having
# difficulty reproducing the
# analysis
#####################################

CXX := g++-8
MKDIR := @mkdir -p
RM := rm -rf

R_OPTIONS = --vanilla
R_COMMAND := Rscript $(R_OPTIONS)

# R creates a blank Rplots.pdf when run
# from the command line to produce
# figures. This removes that file.
FIG_CLEANUP = @$(RM) Rplots.pdf

#####################################
# Installation / dependencies
#
# Rules for prepping analysis
#####################################

.PHONY: depend

depend:
	$(R_COMMAND) $(SRC)/install_needed_packages.R

#####################################
# data locations
#####################################
TITER_DATAFILE = titer_data.csv
RAW_TITER_DATA = $(RAW)/$(TITER_DATAFILE)
CLEANED_TITER_DATA = $(CLEANED)/$(TITER_DATAFILE)

CLEANED_DATA = $(CLEANED_TITER_DATA)
TITER_DATA = $(CLEANED_TITER_DATA)

#####################################
# code locations
#####################################
CLEANING_SCRIPT = $(SRC)/clean_data.R
FITTING_SCRIPT = $(SRC)/fit_stan_model.R
DIAGNOSTIC_SCRIPT = $(SRC)/chain_diagnostics.R
PLOT_PARAMS = $(SRC)/plotting_style.R

TITER_MODEL_SRC = titer_decay_model.stan

#####################################
# model names and output locations
#####################################
TITER_MODEL_NAME = titer
TITER_CHAINS = $(MCMC_CHAINS)/$(TITER_MODEL_NAME)$(CHAINS_SUFFIX)

MODELS = $(TITER_MODEL_NAME)
CHAIN_PATHS = $(TITER_CHAINS)

CHAIN_DIAGNOSTICS = $(OUTPUT)/chain_diagnostics.csv

MATERIALS = Aerosols Copper Cardboard Steel Plastic
FITSUFFIX = .jpg
INDIVIDUAL_FIT_FIGS = $(addprefix figure_individual_fit_, $(addsuffix $(FITSUFFIX), $(MATERIALS)))
FIGURES = figure_main.pdf $(INDIVIDUAL_FIT_FIGS)
FIGURE_PATHS = $(addprefix $(FIGURE_DIR)/, $(FIGURES))

TABLES = table1_half_life_estimates.docx
TABLE_PATHS = $(addprefix $(TABLE_DIR)/, $(TABLES))
PUBLIC_TABLES = $(TABLES)
PUBLIC_TABLE_PATHS = $(addprefix $(TABLE_DIR)/, $(PUBLIC_TABLES)

#####################################
# Rules
#
# definition of dependency
# tree and specification of
# rules for doing stuff
#####################################

##########################
# rules for data cleaning
##########################

$(CLEANED)/%.csv: $(RAW)/%.csv $(CLEANING_SCRIPT)
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $(CLEANING_SCRIPT) $< false $@

#####################################
# rules for model fitting and post-processing
#####################################

$(TITER_CHAINS): $(FITTING_SCRIPT) $(SRC)/$(TITER_MODEL_SRC) $(TITER_DATA)
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@

$(CHAIN_DIAGNOSTICS): $(DIAGNOSTIC_SCRIPT) $(CHAIN_PATHS)
	$(MKDIR) $(OUTPUT)
	$(R_COMMAND) $^ $@


#####################################
# rules for table generation
#####################################

$(TABLE_DIR)/table1_half_life_estimates.docx: $(SRC)/table_half_life.R $(TITER_DATA) $(TITER_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(TABLE_DIR)
	$(R_COMMAND) $^ $@

#####################################
# rules for figure generation
#####################################

$(FIGURE_DIR)/%.pdf: $(SRC)/%.R $(TITER_DATA) $(TITER_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure_individual_fit_%$(FITSUFFIX): $(SRC)/figure_individual_fits.R $(TITER_DATA) $(TITER_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ true $(@:$(FIGURE_DIR)/figure_individual_fit_%$(FITSUFFIX)=%) $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure_pp_check.pdf: $(SRC)/figure_pp_check.R $(TITER_DATA) $(TITER_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


#####################################
# convenience rules for making
# various quantities
#####################################
.PHONY: data
data: $(CLEANED_DATA)

.PHONY: chains
chains: $(CHAIN_PATHS)

.PHONY: diagnostics
diagnostics: $(CHAIN_DIAGNOSTICS)

.PHONY: figures
figures: $(FIGURE_PATHS)

.PHONY: tables
tables: $(TABLE_PATHS)

.PHONY: echo_figures echo_chains
echo_figures:
	echo $(FIGURE_PATHS)
echo_chains:
	echo $(CHAIN_PATHS)

.PHONY: clean

clean:
	$(RM) $(OUTPUT)
	$(RM) $(CLEANED)

all: depend data chains diagnostics figures tables
