# Refresh data to use in test_mrs_all.R

## ----------------------------------------------------------------------------------------------------------------------
# load libraries
library(tidyverse)
library(bio3d)
library(penm)
library(jefuns)
library(here)
library(tictoc)
library(Matrix)


## ----------------------------------------------------------------------------------------------------------------------
load(here("data/wt.rda"))

mrs_all_output <- mrs_all(wt, nmut = 5, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234)
usethis::use_data(mrs_all_output,  overwrite = TRUE)





