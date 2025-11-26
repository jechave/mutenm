## ----------------------------------------------------------------------------------------------------------------------
# Test Mutation Response Scanning functions

# load libraries
library(tidyverse)
library(bio3d)
library(penm)
library(jefuns)
library(here)
library(tictoc)
library(Matrix)

load(here("tests/data/wt.rda"))

load(here("tests/data/mrs_all_output.rda"))
mrs_all_output_test <- mrs_all(wt, nmut = 5, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234)
test_that("mrs with option mean_max is ok", {
  expect_equal(mrs_all_output_test, mrs_all_output)
  })

load(here("tests/data/smrs_all_output.rda"))
smrs_all_output_test <- smrs_all(wt, nmut = 5, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234)
test_that("smrs with option mean_max is ok", {
  expect_equal(smrs_all_output_test, smrs_all_output)
})



