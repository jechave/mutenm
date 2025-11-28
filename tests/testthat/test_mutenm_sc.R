library(here)
library(tidyverse)
library(jefuns)
library(bio3d)


load(here("tests/data/pdb_2acy_A.rda"))
load(here("tests/data/wt_sc.rda"))
load(here("tests/data/mut_sc_lf.rda"))
load(here("tests/data/mut_sc_qf.rda"))


test_that("enm gets wt_sc ok", {
  expect_equal(enm(pdb_2acy_A, node = "sc", model = "ming_wall", d_max = 10.5),
               wt_sc)
})

test_that("mutenm gets mut_sc_lf", {
  expect_equal(
    mutenm(wt_sc, site_mut = 80, mutation = 1,
           mut_model = "lfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_sc_lf)
})

test_that("mutenm gets mut_sc_qf", {
  expect_equal(
    mutenm(wt_sc,  site_mut = 80, mutation = 1,
           mut_model = "sclfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_sc_qf)
})
