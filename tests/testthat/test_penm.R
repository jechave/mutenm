library(here)
library(tidyverse)
library(jefuns)
library(bio3d)


load(here("tests/data/pdb_2acy_A.rda"))
load(here("tests/data/wt.rda"))
load(here("tests/data/mut_lf.rda"))
load(here("tests/data/mut_qf.rda"))


test_that("set_enm gets wt ok", {
  expect_equal(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5), wt)
})

test_that("get_mutant_site gets mut_lf", {
  expect_equal(
    get_mutant_site(wt, site_mut = 80, mutation = 1,
                    mut_model = "lfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_lf)
})

test_that("get_mutant_site gets mut_qf", {
  expect_equal(
    get_mutant_site(wt,  site_mut = 80, mutation = 1,
                    mut_model = "sclfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_qf)
})
