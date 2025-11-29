library(here)

load(here("tests/data/wt.rda"))

# 2acy active site residues (PDB numbering): Arg23, Asn41
pdb_site_active <- c(23, 41)

test_that("enm_v_min returns numeric scalar", {
  result <- enm_v_min(wt)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("enm_g_entropy returns numeric scalar", {
  result <- enm_g_entropy(wt, beta = 1.0)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("dgact_dv returns NA when no active site specified", {
  result <- dgact_dv(wt, wt, pdb_site_active = NA)
  expect_true(is.na(result))
})

test_that("dgact_dv returns numeric when active site specified", {
  result <- dgact_dv(wt, wt, pdb_site_active = pdb_site_active)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("dgact_tds returns NA when no active site specified", {
  result <- dgact_tds(wt, wt, pdb_site_active = NA)
  expect_true(is.na(result))
})

test_that("dgact_tds returns numeric when active site specified", {
  result <- dgact_tds(wt, wt, pdb_site_active = pdb_site_active)
  expect_type(result, "double")
  expect_length(result, 1)
})
