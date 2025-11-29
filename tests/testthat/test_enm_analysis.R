load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "ref_enm_analysis.rda"))

test_that("msfi returns correct values", {
  expect_equal(msfi(wt), ref_msf_site)
})

test_that("msfn returns correct values", {
  expect_equal(msfn(wt), ref_msf_mode)
})
