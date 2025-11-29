load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "ref_enm_analysis.rda"))

test_that("get_msf_site returns correct values", {
  expect_equal(get_msf_site(wt), ref_msf_site)
})

test_that("get_msf_mode returns correct values", {
  expect_equal(get_msf_mode(wt), ref_msf_mode)
})
