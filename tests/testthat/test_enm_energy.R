load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "ref_enm_energy.rda"))

test_that("v_min returns correct value", {
  expect_equal(v_min(wt), ref_enm_v_min)
})
