load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "mut_qf.rda"))
load(test_path("..", "data", "ref_delta_structure_motion.rda"))

test_that("Dr2i returns correct values", {
  expect_equal(Dr2i(wt, mut_qf), ref_dr2i)
})
