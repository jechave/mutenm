load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "mut_qf.rda"))
load(test_path("..", "data", "ref_delta_structure_motion.rda"))

test_that("delta_structure_dr2i returns correct values", {
  expect_equal(delta_structure_dr2i(wt, mut_qf), ref_dr2i)
})

test_that("delta_structure_dr2n returns correct values", {
  expect_equal(delta_structure_dr2n(wt, mut_qf), ref_dr2n)
})

test_that("Dmsfi returns correct values", {
  expect_equal(Dmsfi(wt, mut_qf), ref_dmsfi)
})

test_that("Dmsfn returns correct values", {
  expect_equal(Dmsfn(wt, mut_qf), ref_dmsfn)
})
