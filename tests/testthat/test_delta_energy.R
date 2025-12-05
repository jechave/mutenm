load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "mut_qf.rda"))

test_that("Dv_min returns numeric", {
  result <- Dv_min(wt, mut_qf)
  expect_type(result, "double")
})

test_that("Dg_ent returns non-zero for sclfenm mutant", {
  result <- Dg_ent(wt, mut_qf)
  expect_type(result, "double")
  expect_false(result == 0)
})
