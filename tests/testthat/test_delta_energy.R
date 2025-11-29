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

test_that("delta_energy_dvs returns numeric", {
  result <- delta_energy_dvs(wt, mut_qf)
  expect_type(result, "double")
})

test_that("delta_energy_act_dv returns NA when no active site specified", {
  result <- delta_energy_act_dv(wt, mut_qf)
  expect_true(is.na(result))
})

test_that("delta_energy_act_tds returns NA when no active site specified", {
  result <- delta_energy_act_tds(wt, mut_qf)
  expect_true(is.na(result))
})
