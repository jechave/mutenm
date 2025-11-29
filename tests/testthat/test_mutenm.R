load(test_path("..", "data", "pdb_2acy_A.rda"))
load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "mut_lf.rda"))
load(test_path("..", "data", "mut_qf.rda"))


test_that("enm gets wt ok", {
  expect_equal(enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5), wt)
})

test_that("mutenm gets mut_lf", {
  expect_equal(
    mutenm(wt, site_mut = 80, mutation = 1,
           mut_model = "lfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_lf)
})

test_that("mutenm gets mut_qf", {
  expect_equal(
    mutenm(wt,  site_mut = 80, mutation = 1,
           mut_model = "sclfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_qf)
})
