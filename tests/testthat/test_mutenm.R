load(test_path("..", "data", "pdb_2acy_A.rda"))
load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "mut_lf.rda"))
load(test_path("..", "data", "mut_qf.rda"))


test_that("enm gets wt ok", {
  result <- enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5)
  expect_prot_equal(result, wt)
})

test_that("mutenm gets mut_lf", {
  result <- mutenm(wt, site_mut = 80, mutation = 1,
                   mut_model = "lfenm", mut_sd_min = 1, mut_dl_sigma = 0.3)
  expect_prot_equal(result, mut_lf)
})

test_that("mutenm gets mut_qf", {
  result <- mutenm(wt, site_mut = 80, mutation = 1,
                   mut_model = "sclfenm", mut_sd_min = 1, mut_dl_sigma = 0.3)
  expect_prot_equal(result, mut_qf)
})
