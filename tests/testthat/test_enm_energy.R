load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "ref_enm_energy.rda"))

test_that("v_min returns correct value", {
  expect_equal(v_min(wt), ref_enm_v_min)
})

test_that("enm_g_entropy returns correct value", {
  expect_equal(enm_g_entropy(wt, beta = 1.0), ref_enm_g_entropy)
})

test_that("dgact_dv returns NA when no active site", {
  expect_true(is.na(dgact_dv(wt, wt, pdb_site_active = NA)))
})

test_that("dgact_dv returns correct value with active site", {
  expect_equal(dgact_dv(wt, wt, pdb_site_active), ref_dgact_dv)
})

test_that("dgact_tds returns NA when no active site", {
  expect_true(is.na(dgact_tds(wt, wt, pdb_site_active = NA)))
})

test_that("dgact_tds returns correct value with active site", {
  expect_equal(dgact_tds(wt, wt, pdb_site_active), ref_dgact_tds)
})
