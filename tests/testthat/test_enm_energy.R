load(test_path("..", "data", "wt.rda"))
load(test_path("..", "data", "ref_enm_energy.rda"))

test_that("v_min returns correct value", {
  expect_equal(v_min(wt), ref_enm_v_min)
})

test_that("g_ent returns correct value", {
  expect_equal(g_ent(wt, beta = 1.0), ref_enm_g_entropy)
})

test_that("dv_act returns NA when no active site", {
  expect_true(is.na(dv_act(wt, wt, pdb_site_active = NA)))
})

test_that("dv_act returns correct value with active site", {
  expect_equal(dv_act(wt, wt, pdb_site_active), ref_dgact_dv)
})

test_that("dgact_tds returns NA when no active site", {
  expect_true(is.na(dgact_tds(wt, wt, pdb_site_active = NA)))
})

test_that("dgact_tds returns correct value with active site", {
  expect_equal(dgact_tds(wt, wt, pdb_site_active), ref_dgact_tds)
})
