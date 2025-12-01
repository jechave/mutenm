# Test mrs() function

test_that("mrs produces expected output", {
  # Load reference data
  load(test_path("..", "data", "mrs_all_output.rda"))

  # Load wt
  load(test_path("..", "data", "prot_2acy_A_ming_wall_ca.rda"))
  wt <- prot_2acy_A_ming_wall_ca

  # Generate mrs output with same parameters
  result <- mrs(wt, nmut = 2, mut_model = "sclfenm", seed = 42)

  # Compare site responses
  expect_equal(
    arrange(result$Dr2i, i, j)$Dr2i,
    arrange(mrs_all_output$Dr2i, i, j)$Dr2i,
    tolerance = 1e-10
  )
  expect_equal(
    arrange(result$Dmsfi, i, j)$Dmsfi,
    arrange(mrs_all_output$Dmsfi, i, j)$Dmsfi,
    tolerance = 1e-10
  )

  # Compare mode responses
  expect_equal(
    arrange(result$Dr2n, n, j)$Dr2n,
    arrange(mrs_all_output$Dr2n, n, j)$Dr2n,
    tolerance = 1e-10
  )
  expect_equal(
    arrange(result$Dmsfn, n, j)$Dmsfn,
    arrange(mrs_all_output$Dmsfn, n, j)$Dmsfn,
    tolerance = 1e-10
  )

  # Compare scalar responses
  expect_equal(result$Dv_min$Dv_min, mrs_all_output$Dv_min$Dv_min, tolerance = 1e-10)
  expect_equal(result$Dg_ent$Dg_ent, mrs_all_output$Dg_ent$Dg_ent, tolerance = 1e-10)
  expect_equal(result$Ddv_act$Ddv_act, mrs_all_output$Ddv_act$Ddv_act, tolerance = 1e-10)
  expect_equal(result$Ddg_ent_act$Ddg_ent_act, mrs_all_output$Ddg_ent_act$Ddg_ent_act, tolerance = 1e-10)
})

test_that("mrs with lfenm produces warning", {
  load(test_path("..", "data", "prot_2acy_A_ming_wall_ca.rda"))
  wt <- prot_2acy_A_ming_wall_ca

  expect_warning(
    mrs(wt, nmut = 1, mut_model = "lfenm", seed = 42),
    "lfenm: dynamics unchanged"
  )
})
