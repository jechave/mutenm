# Helper function to compare prot objects excluding eigenvectors (umat)
# Eigenvector signs can differ across platforms/LAPACK versions
expect_prot_equal <- function(object, expected) {
  # Compare eigenvalues (these are deterministic)
  expect_equal(object$nma$evalue, expected$nma$evalue)

  # Compare everything except umat
  obj_no_umat <- object
  exp_no_umat <- expected
  obj_no_umat$nma$umat <- NULL
  exp_no_umat$nma$umat <- NULL

  expect_equal(obj_no_umat, exp_no_umat)
}
