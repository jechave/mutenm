# Pair comparison functions
# Functions that compare two prot objects: Df(prot1, prot2)


# Structure ---------------------------------------------------------------

#' Squared Displacement per Site
#'
#' Calculates the squared displacement at each site between wild-type and mutant.
#'
#' @inheritParams mutenm-params
#' @return a vector of size nsites with squared displacements
#'
#' @family mutation-effect functions
#'
#' @export
Dr2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}

#' Squared Displacement per Mode
#'
#' Calculates the squared displacement projected onto each normal mode of the
#' wild-type.
#'
#' @inheritParams mutenm-params
#' @return a vector of size nmodes with squared displacements
#'
#' @family mutation-effect functions
#'
#' @export
Dr2n <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  drn <- as.vector(crossprod(get_umat(wt), dr))
  stopifnot(length(wt$nma$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}


# Dynamics ----------------------------------------------------------------

#' Change in Mean-Square Fluctuation per Site
#'
#' Calculates the change in mean-square fluctuation at each site between
#' wild-type and mutant.
#'
#' Note: Requires the `sclfenm` mutation model to be meaningful.
#'
#' @inheritParams mutenm-params
#' @return a vector of size nsites
#'
#' @family mutation-effect functions
#'
#' @export
Dmsfi <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dmsf = msfi(mut) - msfi(wt)
  dmsf
}

#' Change in Mean-Square Fluctuation per Mode
#'
#' Calculates the change in MSF projected onto each normal mode of the wild-type.
#'
#' Note: Requires the `sclfenm` mutation model to be meaningful.
#'
#' @inheritParams mutenm-params
#' @return a vector of size nmodes
#'
#' @family mutation-effect functions
#'
#' @export
Dmsfn <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  msf_wt <- msfn(wt)
  msf_mut <- diag(t(get_umat(wt)) %*% (get_cmat(mut) %*% get_umat(wt)))
  dmsf = msf_mut - msf_wt
  dmsf
}


# Energy ------------------------------------------------------------------

#' Change in Minimum Energy
#'
#' Calculates the change in minimum (stress) energy between wild-type and mutant.
#'
#' @inheritParams mutenm-params
#' @return a scalar energy difference
#'
#' @family mutation-effect functions
#'
#' @export
Dv_min <- function(wt, mut)
  v_min(mut) - v_min(wt)

#' Change in Entropic Free Energy
#'
#' Calculates the change in entropic free energy between wild-type and mutant.
#'
#' @inheritParams mutenm-params
#' @return a scalar energy difference
#'
#' @family mutation-effect functions
#'
#' @export
Dg_ent <- function(wt, mut, beta = beta_boltzmann())
  g_ent(mut, beta) - g_ent(wt, beta)


