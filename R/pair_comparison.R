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




