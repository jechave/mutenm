#' Calculate mode-dependent structural difference between two proteins
#'
#' This version works only for wt and mut with no indels
#'
#' @param wt A protein object with \code{xyz} and \code{enm} defined
#' @param mut A second protein object with \code{xyz} defined
#'
#' @return A vector \code{(dr2_n)} of size \code{nmodes}, the squared displacement in mode n.
#'
#' @export
#'
delta_structure_dr2n <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  drn <- as.vector(crossprod(get_umat(wt), dr))
  stopifnot(length(wt$nma$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}
