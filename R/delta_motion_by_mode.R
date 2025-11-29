#' Calculate mode-dependent motion difference between two proteins
#'
#' This version works only for wt and mut with no indels.
#' Calculates fluctuations along normal modes of wt, independent of mode assignment issues.
#'
#' @param wt A protein object with \code{enm} defined
#' @param mut A second protein object with \code{enm} defined
#'
#' @return A vector \code{(dmsf_n)} of size \code{nmodes}, the MSF change in mode n.
#'
#' @export
#'
Dmsfn <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  msf_wt <- msfn(wt)
  msf_mut <- diag(t(get_umat(wt)) %*% (get_cmat(mut) %*% get_umat(wt)))
  dmsf = msf_mut - msf_wt
  dmsf
}
