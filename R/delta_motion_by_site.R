#' Calculate site-dependent motion difference between two proteins
#'
#' This version works only for wt and mut with no indels
#'
#' @param wt A protein object with \code{enm} defined
#' @param mut A second protein object with \code{enm} defined
#'
#' @return A vector \code{(dmsf_i)} of size \code{nsites}, the MSF change at site i.
#'
#' @export
#'
delta_motion_dmsfi <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dmsf = get_msf_site(mut) - get_msf_site(wt)
  dmsf
}
