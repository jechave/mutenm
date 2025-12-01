# Pair comparison functions
# Functions that compare two prot objects: Df(prot1, prot2)


# Structure ---------------------------------------------------------------

#' Calculate site-dependent structural difference between two proteins
#'
#' This version works only for wt and mut with no indels
#'
#' @param wt A protein object with \code{xyz} defined
#' @param mut A second protein object  with \code{xyz} defined
#'
#' @return A vector \code{(dr2_i)} of size \code{nsites}, the squared displacement at site i.
#'
#' @export
#'
Dr2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}


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
Dr2n <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  drn <- as.vector(crossprod(get_umat(wt), dr))
  stopifnot(length(wt$nma$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}


# Dynamics ----------------------------------------------------------------

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
Dmsfi <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dmsf = msfi(mut) - msfi(wt)
  dmsf
}


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


# Energy ------------------------------------------------------------------

#' Calculate energy differences between two proteins
#'
#' @param wt A protein object
#' @param mut A second protein
#'
#' @return A (scalar) energy difference between mutant and wild type.
#'
#' @name delta_energy
#'
NULL

#' @rdname delta_energy
#'
#' @details `Dv_min` calculates the minimum-energy difference between \code{mut} and \code{wt}
#'
#' @export
#'
Dv_min <- function(wt, mut)
  v_min(mut) - v_min(wt)

#' @rdname delta_energy
#'
#' @details `Dg_ent` calculates the entropic free energy difference between \code{mut} and \code{wt}
#'
#' @export
#'
Dg_ent <- function(wt, mut, beta = beta_boltzmann())
  g_ent(mut, beta) - g_ent(wt, beta)


#' @rdname delta_energy
#'
#' @details `Ddv_act` calculates the energy contribution to the change in activation energy between \code{mut} and \code{wt}
#'
#' @export
#'
Ddv_act <- function(wt, mut, ideal = wt, pdb_site_active = NA) {
  result <- dv_act(mut, ideal, pdb_site_active) - dv_act(wt, ideal, pdb_site_active)
  result
}


#' @rdname delta_energy
#'
#' @details `Ddg_ent_act` calculates the entropy contribution to the change in activation energy between \code{mut} and \code{wt}
#'
#' @export
#'
Ddg_ent_act <- function(wt, mut, ideal = wt, pdb_site_active = NA, beta = beta_boltzmann()) {
  result <- dg_ent_act(mut, ideal, pdb_site_active) - dg_ent_act(wt, ideal, pdb_site_active)
  result
}
