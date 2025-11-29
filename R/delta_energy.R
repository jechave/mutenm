## Energy diferences
#' Calculate energy differences between a mutant and wild type
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


## Activation energy


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

