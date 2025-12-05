# Single protein properties
# Functions that take one prot object and return a property


# Structure ---------------------------------------------------------------

#' Contact Number
#'
#' Calculates the Contact Number (CN) of each site - the number of neighbors
#' within the cutoff distance.
#'
#' @inheritParams mutenm-params
#' @returns a vector of size nsites with cn values for each site
#'
#' @family single-protein functions
#'
#' @export
cn <- function(prot) cn_xyz(get_xyz(prot), get_d_max(prot))

#' Weighted Contact Number
#'
#' Calculates the Weighted Contact Number (WCN) of each site - the sum of 1/d²
#' for all pairs.
#'
#' @inheritParams mutenm-params
#' @returns a vector of size nsites with wcn values for each site
#'
#' @family single-protein functions
#'
#' @export
wcn <- function(prot) wcn_xyz(get_xyz(prot))

# Energy ------------------------------------------------------------------

#' Minimum Energy
#'
#' Calculates the minimum (stress) energy of the protein at its current
#' conformation.
#'
#' @inheritParams mutenm-params
#' @return a scalar energy value
#'
#' @family single-protein functions
#'
#' @export
v_min <- function(prot) {
  graph <- get_graph(prot)
  v <- with(graph, {
    v_dij(dij, v0ij, kij, lij)
  })
  v
}

# Internal helpers --------------------------------------------------------

#' Calculate energy of a single spring
#'
#' @noRd
#'
v_dij <- function(dij, v0ij, kij, lij) {
  # Calculates energy of a given conformation (dij).
  sum(v0ij + .5 * kij * (dij - lij) ^ 2)
}


