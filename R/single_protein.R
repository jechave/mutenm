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

# Dynamics ----------------------------------------------------------------

#' Mean-Square Fluctuation per Site
#'
#' Calculates the mean-square fluctuation of each site (diagonal of the
#' reduced covariance matrix).
#'
#' @inheritParams mutenm-params
#' @returns a vector of size nsites with msf values for each site
#'
#' @family single-protein functions
#'
#' @export
msfi <- function(prot) {
  diag(get_reduced_cmat(prot))
}

#' Mean-Square Fluctuation per Mode
#'
#' Calculates the mean-square fluctuation contribution from each normal mode
#' (1/eigenvalue).
#'
#' @inheritParams mutenm-params
#' @returns a vector of size nmodes with msf values for each mode
#'
#' @family single-protein functions
#'
#' @export
msfn <- function(prot) 1 / get_evalue(prot)


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

#' Entropic Free Energy
#'
#' Calculates the entropic free energy contribution from vibrational modes.
#'
#' @inheritParams mutenm-params
#' @return a scalar energy value
#'
#' @family single-protein functions
#'
#' @export
g_ent <- function(prot, beta = beta_boltzmann()) {
  # Calculate T*S from the energy spectrum
  energy <- get_evalue(prot)
  sum(g_ent_mode(energy, beta))
}


# Internal helpers --------------------------------------------------------

#' Calculate reduced covariance matrix
#'
#' @noRd
#'
get_reduced_cmat <- function(prot) {
  get_cmat(prot) %>%
    reduce_matrix()
}


#' Calculate energy of a single spring
#'
#' @noRd
#'
v_dij <- function(dij, v0ij, kij, lij) {
  # Calculates energy of a given conformation (dij).
  sum(v0ij + .5 * kij * (dij - lij) ^ 2)
}


#' Entropic contribution of a single mode
#'
#' @noRd
#'
g_ent_mode <- function(energy, beta) {
  # returns vector of entropic terms given vector of mode energies
  g_entropy_mode <- 1 / (2 * beta) * log((beta * energy) / (2 * pi))
  g_entropy_mode
}


