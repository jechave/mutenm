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

#' Distance to Active Site
#'
#' Calculates the distance from each site to the closest active site residue.
#'
#' @inheritParams mutenm-params
#' @returns a vector of size nsites with dactive values for each site
#'
#' @family single-protein functions
#'
#' @export
dactive <- function(prot, pdb_site_active) {
  xyz <- get_xyz(prot)
  asite <- active_site_indexes(prot, pdb_site_active)
  site_active <- asite$site_active
  result <- dactive.xyz(xyz, site_active)
  result
}


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


#' Activation Energy (Internal)
#'
#' Calculates the internal energy contribution to activation - the cost to
#' deform the active site to the ideal conformation.
#'
#' @inheritParams mutenm-params
#' @return a scalar energy value (NA if pdb_site_active is NA)
#'
#' @family single-protein functions
#'
#' @export
dv_act <- function(prot, ideal, pdb_site_active = NA) {
  if (anyNA(pdb_site_active)) {
    result <- NA
  } else {
    kmat_asite <- kmat_asite(prot, pdb_site_active)
    dxyz_asite <- dxyz_asite(prot, ideal, pdb_site_active)
    result <- .5 * my_quad_form(dxyz_asite, kmat_asite, dxyz_asite)
  }
  result
}

#' Activation Energy (Entropic)
#'
#' Calculates the entropic contribution to activation energy.
#'
#' @inheritParams mutenm-params
#' @return a scalar energy value (NA if pdb_site_active is NA)
#'
#' @family single-protein functions
#'
#' @export
dg_ent_act <- function(prot, ideal, pdb_site_active = NA, beta = beta_boltzmann()) {
  # Calculate entropic contribution to dg_activation
  if (anyNA(pdb_site_active)) {
    result <- NA
  } else {
    # prot, with its active site in a "relaxed" state
    kmat_asite <- kmat_asite(prot, pdb_site_active)
    eig <- eigen(kmat_asite, symmetric = TRUE, only.values = TRUE)
    evalue <- eig$values
    gact_prot <- sum(g_ent_mode(evalue, beta))
    # prot, with its active site in the ideal conformation
    gact_ideal <- 0 # assume in the TS the active-site is rigid

    result <- gact_ideal - gact_prot
  }
  result
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


#' Active site indexes
#'
#' @noRd
#'
active_site_indexes <- function(prot, pdb_site_active) {
  pdb_site_active <- pdb_site_active  # active sites in pdb numbering
  site <- get_site(prot)
  pdb_site <- get_pdb_site(prot)
  site_active <- site[pdb_site %in% pdb_site_active] # active site in protein numbering
  ind_active <- xyz_indices_site(site_active) # indices of xyz coordinates of active sites
  lst(pdb_site_active, site_active, ind_active)
}


#' Caculate effective K matrix of active site
#'
#' @noRd
#'
kmat_asite <- function(prot, pdb_site_active) {
  asite <- active_site_indexes(prot, pdb_site_active)
  cmat <- get_cmat(prot)
  cmat_asite <- cmat[asite$ind_active, asite$ind_active]
  kmat_asite <- solve(cmat_asite)
  kmat_asite
}


#' Activation energy, internal energy term
#'
#' @noRd
#'
dxyz_asite <- function(prot, ideal, pdb_site_active = NA) {
  if (anyNA(pdb_site_active)) {
    result <- NA
  } else {
    asite <- active_site_indexes(prot, pdb_site_active)

    dxyz <- get_xyz(prot) - get_xyz(ideal)
    site_active <- asite$site_active
    dxyz <- my_as_xyz(dxyz)
    result <- as.vector(dxyz[, site_active])
  }
  result
}
