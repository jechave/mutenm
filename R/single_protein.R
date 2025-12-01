# Single protein properties
# Functions that take one prot object and return a property


# Structure ---------------------------------------------------------------

#' Calculate CN site-dependent profile
#'
#' Calculates the Contact Number (CN) of each site
#'
#' @param prot is a protein object obtained using enm()
#' @returns a vector of size nsites with cn values for each site
#'
#' @export
#'
#' @family site profiles
#'
#'
cn <- function(prot) cn_xyz(get_xyz(prot), get_d_max(prot))

#' Calculate WCN site-dependent profile
#'
#' Calculates the Weighted Contact Number (WCN) of each site
#'
#' @param prot is a protein object obtained using enm()
#' @returns a vector of size nsites with wcn values for each site
#'
#' @export
#'
#' @family site profiles
#'
wcn <- function(prot) wcn_xyz(get_xyz(prot))


#' Calculate distance to active site
#'
#' Calculates the distance from each site to the closest active residue
#'
#' @param prot is a protein object obtained using enm()
#' @param pdb_site_active is a vector of pdb resno of active residues.
#'
#' @returns a vector of size nsites with dactive values for each site
#'
#' @export
#'
#' @family site profiles
#'
dactive <- function(prot, pdb_site_active) {
  xyz <- get_xyz(prot)
  asite <- active_site_indexes(prot, pdb_site_active)
  site_active <- asite$site_active
  result <- dactive.xyz(xyz, site_active)
  result
}


# Dynamics ----------------------------------------------------------------

#' Calculate MSF site-dependent profile
#'
#' Calculates the mean-square-fluctuation of each site
#'
#' @param prot is a protein object obtained using enm()
#' @returns a vector of size nsites with msf values for each site
#'
#' @export
#'
#'
#' @family site profiles
#'
msfi <- function(prot) {
  diag(get_reduced_cmat(prot))
}


#' Calculate MSF mode-dependent profile
#'
#' Calculates the mean-square-fluctuation in the direction of each normal mode
#'
#' @param prot is a protein object obtained using enm()
#' @returns a vector of size nsites with msf values for each mode
#'
#' @export
#'
#'
#' @family mode profiles
#'
msfn <-  function(prot) 1 / get_evalue(prot)


# Energy ------------------------------------------------------------------

#' Calculate minimum energy of a given prot object
#' @param prot is a prot object, wit a component graph tibble
#' where v0ij, kij, lij and the dij for the minimum conformation are found.
#' @return a scalar: the energy at the minimum-enegy conformation
#'
#' @export
#' @family enm_energy
#'
v_min <- function(prot) {
  graph <- get_graph(prot)
  v <- with(graph, {
    v_dij(dij, v0ij, kij, lij)
  })
  v
}


#' Calculate entropic total free energy of prot object
#'
#' @param prot is a prot object with known eigenvalues of enm model
#' @param beta is 1 / kT
#' @family enm_energy
#' @export
#'
g_ent <- function(prot, beta) {
  # Calculate T*S from the energy spectrum
  energy <- get_evalue(prot)
  sum(g_ent_mode(energy, beta))
}


#' Activation free energy, internal energy contribution
#'
#' @export
#' @family enm_energy
#'
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


#' Activation free energy, entropic contribution
#'
#' @export
#' @family enm_energy
#'
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
