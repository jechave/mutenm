# Calculate various protein properties


# site profiles ----------------------------------------------------

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



#' Calculate MLMS site-dependent profile
#'
#' Calculates the Mean Local Mutational Stress (MLMS) profile using graph of prot object
#'
#' @param prot is a protein object obtained using enm()
#' @param sdij_cut An integer cutoff of sequence distance to include in calculation
#' @returns the profile of mean-local-mutational-stress (mlms) values
#'
#' @export
#'
#'
#' @family site profiles
#'
get_mlms <- function(prot, sdij_cut = 2) {
  g1 <- get_graph(prot)
  g2 <- g1 %>%
    select(edge, j, i, v0ij, sdij, lij, kij, dij)
  names(g2) <- names(g1)
  g <- rbind(g1, g2)

  g <- g %>%
    filter(sdij >= sdij_cut) %>%
    group_by(i) %>%
    summarise(mlms = sum(kij))  %>%
    select(mlms)

  as.vector(g$mlms)

}




# Get mode profiles -------------------------------------------------------

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




# internal helper -----------------------------------------------

#' Calculate reduced covariance matrix
#'
#' @noRd
#'
get_reduced_cmat <- function(prot) {
  get_cmat(prot) %>%
    reduce_matrix()
}




