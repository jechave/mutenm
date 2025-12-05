# Pair comparison functions
# Functions that compare two prot objects: Df(prot1, prot2)


# Structure ---------------------------------------------------------------

#' Squared Displacement per Site
#'
#' Calculates the squared displacement at each site between wild-type and mutant.
#'
#' @inheritParams mutenm-params
#' @return a vector of size nsites with squared displacements
#'
#' @family mutation-effect functions
#'
#' @export
Dr2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}





