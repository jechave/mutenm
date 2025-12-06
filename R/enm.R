
# Suppress R CMD check notes for tidyverse non-standard evaluation
utils::globalVariables(c("dij", "sdij", "edge", "v0ij", "lij", "kij", "i", "j"))

# Create and set prot object ----------------------------------------------


#' Create an Elastic Network Model
#'
#' @description
#' Creates a `prot` object containing information on ENM structure, parameters,
#' and normal modes.
#'
#' @details
#' An Elastic Network Model (ENM) represents a protein as a network of nodes
#' connected by harmonic springs. The network's vibrational modes capture the
#' protein's intrinsic flexibility and are used to predict mutation effects.
#'
#' \subsection{Node types}{
#' One node per residue, placed at:
#' \itemize{
#'   \item \code{"ca"}: Calpha coordinates
#'   \item \code{"cb"}: Cbeta coordinates
#'   \item \code{"sc"}: Side-chain centroid coordinates
#' }
#' }
#'
#' \subsection{Models}{
#' \itemize{
#'   \item \code{"anm"}: Anisotropic Network Model (Atilgan et al., 2001) - uniform spring constants (k=1 for d <= d_max)
#'   \item \code{"ming_wall"}: Ming-Wall model (Ming & Wall, 2005) - uniform k with strengthened backbone (k_backbone = 42k)
#'   \item \code{"pfanm"}: Parameter-free ANM (Yang et al., 2009) - distance-dependent springs (k ~ 1/d^2)
#'   \item \code{"hnm"}: Hinsen model (Hinsen, 1998) - piecewise function: linear for d <= 4A, decays as 1/d^6 otherwise
#'   \item \code{"hnm0"}: Hinsen exponential model (Hinsen et al., 2000) - Gaussian decay exp(-(d/c)^2)
#'   \item \code{"reach"}: REACH model (Moritsugu & Smith, 2007) - sequence-distance dependent with exponential decay
#' }
#' }
#'
#' \subsection{The prot object}{
#' The returned \code{prot} object is a list of class "prot" containing:
#' \describe{
#'   \item{\code{param}}{List of ENM parameters: \code{node}, \code{model}, \code{d_max}}
#'   \item{\code{nodes}}{List with node information:
#'     \itemize{
#'       \item \code{nsites}: number of residues/nodes
#'       \item \code{site}: sequential site indices (1 to nsites)
#'       \item \code{pdb_site}: residue numbers from PDB file
#'       \item \code{bfactor}: crystallographic B-factors
#'       \item \code{xyz}: 3 x nsites matrix of node coordinates
#'     }
#'   }
#'   \item{\code{graph}}{Tibble with network edges (one row per spring):
#'     \itemize{
#'       \item \code{edge}: edge label "i-j"
#'       \item \code{i}, \code{j}: connected node indices
#'       \item \code{sdij}: sequence distance |j-i|
#'       \item \code{dij}: spatial distance (Angstroms)
#'       \item \code{lij}: equilibrium spring length (= dij for wild-type)
#'       \item \code{kij}: spring constant
#'       \item \code{v0ij}: spring potential energy offset (0 for wild-type)
#'     }
#'   }
#'   \item{\code{eij}}{Matrix (n_edges x 3) of unit vectors along each spring}
#'   \item{\code{kmat}}{Hessian matrix (3N x 3N) of second derivatives of potential energy}
#'   \item{\code{nma}}{Normal mode analysis results:
#'     \itemize{
#'       \item \code{mode}: mode indices
#'       \item \code{evalue}: eigenvalues (squared frequencies)
#'       \item \code{umat}: eigenvectors (normal modes), 3N x n_modes matrix
#'       \item \code{cmat}: covariance matrix (3N x 3N), used for linear response
#'     }
#'   }
#' }
#' }
#'
#' @param pdb   pdb object obtained using bio3d::read.pdb
#' @param node  node type: "ca", "cb", or "sc" (default: "ca")
#' @param model ENM model: "anm", "ming_wall", "pfanm", "hnm", "hnm0", or "reach" (default: "anm")
#' @param d_max distance cutoff in Angstroms (default: 10.5)
#'
#' @returns an object of class `prot`, which is a list `lst(param, node, graph, eij, kmat, nma)`
#'
#' @references
#' Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin O, Bahar I (2001).
#' Anisotropy of fluctuation dynamics of proteins with an elastic network model.
#' \emph{Biophysical Journal}, 80:505-515. \doi{10.1016/S0006-3495(01)76033-X}
#'
#' Hinsen K (1998). Analysis of domain motions by approximate normal mode calculations.
#' \emph{Proteins}, 33:417-429. \doi{10.1002/(SICI)1097-0134(19981115)33:3<417::AID-PROT10>3.0.CO;2-8}
#'
#' Hinsen K, Petrescu AJ, Dellerue S, Bellissent-Funel MC, Kneller GR (2000).
#' Harmonicity in slow protein dynamics.
#' \emph{Chemical Physics}, 261:25-37. \doi{10.1016/S0301-0104(00)00222-6}
#'
#' Ming D, Wall ME (2005). Allostery in a coarse-grained model of protein dynamics.
#' \emph{Physical Review Letters}, 95:198103. \doi{10.1103/PhysRevLett.95.198103}
#'
#' Moritsugu K, Smith JC (2007). Coarse-grained biomolecular simulation with REACH:
#' Realistic extension algorithm via covariance Hessian.
#' \emph{Biophysical Journal}, 93:3460-3469. \doi{10.1529/biophysj.107.111898}
#'
#' Yang L, Song G, Jernigan RL (2009). Protein elastic network models and the ranges of cooperativity.
#' \emph{PNAS}, 106:12347-12352. \doi{10.1073/pnas.0902159106}
#'
#' @family core functions
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5)
#' enm(pdb, node = "sc", model = "anm", d_max = 12.5)
#' enm(pdb, node = "cb", model = "anm", d_max = 12.0)
#' }
enm <- function(pdb, node = "ca", model = "anm", d_max = 10.5) {

  prot <- create_enm() %>%
    set_enm_param(node = node, model = model, d_max = d_max) %>%
    set_enm_nodes(pdb = pdb) %>%
    set_enm_graph() %>%
    set_enm_eij() %>%
    set_enm_kmat() %>%
    set_enm_nma()

  prot
}



# Set prot components ------------------------------------------------

#' Create an empty prot object
#'
#' @noRd
#'

create_enm <- function() {
  prot <- lst(param = NA, nodes = NA, graph = NA, eij = NA, kmat = NA, nma = NA)
  class(prot) <- c("prot", class(prot))
  prot
}

#' Set param of prot object
#'
#' @noRd
#'
set_enm_param <- function(prot, node, model, d_max) {
  prot$param <- lst(node, model, d_max)
  prot
}


#' Set nodes of prot object
#'
#' @noRd
#'
set_enm_nodes <- function(prot, pdb) {
  prot$nodes <- calculate_enm_nodes(pdb, get_enm_node(prot))
  return(prot)
}


#' Set graph of prot object
#'
#' @noRd
#'
set_enm_graph <- function(prot) {
  prot$graph <- calculate_enm_graph(get_xyz(prot), get_pdb_site(prot), get_enm_model(prot), get_d_max(prot))
  prot
}

#' Set eij unit vectors of prot object
#'
#' @noRd
#'
set_enm_eij <- function(prot) {
  prot$eij <- calculate_enm_eij(get_xyz(prot), get_graph(prot)$i, get_graph(prot)$j)
  prot
}

#' Set enm's kmat of prot object
#'
#' @noRd
#'
set_enm_kmat <- function(prot) {
  prot$kmat <- calculate_enm_kmat(get_graph(prot), get_eij(prot), get_nsites(prot))
  prot
}

#' Set normal-mode-analysis component of prot object
#'
#' @noRd
#'
set_enm_nma <- function(prot) {
  prot$nma <- calculate_enm_nma(get_kmat(prot))
  prot
}


# Calculate prot components -----------------------------------------------


#' Calculate nodes of prot object
#'
#' @param pdb pdb object obtained using bio3d::read.pdb()
#' @param node type, either "ca" or "sc"
#'
#' @returns a list of node properties:  \code{lst(nsites, site, pdb_site, bfactor, xyz)}
#'
#'@family enm builders
#' @noRd
#'
calculate_enm_nodes <- function(pdb, node) {
  if (node == "calpha" | node == "ca") {
    nodes <- prot_ca(pdb)
    return(nodes)
  }
  if (node == "side_chain" | node == "sc") {
    nodes <- prot_sc(pdb)
    return(nodes)
  }
  if (node == "cb" | node == "beta") {
    nodes <- prot_cb(pdb)
    return(nodes)
  }
  stop("Error: node must be ca, calpha, sc, side_chain, cb, or beta")
}


#' Calculate ENM graph
#'
#' Calculates graph representation of Elastic Network Model (ENM), the typical relaxed case (lij = dij)
#'
#' @param xyz matrix of size \code{c(3,N)} containing each column the \code{x, y, z} coordinates of each of N nodes
#' @param pdb_site integer vector of size N containing the number of each node (pdb residue number)
#' @param model  character variable specifying the ENM model variant, default is \code{"gnm"}, options:
#'     \code{gnm, anm, ming_wall, hnm, hnm0, pfgnm, reach}.
#' @param d_max distance-cutoff to define network contacts
#' @return a tibble that contains the graph representation of the network
#'
#' @examples
#' \dontrun{
#'  calculate_enm_graph(xyz, pdb_site, model, d_max)
#' }
#'
#'@family enm builders
#' @noRd
#'
calculate_enm_graph <- function(xyz, pdb_site, model, d_max, ...) {
    # Calculate (relaxed) enm graph from xyz
    # Returns graph for the relaxed case

    # put xyz in the right format and check size
    xyz <- my_as_xyz(xyz)
    nsites <- length(pdb_site)
    stopifnot(ncol(xyz) == nsites)

    # set function to calculate i-j spring constants
    kij_fun <- match.fun(paste0("kij_", model))
    kij_par <- lst(d_max = d_max)

    site <- seq(nsites)
    # calculate graph
    graph <- as_tibble(expand_grid(i = site, j = site)) %>%
      filter(j > i) %>%
      arrange(i, j) %>%
      mutate(dij = dij_edge(xyz, i, j)) %>%
      mutate(sdij = sdij_edge(pdb_site, i, j)) %>%
      filter(dij <= d_max | sdij == 1) %>%
      mutate(lij = dij)

    graph$kij <- do.call(kij_fun,
                         c(lst(
                           dij = graph$dij, sdij = graph$sdij
                         ), kij_par))

    graph <- graph %>%
      mutate(edge = paste(i, j, sep = "-"),
             lij = dij) %>%
      mutate(v0ij = 0) %>%
      dplyr::select(edge, i, j, v0ij, sdij, lij, kij, dij)

    graph
  }

#' Calculate distance of edges
#'
#' @noRd
#'
dij_edge <- function(xyz, i, j) {
  stopifnot(length(i) == length(j))
  xyz <- my_as_xyz(xyz)
  dij <- rep(NA, length(i))
  for (k in seq(length(i)))  {
    ik <- i[k]
    jk <- j[k]
    rij <- xyz[, jk] - xyz[, ik]
    dij[k] <- sqrt(sum(rij^2))
  }
  dij
}

#' Calculate edge sequence distance
#'
#' @noRd
#'
sdij_edge <- function(pdb_site, i, j) {
  # sequence distance
  stopifnot(length(i) == length(j))
  sdij <- abs(pdb_site[j] - pdb_site[i])
  sdij
}


#' Calculate unit vectors of edges
#'
#' @param i,j integer vectors of nodes connected in each edge
#' @param xyz vector of xyz coordinates
#' @return matrix with n_edge rows and 3 columns (x, y, z)
#'
#' @family enm builders
#' @noRd
#'
calculate_enm_eij <- function(xyz, i, j) {

  stopifnot(length(i) == length(j))
  n_edges <- length(i)
  xyz <- my_as_xyz(xyz)
  # eij <- tibble(eij_x = rep(NA,n_edges),
  #               eij_y = rep(NA,n_edges),
  #               eij_z = rep(NA,n_edges))

  eij <- matrix(NA, n_edges, 3)
  for (k in seq(n_edges)) {
    ik <- i[k]
    jk <- j[k]
    rij <- xyz[, jk] - xyz[, ik]
    eij[k, ] <- rij / sqrt(sum(rij^2))
  }
  eij
}


#' Calculate kmat given the ENM graph
#'
#' @param graph A tibble representing the ENM graph (with edge information, especially \code{kij}
#' @param eij A matrix of size \code{n_edges x 3} of \code{eij} versors directed along ENM contacts
#' @param nsites The number of nodes of the ENM network
#'
#' @return The \code{3 nsites x 3 nsites} stiffness matrix of the ENM
#'
#' @examples
#' \dontrun{
#' pdb <- read.pdb("2acy")
#' nodes <- calculate_enm_nodes <- function(pdb, node = "ca")
#' graph <- calculate_enm_graph(nodes$xyz, nodes$pdb_site, model = "anm", d_max = 10.5)
#' eij <- calculate_enm_eij(nodes$xyz, graph$i, graph$j)
#' kmat <- calculate_enm_kmat(graph, eij, nsites)
#' }
#'
#' @family enm builders
#' @noRd
#'
calculate_enm_kmat <- function(graph, eij, nsites) {
  stopifnot(max(graph$i, graph$j) <= nsites,
            nrow(graph) == nrow(eij))
  kmat <- array(0, dim = c(3, nsites, 3, nsites))
  for (edge in seq(nrow(graph))) {
    i <- graph$i[[edge]]
    j <- graph$j[[edge]]
    kij <- graph$kij[[edge]]
    eij_v <- eij[edge, ]
    eij_mat <- tcrossprod(eij_v, eij_v)
    kij_mat <- -kij * eij_mat
    kmat[, j, , i] <- kmat[, i, , j] <- kij_mat
  }
  for (i in seq(nsites)) {
    kmat[, i, , i] <- -apply(kmat[, i, , -i], c(1, 2), sum)
  }

  dim(kmat) <- c(3 * nsites, 3 * nsites)
  kmat
}


#' Perform Normal Mode Analysis
#'
#' Given an enm `kmat`, perform NMA
#'
#' @param kmat The K matrix to diagonalize
#' @param too_small=1.e-5 A small value, eigenvectors with eigenvalues larger than `too_small` are discarded
#'
#' @return A list with elements \code{lst(mode,evalue,cmat,umat)}
#'
#' @examples
#' \dontrun{
#' calculate_enm_anm(kmat, too_small = 1.e-10)
#' }
#'
#'@family enm builders
#' @noRd
#'
calculate_enm_nma <- function(kmat, too_small = 1.e-5) {
  eig <- eigen(kmat, symmetric = TRUE)
  evalue <- eig$values
  umat <- eig$vectors
  modes <- evalue > too_small
  evalue <- evalue[modes]
  umat  <- umat[, modes]

  nmodes <- sum(modes)
  mode <- order(seq(nmodes), decreasing = T)
  evalue <- evalue[mode]
  umat <- umat[, mode]
  mode <- mode[mode]


  cmat <-  umat %*% ((1 / evalue) * t(umat))

  nma <- list(
    mode = mode,
    evalue = evalue,
    cmat = cmat,
    umat = umat
  )
  nma
}
