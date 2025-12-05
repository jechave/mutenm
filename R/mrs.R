#' Mutation Response Scanning
#'
#' One-stop function for mutation response analysis. Takes a PDB structure,
#' builds an ENM, scans all sites with mutations, and returns the response
#' matrix along with influence and sensitivity profiles.
#'
#' @details
#' Mutation Response Scanning (MRS) provides a complete workflow:
#' \enumerate{
#'   \item Build an Elastic Network Model from the PDB structure
#'   \item For each site, simulate \code{nmut} random mutations
#'   \item Calculate the squared displacement at every site for each mutation
#'   \item Average over mutations to get the response matrix
#'   \item Compute influence (column means) and sensitivity (row means) profiles
#' }
#'
#' The response matrix \code{dr2ij[i,j]} represents the mean squared displacement
#' at site \code{i} caused by mutations at site \code{j}.
#'
#' \subsection{Influence and Sensitivity}{
#' \itemize{
#'   \item \strong{Influence}: How much does mutating site j affect the whole protein?
#'         (column means of dr2ij)
#'   \item \strong{Sensitivity}: How much does site i respond to mutations anywhere?
#'         (row means of dr2ij)
#' }
#' }
#'
#' @param pdb A pdb object from \code{bio3d::read.pdb()}
#' @param node Node representation: "ca" (C-alpha), "cb" (C-beta), or "sc" (side-chain centroid)
#' @param model ENM model: "anm", "ming_wall", "pfanm", "reach", "hnm", or "hnm0"
#' @param d_max Distance cutoff for contacts (Angstroms)
#' @param mut_model Mutation model: "lfenm" (fast) or "sclfenm" (recalculates ENM)
#' @param nmut Number of mutations to simulate per site
#' @param mut_dl_sigma Standard deviation for bond length perturbations
#' @param mut_sd_min Minimum sequence distance of contacts to perturb
#' @param seed Random seed for reproducibility
#'
#' @return A list with components:
#'   \describe{
#'     \item{dr2ij}{Response matrix (nsites x nsites). Element [i,j] is the mean
#'           squared displacement at site i due to mutation at site j.}
#'     \item{influence}{Vector of length nsites. Column means of dr2ij —
#'           the average effect of mutating each site.}
#'     \item{sensitivity}{Vector of length nsites. Row means of dr2ij —
#'           the average response of each site to mutations.}
#'     \item{params}{List of parameters used (for reproducibility and plotting).}
#'   }
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' result <- mrs(pdb)
#' # View response matrix
#' image(log10(result$dr2ij))
#' # Plot profiles
#' plot(result$influence, type = "l")
#' }
#'
#' @family core functions
#' @seealso \code{\link{enm}}, \code{\link{mutenm}}, \code{\link{plot_mrs}}
#'
#' @export
mrs <- function(pdb,
                node = "ca",
                model = "anm",
                d_max = 10.5,
                mut_model = "lfenm",
                nmut = 1,
                mut_dl_sigma = 0.3,
                mut_sd_min = 2,
                seed = NULL) {

  # Build ENM

  wt <- enm(pdb, node = node, model = model, d_max = d_max)

  # Setup
  if (!is.null(seed)) set.seed(seed)
  nsites <- get_nsites(wt)
  sites <- get_site(wt)


  # Initialize response matrix
  dr2ij <- matrix(0, nsites, nsites)

  # Main loop: for each mutation site j

  for (jj in seq_along(sites)) {
    j <- sites[jj]

    # Accumulator for mutations at site j
    dr2i_accum <- matrix(0, nsites, nmut)

    for (m in seq_len(nmut)) {
      mut <- mutenm(wt, j, m, mut_model, mut_dl_sigma, mut_sd_min, seed)

      # Calculate squared displacement at each site
      dxyz <- mut$nodes$xyz - wt$nodes$xyz
      dr2i_accum[, m] <- colSums(matrix(dxyz^2, nrow = 3))

      rm(mut)
    }

    # Average over mutations for site j
    dr2ij[, jj] <- rowMeans(dr2i_accum)

    # Periodic garbage collection
    if (jj %% 10 == 0) gc()
  }

  # Calculate profiles (means, not sums)
  influence <- colMeans(dr2ij)
  sensitivity <- rowMeans(dr2ij)

  # Return results
  list(
    dr2ij = dr2ij,
    influence = influence,
    sensitivity = sensitivity,
    params = list(
      node = node,
      model = model,
      d_max = d_max,
      mut_model = mut_model,
      nmut = nmut,
      mut_dl_sigma = mut_dl_sigma,
      mut_sd_min = mut_sd_min,
      seed = seed,
      nsites = nsites
    )
  )
}
