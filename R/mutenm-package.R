#' @description
#' Build elastic network models (ENMs) of proteins and simulate mutation effects
#' using linear response theory.
#'
#' @details
#' \subsection{Workflow}{
#' A typical analysis has three steps:
#' \enumerate{
#'   \item Build the ENM from a PDB structure using \code{\link{enm}}
#'   \item Simulate mutations using \code{\link{mutenm}} for a single mutation,
#'         or \code{\link{mrs}} for scanning all sites
#'   \item Visualize results using \code{\link{plot_mrs}}
#' }
#' }
#'
#' \subsection{Functions}{
#' \itemize{
#'   \item \code{\link{enm}}: Build an elastic network model from a PDB structure
#'   \item \code{\link{mutenm}}: Simulate a mutation at a specific site
#'   \item \code{\link{mrs}}: Mutation Response Scanning - scan all sites
#'   \item \code{\link{plot_mrs}}: Visualize MRS results
#' }
#' }
#'
#' @aliases NULL
"_PACKAGE"

#' Common parameters for mutenm functions
#'
#' @param prot A protein object created by [enm()]
#' @param wt Wild-type protein object created by [enm()]
#' @param mut Mutant protein object created by [mutenm()]
#'
#' @name mutenm-params
#' @keywords internal
NULL
