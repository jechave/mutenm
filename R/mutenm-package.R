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
#'   \item Analyze the effects using property functions
#' }
#' }
#'
#' \subsection{Core functions}{
#' \itemize{
#'   \item \code{\link{enm}}: Build an elastic network model from a PDB structure
#'   \item \code{\link{mutenm}}: Simulate a mutation at a specific site
#'   \item \code{\link{mrs}}: Mutation Response Scanning - scan all sites
#' }
#' }
#'
#' \subsection{Single-protein properties}{
#' Functions that calculate properties of one protein:
#' \itemize{
#'   \item Structure: \code{\link{cn}}, \code{\link{wcn}}, \code{\link{dactive}}
#'   \item Dynamics: \code{\link{msfi}}, \code{\link{msfn}}
#'   \item Energy: \code{\link{v_min}}, \code{\link{g_ent}}, \code{\link{dv_act}},
#'         \code{\link{dg_ent_act}}
#' }
#' }
#'
#' \subsection{Mutation effects (pair comparison)}{
#' Functions that compare wild-type and mutant. The \code{D} prefix stands for
#' "delta" (change):
#' \itemize{
#'   \item Structure: \code{\link{Dr2i}}, \code{\link{Dr2n}}
#'   \item Dynamics: \code{\link{Dmsfi}}, \code{\link{Dmsfn}}
#'   \item Energy: \code{\link{Dv_min}}, \code{\link{Dg_ent}}, \code{\link{Ddv_act}},
#'         \code{\link{Ddg_ent_act}}
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
#' @param beta Inverse temperature (1/kT). Default: beta_boltzmann()
#' @param pdb_site_active Vector of PDB residue numbers defining the active site
#' @param ideal Reference protein for ideal active site conformation
#'
#' @name mutenm-params
#' @keywords internal
NULL

