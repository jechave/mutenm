#' Mutation Response Scanning
#'
#' Calculate mutation response matrices using process-and-discard approach.
#' Memory efficient: generates one mutant at a time, calculates responses, discards.
#'
#' @param wt Wild-type protein object from enm()
#' @param nmut Number of mutations per site
#' @param mut_model Mutation model: "lfenm" (fast) or "sclfenm" (recalculates ENM)
#' @param mut_dl_sigma Sigma of normal distribution for bond length perturbations
#' @param mut_sd_min Minimum sequence distance of contacts to perturb
#' @param seed Random seed for reproducibility
#'
#' @return List with response tibbles and params:
#'   \itemize{
#'     \item Dr2i: tibble (i, j, Dr2i) - structural change at site i due to mutation at site j
#'     \item Dmsfi: tibble (i, j, Dmsfi) - dynamics change at site i due to mutation at site j
#'     \item Dr2n: tibble (n, j, Dr2n) - structural change in mode n due to mutation at site j
#'     \item Dmsfn: tibble (n, j, Dmsfn) - dynamics change in mode n due to mutation at site j
#'     \item Dv_min: tibble (j, Dv_min) - minimum energy change due to mutation at site j
#'     \item Dg_ent: tibble (j, Dg_ent) - entropic free energy change due to mutation at site j
#'     \item Ddv_act: tibble (j, Ddv_act) - activation energy change due to mutation at site j
#'     \item Ddg_ent_act: tibble (j, Ddg_ent_act) - activation entropy change due to mutation at site j
#'     \item params: list with enm_param and mut_param
#'   }
#'
#' @export
mrs <- function(wt, nmut, mut_model = "lfenm", mut_dl_sigma = 0.3,
                mut_sd_min = 2, seed = NULL) {

  # Warning for lfenm model
  if (mut_model == "lfenm") {
    warning("lfenm: dynamics unchanged, Dmsfi/Dmsfn/Dg_ent/Ddg_ent_act will be zero")
  }

  # Setup
  if (!is.null(seed)) set.seed(seed)
  nsites <- get_nsites(wt)
  sites <- get_site(wt)
  nmodes <- get_nmodes(wt)

  # Initialize accumulators (matrices for site/mode responses, vectors for scalars)
  Dr2i_mat <- matrix(0, nsites, nsites)
  Dmsfi_mat <- matrix(0, nsites, nsites)
  Dr2n_mat <- matrix(0, nmodes, nsites)
  Dmsfn_mat <- matrix(0, nmodes, nsites)
  Dv_min_vec <- numeric(nsites)
  Dg_ent_vec <- numeric(nsites)
  Ddv_act_vec <- numeric(nsites)
  Ddg_ent_act_vec <- numeric(nsites)

  # Main loop: process and discard
  for (jj in seq_along(sites)) {
    j <- sites[jj]

    # Accumulators for mutations at site j
    Dr2i_accum <- matrix(0, nsites, nmut)
    Dmsfi_accum <- matrix(0, nsites, nmut)
    Dr2n_accum <- matrix(0, nmodes, nmut)
    Dmsfn_accum <- matrix(0, nmodes, nmut)
    Dv_min_accum <- numeric(nmut)
    Dg_ent_accum <- numeric(nmut)
    Ddv_act_accum <- numeric(nmut)
    Ddg_ent_act_accum <- numeric(nmut)

    for (m in seq_len(nmut)) {
      mut <- mutenm(wt, j, m, mut_model, mut_dl_sigma, mut_sd_min, seed)

      # Site responses
      Dr2i_accum[, m] <- Dr2i(wt, mut)
      Dmsfi_accum[, m] <- Dmsfi(wt, mut)

      # Mode responses
      Dr2n_accum[, m] <- Dr2n(wt, mut)
      Dmsfn_accum[, m] <- Dmsfn(wt, mut)

      # Scalar responses
      Dv_min_accum[m] <- Dv_min(wt, mut)
      Dg_ent_accum[m] <- Dg_ent(wt, mut)
      Ddv_act_accum[m] <- Ddv_act(wt, mut)
      Ddg_ent_act_accum[m] <- Ddg_ent_act(wt, mut)

      rm(mut)
    }

    # Average over mutations for site j
    Dr2i_mat[, jj] <- rowMeans(Dr2i_accum)
    Dmsfi_mat[, jj] <- rowMeans(Dmsfi_accum)
    Dr2n_mat[, jj] <- rowMeans(Dr2n_accum)
    Dmsfn_mat[, jj] <- rowMeans(Dmsfn_accum)
    Dv_min_vec[jj] <- mean(Dv_min_accum)
    Dg_ent_vec[jj] <- mean(Dg_ent_accum)
    Ddv_act_vec[jj] <- mean(Ddv_act_accum)
    Ddg_ent_act_vec[jj] <- mean(Ddg_ent_act_accum)

    if (jj %% 10 == 0) gc()
  }

  # Convert to tibbles and return
  list(
    Dr2i = matrix_to_tibble(Dr2i_mat, "i", "j", "Dr2i"),
    Dmsfi = matrix_to_tibble(Dmsfi_mat, "i", "j", "Dmsfi"),
    Dr2n = matrix_to_tibble(Dr2n_mat, "n", "j", "Dr2n"),
    Dmsfn = matrix_to_tibble(Dmsfn_mat, "n", "j", "Dmsfn"),
    Dv_min = tibble(j = sites, Dv_min = Dv_min_vec),
    Dg_ent = tibble(j = sites, Dg_ent = Dg_ent_vec),
    Ddv_act = tibble(j = sites, Ddv_act = Ddv_act_vec),
    Ddg_ent_act = tibble(j = sites, Ddg_ent_act = Ddg_ent_act_vec),
    params = list(
      enm_param = get_enm_param(wt),
      mut_param = list(nmut = nmut, mut_model = mut_model,
                       mut_dl_sigma = mut_dl_sigma, mut_sd_min = mut_sd_min)
    )
  )
}
