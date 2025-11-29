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
#' @param responses Character vector of responses to calculate
#' @param seed Random seed for reproducibility
#'
#' @return List with requested response tibbles and params
#'
#' @importFrom tidyr pivot_longer
#'
#' @export
mrs <- function(wt, nmut, mut_model = "lfenm", mut_dl_sigma = 0.3,
                mut_sd_min = 2, responses = c("dr2ij"), seed = NULL) {

  # Response categories
  site_responses <- c("dr2ij", "dmsfij")
  mode_responses <- c("dr2nj", "dmsfnj")
  scalar_responses <- c("dv", "tds", "dvs", "act_dv", "act_tds")
  motion_responses <- c("dmsfij", "dmsfnj")
  all_valid <- c(site_responses, mode_responses, scalar_responses)

  # Validate responses
  invalid <- setdiff(responses, all_valid)
  if (length(invalid) > 0) {
    stop("Invalid responses: ", paste(invalid, collapse = ", "),
         "\nValid: ", paste(all_valid, collapse = ", "))
  }

  # Warn if motion requested with lfenm
  requested_motion <- intersect(responses, motion_responses)
  if (mut_model == "lfenm" && length(requested_motion) > 0) {
    warning("Motion responses require mut_model='sclfenm'. Skipping: ",
            paste(requested_motion, collapse = ", "))
    responses <- setdiff(responses, motion_responses)
  }

  if (length(responses) == 0) {
    stop("No valid responses to calculate")
  }

  # Setup
  if (!is.null(seed)) set.seed(seed)
  nsites <- get_nsites(wt)
  sites <- get_site(wt)
  nmodes <- get_nmodes(wt)

  # Identify which response types are requested
  req_site <- intersect(responses, site_responses)
  req_mode <- intersect(responses, mode_responses)
  req_scalar <- intersect(responses, scalar_responses)

  # Initialize accumulators
  site_results <- lapply(req_site, function(r) matrix(0, nrow = nsites, ncol = nsites))
  names(site_results) <- req_site

  mode_results <- lapply(req_mode, function(r) matrix(0, nrow = nmodes, ncol = nsites))
  names(mode_results) <- req_mode

  scalar_results <- lapply(req_scalar, function(r) numeric(nsites))
  names(scalar_results) <- req_scalar

  # Main loop: process and discard
  for (jj in seq_along(sites)) {
    j <- sites[jj]

    # Accumulate over mutations at site j
    site_accum <- lapply(req_site, function(r) matrix(0, nrow = nsites, ncol = nmut))
    names(site_accum) <- req_site

    mode_accum <- lapply(req_mode, function(r) matrix(0, nrow = nmodes, ncol = nmut))
    names(mode_accum) <- req_mode

    scalar_accum <- lapply(req_scalar, function(r) numeric(nmut))
    names(scalar_accum) <- req_scalar

    for (m in seq_len(nmut)) {
      mut <- mutenm(wt, j, m, mut_model, mut_dl_sigma, mut_sd_min, seed)

      # Site responses
      if ("dr2ij" %in% req_site) site_accum$dr2ij[, m] <- delta_structure_dr2i(wt, mut)
      if ("dmsfij" %in% req_site) site_accum$dmsfij[, m] <- delta_motion_dmsfi(wt, mut)

      # Mode responses
      if ("dr2nj" %in% req_mode) mode_accum$dr2nj[, m] <- delta_structure_dr2n(wt, mut)
      if ("dmsfnj" %in% req_mode) mode_accum$dmsfnj[, m] <- delta_motion_dmsfn(wt, mut)

      # Scalar responses
      if ("dv" %in% req_scalar) scalar_accum$dv[m] <- Dv_min(wt, mut)
      if ("tds" %in% req_scalar) scalar_accum$tds[m] <- Dg_ent(wt, mut)
      if ("dvs" %in% req_scalar) scalar_accum$dvs[m] <- delta_energy_dvs(wt, mut)
      if ("act_dv" %in% req_scalar) scalar_accum$act_dv[m] <- delta_energy_act_dv(wt, mut)
      if ("act_tds" %in% req_scalar) scalar_accum$act_tds[m] <- delta_energy_act_tds(wt, mut)

      rm(mut)
    }

    # Average over mutations for site j
    for (resp in req_site) site_results[[resp]][, jj] <- rowMeans(site_accum[[resp]])
    for (resp in req_mode) mode_results[[resp]][, jj] <- rowMeans(mode_accum[[resp]])
    for (resp in req_scalar) scalar_results[[resp]][jj] <- mean(scalar_accum[[resp]])

    if (jj %% 10 == 0) gc()
  }

  # Convert to tibbles
  output <- list()

  for (resp in req_site) {
    colnames(site_results[[resp]]) <- sites
    output[[resp]] <- site_results[[resp]] %>%
      as_tibble() %>%
      mutate(i = sites) %>%
      pivot_longer(-i, names_to = "j", values_to = resp) %>%
      mutate(j = as.integer(j))
  }

  for (resp in req_mode) {
    colnames(mode_results[[resp]]) <- sites
    output[[resp]] <- mode_results[[resp]] %>%
      as_tibble() %>%
      mutate(n = seq_len(nmodes)) %>%
      pivot_longer(-n, names_to = "j", values_to = resp) %>%
      mutate(j = as.integer(j))
  }

  for (resp in req_scalar) {
    output[[resp]] <- tibble(j = sites, !!resp := scalar_results[[resp]])
  }

  output$params <- list(
    enm_param = get_enm_param(wt),
    mut_param = list(nmut = nmut, mut_model = mut_model,
                     mut_dl_sigma = mut_dl_sigma, mut_sd_min = mut_sd_min)
  )

  output
}
