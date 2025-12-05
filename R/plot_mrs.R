# Suppress R CMD check notes for ggplot2 non-standard evaluation
utils::globalVariables(c("i", "j", "value", "site"))

#' Plot Mutation Response Scanning Results
#'
#' Creates a composite plot showing the response matrix as a heatmap with
#' influence and sensitivity profiles in separate panels.
#'
#' @details
#' The plot consists of three panels arranged as: heatmap | (influence / sensitivity)
#' \itemize{
#'   \item \strong{Left}: Heatmap of dr2ij (log10 color scale) showing the response matrix
#'   \item \strong{Top right}: Influence profile (log10 scale) — column means (effect of mutating each site)
#'   \item \strong{Bottom right}: Sensitivity profile (log10 scale) — row means (response of each site)
#' }
#'
#' Log10 scales are used for visualization because the response values
#' span several orders of magnitude.
#'
#' @param mrs_result Output from \code{\link{mrs}}
#'
#' @return A ggplot object (assembled with patchwork)
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' result <- mrs(pdb, pdb_id = "2ACY")
#' plot_mrs(result)
#'
#' # Save as PDF (recommended for publications)
#' p <- plot_mrs(result)
#' ggplot2::ggsave("mrs_figure.pdf", p, width = 8, height = 5)
#' }
#'
#' @family core functions
#' @seealso \code{\link{mrs}}
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_line coord_fixed
#' @importFrom ggplot2 scale_fill_viridis_c scale_x_continuous scale_y_continuous scale_y_reverse scale_y_log10
#' @importFrom ggplot2 theme_minimal theme_bw theme element_blank element_text element_line margin unit
#' @importFrom ggplot2 labs guide_colorbar
#'
#' @export
plot_mrs <- function(mrs_result) {

  # Check for patchwork
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for plot_mrs(). Install it with: install.packages('patchwork')")
  }
  # Check for cowplot
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for plot_mrs(). Install it with: install.packages('cowplot')")
  }

  # Extract data
  dr2ij <- mrs_result$dr2ij
  influence <- mrs_result$influence
  sensitivity <- mrs_result$sensitivity
  nsites <- mrs_result$params$nsites

  # Normalize by mean for better scale
 mean_dr2 <- mean(dr2ij)
  dr2ij_rel <- dr2ij / mean_dr2
  influence_rel <- influence / mean_dr2
  sensitivity_rel <- sensitivity / mean_dr2

  # Create data frames for ggplot
  sites <- seq_len(nsites)

  # Heatmap data (relative values)
  heatmap_df <- expand.grid(i = sites, j = sites)
  heatmap_df$value <- as.vector(dr2ij_rel)

  # Profile data (relative values, will use log scale)
  influence_df <- data.frame(site = sites, value = influence_rel)
  sensitivity_df <- data.frame(site = sites, value = sensitivity_rel)

  # --- Heatmap (left) ---
  p_heatmap <- ggplot(heatmap_df, aes(x = j, y = i, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(
      name = expression(Delta*r^2 / mean),
      trans = "log10",
      guide = guide_colorbar(
        title.position = "top",
        barheight = unit(4, "cm")
      )
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    labs(x = "mutation site (j)", y = "response site (i)", title = "Response matrix") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 11, face = "bold")
    )

  # Profile theme: clean with major gridlines
  profile_theme <- theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      plot.title = element_text(size = 11, face = "bold")
    )

  # --- Influence profile (top right) ---
  p_influence <- ggplot(influence_df, aes(x = site, y = value)) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_log10() +
    labs(x = "site", y = "relative influence", title = "Influence") +
    profile_theme

  # --- Sensitivity profile (bottom right) ---
  p_sensitivity <- ggplot(sensitivity_df, aes(x = site, y = value)) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_log10() +
    labs(x = "site", y = "relative sensitivity", title = "Sensitivity") +
    profile_theme


  # --- Combine with patchwork: heatmap | (influence / sensitivity) ---
  p_combined <- p_heatmap | (p_influence / p_sensitivity)

  # Add title (top, centered) and caption with key model parameters (bottom)
  params <- mrs_result$params
  plot_title <- paste0("Mutation response analysis of ", params$pdb_id)

  # Key model parameters only
  plot_caption <- paste0(
    "ENM: ", params$model, " | nodes: ", params$node,
    " | mutation model: ", params$mut_model
  )

  p_combined <- p_combined +
    patchwork::plot_annotation(
      title = plot_title,
      caption = plot_caption,
      theme = ggplot2::theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.caption = element_text(size = 9, color = "grey30", hjust = 0.5)
      )
    )

  p_combined
}
