# Generate figure for JOSS paper
# Focus on dynamics changes - what PRS cannot do

library(mutenm)
library(ggplot2)
library(patchwork)

# Build ENM from example protein
pdb <- bio3d::read.pdb("inst/extdata/2acy.pdb")
wt <- enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5)
sites <- wt$nodes$pdb_site
nsites <- length(sites)

# -----------------------------------------------------------------------------
# Panel A: Wild-type fluctuation profile (log scale)
# -----------------------------------------------------------------------------
msfi_wt <- msfi(wt)
df_a <- data.frame(site = sites, msfi = msfi_wt)

panel_a <- ggplot(df_a, aes(x = site, y = log10(msfi))) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  labs(x = "Residue", y = expression(log[10](MSF[i])),
       title = "A) Wild-type fluctuations") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10))

# -----------------------------------------------------------------------------
# Panel B: Single mutation effect - show Dmsfi (difference)
# Pick a site with intermediate CN for visible effect
# -----------------------------------------------------------------------------
cn_vals <- cn(wt)
# Find site with CN closest to median
median_cn <- median(cn_vals)
site_mut <- which.min(abs(cn_vals - median_cn))
cat("Mutating site:", sites[site_mut], "(intermediate CN =", cn_vals[site_mut], ")\n")

mut <- mutenm(wt, site_mut = site_mut, mutation = 1,
              mut_model = "sclfenm", seed = 42)

dmsfi_single <- Dmsfi(wt, mut)
# Signed log transform
dmsfi_signed_log <- sign(dmsfi_single) * log10(abs(dmsfi_single) + 1e-10)

df_b <- data.frame(site = sites, value = dmsfi_signed_log)

panel_b <- ggplot(df_b, aes(x = site, y = value)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_vline(xintercept = sites[site_mut], linetype = "dashed", color = "darkred") +
  labs(x = "Residue", y = expression(sign %.% log[10](abs(Delta*MSF[i]))),
       title = paste0("B) Mutation at site ", sites[site_mut])) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10))

# -----------------------------------------------------------------------------
# Panel C: MRS dynamics response (Dmsfi matrix + profiles)
# This is slow - run once and save
# -----------------------------------------------------------------------------
mrs_file <- "paper_mrs_results.rds"
if (file.exists(mrs_file)) {
  cat("Loading cached MRS results...\n")
  responses <- readRDS(mrs_file)
} else {
  cat("Running MRS (nmut=10, sclfenm) - this will take a while...\n")
  responses <- mrs(wt, nmut = 10, mut_model = "sclfenm", seed = 42)
  saveRDS(responses, mrs_file)
  cat("MRS results saved to", mrs_file, "\n")
}

# Dmsfi matrix
dmsfi_df <- responses$Dmsfi

# Heatmap (log scale, preserve sign)
# Use sign(Dmsfi) * log10(abs(Dmsfi))
dmsfi_df$signed_log <- sign(dmsfi_df$Dmsfi) * log10(abs(dmsfi_df$Dmsfi) + 1e-10)

panel_c1 <- ggplot(dmsfi_df, aes(x = j, y = i, fill = signed_log)) +
  geom_tile() +
  scale_fill_gradient2(name = expression(sign %.% log[10](abs(Delta*MSF[i]))),
                       low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_continuous(breaks = seq(1, nsites, by = 20),
                     labels = sites[seq(1, nsites, by = 20)]) +
  scale_y_continuous(breaks = seq(1, nsites, by = 20),
                     labels = sites[seq(1, nsites, by = 20)]) +
  coord_fixed() +
  labs(x = "Mutation site", y = "Response site",
       title = "C) Dynamics response matrix") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 6),
        legend.key.height = unit(0.4, "cm"))

# Sensitivity (column mean) and Influence (row mean) profiles
# First average Dmsfi, then apply sign * log10(abs())
sensitivity <- aggregate(Dmsfi ~ i, data = dmsfi_df, FUN = mean)
influence <- aggregate(Dmsfi ~ j, data = dmsfi_df, FUN = mean)

# Apply signed log transformation
sens_signed_log <- sign(sensitivity$Dmsfi) * log10(abs(sensitivity$Dmsfi) + 1e-10)
infl_signed_log <- sign(influence$Dmsfi) * log10(abs(influence$Dmsfi) + 1e-10)

df_profiles <- data.frame(
  site = c(sites, sites),
  value = c(sens_signed_log, infl_signed_log),
  type = rep(c("Sensitivity", "Influence"), each = nsites)
)

panel_c2 <- ggplot(df_profiles, aes(x = site, y = value)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  labs(x = "Residue", y = expression(sign %.% log[10](abs(Mean~Delta*MSF))),
       title = "D) Sensitivity & Influence") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10),
        strip.text = element_text(face = "bold"))

# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------
fig <- (panel_a | panel_b) / (panel_c1 | panel_c2) +
  plot_annotation(
    caption = "Acylphosphatase (PDB: 2ACY)"
  )

ggsave("paper_figure.png", fig, width = 8, height = 7, dpi = 300)
cat("Figure saved to paper_figure.png\n")
