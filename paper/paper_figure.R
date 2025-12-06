# Generate figure for JOSS paper
# Uses the new 4-function API: enm(), mutenm(), mrs(), plot_mrs()

library(mutenm)
library(ggplot2)

# Run MRS on acylphosphatase
pdb <- bio3d::read.pdb("inst/extdata/2acy.pdb")
result <- mrs(pdb, pdb_id = "2ACY", nmut = 10, seed = 42)

# Generate figure using plot_mrs()
p <- plot_mrs(result)

# Save as PNG for paper
ggsave("paper_figure.png", p, width = 8, height = 5, dpi = 300)
cat("Figure saved to paper_figure.png\n")

# Also save as PDF
ggsave("paper_figure.pdf", p, width = 8, height = 5)
cat("Figure saved to paper_figure.pdf\n")
