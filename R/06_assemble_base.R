# EXTENDED DATA FIG. 7 - the title, subtitle and legend layer.
#
# The finished figure is this text layer with the spatial-map panel merged onto it by python/compose.py.
# The panel is the `pmap` figure written by R/04_abundance_proximity.R as `04_spatial_maps.pdf`.
#
# Usage: Rscript R/06_assemble_base.R
suppressPackageStartupMessages({library(grid)})
INK <- "#1A1A1A"; MUTED <- "#5A5A5A"
OUT <- Sys.getenv("ED7_OUT", file.path(normalizePath(file.path(
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), mustWork = FALSE), "figures"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
pt  <- function(x) unit(x, "bigpts")
txt <- function(s, x, y, size = 7, face = "plain", col = MUTED, hjust = "left", lineheight = 1.25)
  grid.text(s, x = pt(x), y = pt(y), just = c(hjust, "bottom"),
            gp = gpar(fontfamily = "Helvetica", fontsize = size, fontface = face, col = col, lineheight = lineheight))

f <- file.path(OUT, "ED7 base.pdf")
cairo_pdf(f, width = 612 / 72, height = 430 / 72, onefile = FALSE); grid.newpage()
txt("Extended Data Fig. 7 | PI16+ stroma at the tumour-epithelial interface", 30, 406, 9.5, "bold", INK)
txt("Spatial support for the interface section of the Results.", 30, 395, 6.6)
txt(paste(
  "a, Spatial distribution of the adipose-progenitor-like stromal state and of the tumour populations across all available slides, three adenocarcinoma and two neuroendocrine. Each panel is",
  "translated to its own origin and all panels share one millimetre scale, so area and spacing are comparable between slides; absolute positions are not. Colour is the per-cell population call:",
  "ASPC, the adipose-progenitor-like state; EPI_CA, adenocarcinoma epithelium; NE_CA, neuroendocrine epithelium; Multi, cells assigned to more than one population; Other, the remainder.", sep = "\n"),
  30, 12, 5.8)
dev.off()
cat("written", f, "\n")
