###############################################################################
# 05_visualize.R  v12
#
# Per-biopsy visualisation of ASPC, epithelial cancer (EPI_CA), neuroendocrine
# cancer (NE_CA) and the ASPC neighbourhood, drawn on the actual image.
#
# For each slide one figure (PNG + page in 05_biopsy_atlas.pdf):
#   A  whole-biopsy cell map over the DAPI pyramid
#   B  ASPC neighbourhood map: each ASPC coloured by what dominates its
#      RADIUS_UM neighbourhood (from 04 aspc_neighbourhood_cells.rds)
#   C  full-resolution false-colour inset at the ASPC hotspot
#      (R = PDGFa, G = panCK, B = DAPI, magenta = DLL3) with
#      cell calls overlaid
#   D  neighbourhood composition around ASPC, observed vs expected
# Every panel title carries the 03 marker confidence (marker vs matched blank
# at the same threshold), so low-confidence calls are visible on the figure.
#
# INPUT  results_ASPC_NEPC/cells_classified.rds, channels.rds,
#        aspc_neighbourhood_cells.rds, aspc_neighbourhood.csv,
#        qc_marker_usability_primary_6um.csv
# OUTPUT results_ASPC_NEPC/plots/05_biopsy_<slide>.png, 05_biopsy_atlas.pdf
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(reticulate)
})

if (!exists("OUT")) OUT <- "results_ASPC_NEPC"
PLOTDIR <- file.path(OUT, "plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

if (!exists("RADIUS_UM")) RADIUS_UM <- 50
if (!exists("INSET_UM")) INSET_UM <- 500
if (!exists("BG_LEVEL")) BG_LEVEL <- 3L          # pyramid level for map backgrounds

use_virtualenv("r-reticulate", required = TRUE)
tf <- import("tifffile", convert = FALSE)
np <- import("numpy", convert = FALSE)
zr <- import("zarr", convert = FALSE)
py <- import_builtins()

GROUP_OF <- c(PM0T32 = "NEPC", PM677_A8 = "NEPC", PM12_B3 = "NEPC",
              PM1024_B = "ADENO", PM423 = "ADENO", PM159T5 = "ADENO")

POP_COLS <- c(
  ASPC   = "#FF3B5C",
  EPI_CA = "#2FB36B",
  NE_CA  = "#2E86DE",
  Multi  = "#A55EEA",
  Other  = "grey70"
)

NBR_COLS <- c(
  "EPI_CA-adjacent" = "#2FB36B",
  "NE_CA-adjacent"  = "#2E86DE",
  "ASPC cluster"    = "#E8384F",
  "mixed tumour"    = "#A55EEA",
  "stroma / other"  = "#F5B041"
)

###############################################################################
## INPUTS
###############################################################################

cells <- readRDS(file.path(OUT, "cells_classified.rds"))
z <- readRDS(file.path(OUT, "channels.rds"))
PX_UM <- if (!is.null(z$px_um)) z$px_um else 0.28
slides <- z$slides
info <- z$info

nbr_cells_file <- file.path(OUT, "aspc_neighbourhood_cells.rds")
nbr_cells <- if (file.exists(nbr_cells_file)) readRDS(nbr_cells_file) else NULL

nbr_file <- file.path(OUT, "aspc_neighbourhood.csv")
nbr <- if (file.exists(nbr_file)) read.csv(nbr_file) else NULL

usab_file <- file.path(OUT, "qc_marker_usability_primary_6um.csv")
usab <- if (file.exists(usab_file)) read.csv(usab_file) else NULL

cells <- cells |>
  mutate(
    pop5 = case_when(
      multi_pop %in% TRUE ~ "Multi",
      ASPC %in% TRUE ~ "ASPC",
      NE_CA %in% TRUE ~ "NE_CA",
      EPI_CA %in% TRUE ~ "EPI_CA",
      TRUE ~ "Other"
    ),
    pop5 = factor(pop5, levels = names(POP_COLS))
  )

###############################################################################
## IMAGE HELPERS
###############################################################################

.open_level <- function(tiff, level) {
  st <- tf$TiffFile(tiff)$series[0]$aszarr(level = as.integer(level))
  zr$open(store = st, mode = "r")
}

.read_plane <- function(arr, c_index, x0, y0, w, h) {
  idx <- tuple(as.integer(c_index - 1L),
               py$slice(as.integer(y0), as.integer(y0 + h)),
               py$slice(as.integer(x0), as.integer(x0 + w)))
  m <- py_to_r(np$asarray(arr$`__getitem__`(idx)))
  storage.mode(m) <- "double"
  m                                     # rows = y, cols = x
}

.stretch <- function(m, lo_q = 0.01, hi_q = 0.998, gamma = 1) {
  v <- m[m > 0]
  if (!length(v)) return(m * 0)
  lo <- stats::quantile(v, lo_q, names = FALSE)
  hi <- stats::quantile(v, hi_q, names = FALSE)
  pmin(pmax((m - lo) / max(hi - lo, 1e-9), 0), 1)^gamma
}

.conf_of <- function(sid, m) {
  if (is.null(usab)) return("")
  r <- usab[usab$slide_id == sid & usab$marker == m, , drop = FALSE]
  if (!nrow(r)) return("")
  conf <- if ("confidence" %in% names(r)) r$confidence[1] else ""
  sprintf("%s %s (%.1fx)", m, conf, r$enrichment_over_blank[1])
}

###############################################################################
## PER-SLIDE FIGURE
###############################################################################

make_slide_figure <- function(sid) {
  k <- which(slides$slide_id == sid)
  tiff <- slides$tiff[k]
  cm <- info[[sid]]$chmap
  d <- cells |> filter(slide_id == sid)
  if (!nrow(d)) return(NULL)

  grp <- GROUP_OF[[sid]]

  ## ---------------- A: whole-biopsy map over DAPI pyramid ------------------
  lv <- .open_level(tiff, BG_LEVEL)
  shp <- as.integer(unlist(py_to_r(lv$shape)))
  f <- slides$sizeY[k] / shp[2]               # full-res px per level px
  bg <- .read_plane(lv, cm$DAPI_AR, 0, 0, shp[3], shp[2])
  bg <- .stretch(bg, 0.02, 0.995, gamma = 0.7)
  bg_rgb <- grDevices::as.raster(array(c(bg * 0.32, bg * 0.32, bg * 0.38), dim = c(dim(bg), 3)))

  W_mm <- slides$sizeX[k] * PX_UM / 1000
  H_mm <- slides$sizeY[k] * PX_UM / 1000

  ## crop the map to the tissue bounding box (+0.3 mm)
  pad <- 0.3
  xr <- range(d$x_px) * PX_UM / 1000 + c(-pad, pad)
  yr <- range(d$y_px) * PX_UM / 1000 + c(-pad, pad)
  xr <- c(max(0, xr[1]), min(W_mm, xr[2]))
  yr <- c(max(0, yr[1]), min(H_mm, yr[2]))

  d <- d |> mutate(x_mm = x_px * PX_UM / 1000, y_mm = y_px * PX_UM / 1000)

  fg <- d |> filter(pop5 != "Other") |> arrange(match(pop5, c("EPI_CA", "NE_CA", "Multi", "ASPC")))

  n_tab <- table(d$pop5)
  ## dot size scales with the mapped extent so calls stay visible on both a
  ## 1 mm core and a 12 mm section
  ext_mm <- max(diff(xr), diff(yr))
  pts <- max(0.12, min(0.9, 2.2 / ext_mm))

  base_map <- function() {
    ggplot() +
      annotation_raster(bg_rgb, xmin = 0, xmax = W_mm, ymin = -H_mm, ymax = 0) +
      coord_equal(xlim = xr, ylim = -rev(yr), expand = FALSE) +
      scale_x_continuous("mm") +
      scale_y_continuous("mm", labels = function(v) -v) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(), legend.position = "bottom",
            plot.title = element_text(face = "bold"))
  }

  pA <- base_map() +
    geom_point(data = fg, aes(x_mm, -y_mm, colour = pop5), size = pts, alpha = 0.9, stroke = 0) +
    scale_colour_manual(NULL, values = POP_COLS, drop = FALSE,
                        labels = function(l) paste0(l, " (", n_tab[l], ")")) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 1)) +
    labs(title = "A  Cell calls (unlabelled cells = DAPI background)",
         subtitle = paste(.conf_of(sid, "PDGF"), "|", .conf_of(sid, "panCK"), "|",
                          .conf_of(sid, "DLL3")))

  ## ---------------- B: ASPC neighbourhood map ------------------------------
  nb <- if (!is.null(nbr_cells)) nbr_cells |> filter(slide_id == sid) else NULL

  if (!is.null(nb) && nrow(nb)) {
    nb <- nb |>
      mutate(
        x_mm = x_px * PX_UM / 1000, y_mm = y_px * PX_UM / 1000,
        tum = frac_EPI_CA + frac_NE_CA + frac_Multi,
        nclass = case_when(
          frac_EPI_CA >= 0.2 & frac_NE_CA >= 0.2 ~ "mixed tumour",
          frac_EPI_CA >= 0.2 ~ "EPI_CA-adjacent",
          frac_NE_CA  >= 0.2 ~ "NE_CA-adjacent",
          frac_ASPC   >= 0.3 ~ "ASPC cluster",
          TRUE ~ "stroma / other"
        ),
        nclass = factor(nclass, levels = names(NBR_COLS))
      )
    nb_tab <- table(nb$nclass)

    tumour <- d |> filter(pop5 %in% c("EPI_CA", "NE_CA", "Multi"))

    pB <- base_map() +
      geom_point(data = tumour, aes(x_mm, -y_mm), colour = "white", size = pts * 0.6, alpha = 0.35, stroke = 0) +
      geom_point(data = nb |> arrange(nclass != "stroma / other"),
                 aes(x_mm, -y_mm, colour = nclass), size = pts * 1.3, alpha = 0.95, stroke = 0) +
      scale_colour_manual(NULL, values = NBR_COLS, drop = FALSE,
                          labels = function(l) paste0(l, " (", nb_tab[l], ")")) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 2)) +
      labs(title = paste0("B  ASPC neighbourhood (", RADIUS_UM, " um)"),
           subtitle = "each ASPC coloured by its neighbours: >=20% EPI_CA or NE_CA, >=30% ASPC\nwhite = tumour cells")
  } else {
    pB <- ggplot() + theme_void() +
      labs(title = "B  ASPC neighbourhood",
           subtitle = "fewer than MIN_N ASPC -- neighbourhood not computed")
    nb <- NULL
  }

  ## ---------------- C: full-resolution inset at the ASPC hotspot -----------
  win_px <- round(INSET_UM / PX_UM)
  asp <- d |> filter(pop5 == "ASPC")
  if (nrow(asp) >= 5) {
    ## hotspot = 250 um bin with most ASPC that also contains tumour cells
    b <- round(250 / PX_UM)
    hs <- d |>
      mutate(bx = floor(x_px / b), by = floor(y_px / b)) |>
      group_by(bx, by) |>
      summarise(n_aspc = sum(pop5 == "ASPC"),
                n_tum = sum(pop5 %in% c("EPI_CA", "NE_CA", "Multi")),
                .groups = "drop") |>
      mutate(score = n_aspc * (1 + pmin(n_tum, n_aspc))) |>
      arrange(desc(score)) |> slice(1)
    cx <- (hs$bx + 0.5) * b; cy <- (hs$by + 0.5) * b
  } else {
    cx <- stats::median(d$x_px); cy <- stats::median(d$y_px)
  }
  x0 <- as.integer(max(0, min(slides$sizeX[k] - win_px, cx - win_px / 2)))
  y0 <- as.integer(max(0, min(slides$sizeY[k] - win_px, cy - win_px / 2)))
  w <- min(win_px, slides$sizeX[k]); h <- min(win_px, slides$sizeY[k])

  a0 <- .open_level(tiff, 0L)
  R <- .stretch(.read_plane(a0, cm$PDGF, x0, y0, w, h), 0.05, 0.998)
  G <- .stretch(.read_plane(a0, cm$panCK, x0, y0, w, h), 0.05, 0.998)
  B <- .stretch(.read_plane(a0, cm$DAPI_AR, x0, y0, w, h), 0.02, 0.995)
  M <- .stretch(.read_plane(a0, cm$DLL3, x0, y0, w, h), 0.05, 0.998)
  ## DLL3 added as magenta (R+B) so all four fit in one RGB image
  rgb_arr <- array(c(pmin(R + 0.6 * M, 1), G, pmin(0.7 * B + 0.6 * M, 1)), dim = c(h, w, 3))
  inset <- grDevices::as.raster(rgb_arr)

  di <- d |>
    filter(x_px >= x0, x_px < x0 + w, y_px >= y0, y_px < y0 + h, pop5 != "Other") |>
    mutate(xi = (x_px - x0) * PX_UM, yi = (y_px - y0) * PX_UM)

  pC <- ggplot() +
    annotation_raster(inset, xmin = 0, xmax = w * PX_UM, ymin = -h * PX_UM, ymax = 0) +
    geom_point(data = di, aes(xi, -yi, colour = pop5), shape = 1, size = 2.4, stroke = 0.8) +
    scale_colour_manual(NULL, values = POP_COLS, drop = FALSE) +
    coord_equal(xlim = c(0, w * PX_UM), ylim = c(-h * PX_UM, 0), expand = FALSE) +
    scale_x_continuous("um") + scale_y_continuous("um", labels = function(v) -v) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          plot.title = element_text(face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size = 3, stroke = 1), nrow = 1)) +
    labs(title = sprintf("C  ASPC hotspot, full resolution (x=%.2f, y=%.2f mm)",
                         x0 * PX_UM / 1000, y0 * PX_UM / 1000),
         subtitle = "red PDGFa, green panCK, magenta DLL3, blue DAPI\nyellow/white = autofluorescence (bright in all channels); circles = cell calls")

  ## ---------------- D: neighbourhood composition ---------------------------
  nd <- if (!is.null(nbr)) nbr |> filter(slide_id == sid,
                                          neighbour %in% c("ASPC", "EPI_CA", "NE_CA", "Multi", "Other")) else NULL
  if (!is.null(nd) && nrow(nd)) {
    ndl <- rbind(
      data.frame(neighbour = nd$neighbour, what = "around ASPC (observed)", frac = nd$obs_frac),
      data.frame(neighbour = nd$neighbour, what = "whole biopsy (expected)", frac = nd$exp_frac)
    )
    ndl$neighbour <- factor(ndl$neighbour, levels = c("ASPC", "EPI_CA", "NE_CA", "Multi", "Other"))
    pD <- ggplot(ndl, aes(neighbour, frac, fill = what)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.7) +
      geom_text(data = nd, inherit.aes = FALSE,
                aes(x = neighbour, y = pmax(obs_frac, exp_frac) + 0.03,
                    label = sprintf("%+.2f", log2_enrich)), size = 3) +
      scale_fill_manual(NULL, values = c("around ASPC (observed)" = "#E8384F",
                                         "whole biopsy (expected)" = "grey60")) +
      scale_y_continuous("fraction of cells", limits = c(0, NA)) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
      labs(x = NULL, title = paste0("D  Composition within ", RADIUS_UM, " um of ASPC"),
           subtitle = sprintf("%d ASPC, median %g neighbours; label = log2 obs/exp",
                              nd$n_aspc[1], nd$median_neighbours[1]))
  } else {
    pD <- ggplot() + theme_void() + labs(title = "D  Composition around ASPC", subtitle = "not available")
  }

  fig <- (pA | pB) / (pC | pD) +
    plot_annotation(
      title = sprintf("%s  [%s]   %s cells | ASPC %d | EPI_CA %d | NE_CA %d | Multi %d",
                      sid, grp, format(nrow(d), big.mark = ","),
                      sum(d$pop5 == "ASPC"), sum(d$pop5 == "EPI_CA"),
                      sum(d$pop5 == "NE_CA"), sum(d$pop5 == "Multi")),
      subtitle = paste("Marker confidence (marker vs matched blank at the same threshold):",
                       paste(vapply(c("PDGF", "panCK", "FOLH1", "AR", "DLL3", "CEBPg"),
                                    function(m) .conf_of(sid, m), character(1)),
                             collapse = " | ")),
      theme = theme(plot.title = element_text(face = "bold", size = 15))
    )

  fig
}

###############################################################################
## RUN
###############################################################################

sids <- intersect(names(GROUP_OF), unique(cells$slide_id))
figs <- list()
for (sid in sids) {
  message("05: ", sid)
  fig <- tryCatch(make_slide_figure(sid), error = function(e) {
    warning(sid, ": ", conditionMessage(e), call. = FALSE); NULL
  })
  if (is.null(fig)) next
  figs[[sid]] <- fig
  ggsave(file.path(PLOTDIR, paste0("05_biopsy_", sid, ".png")), fig,
         width = 16, height = 15, dpi = 130, bg = "white")
}

if (length(figs)) {
  grDevices::cairo_pdf(file.path(PLOTDIR, "05_biopsy_atlas.pdf"), width = 16, height = 15, onefile = TRUE)
  for (f in figs) print(f)
  invisible(dev.off())
}

message("05_visualize.R complete -> ", PLOTDIR)
