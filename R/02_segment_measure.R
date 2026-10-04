###############################################################################
# 02_segment_measure.R  v12
#
# Tiled segmentation + corrected per-compartment marker feature extraction.
#
# CHANGES IN v12 (review against the COMET run logs, 2026-09-23):
#   1. Blank model. Each marker is now regressed on ALL same-fluorophore blank
#      planes acquired BEFORE it (PRIOR_BLANKS), not a single NC plane.
#      v11 corrected CEBPg (cycle 2) with NC3_Cy5, which is acquired AFTER the
#      CEBPg elution and therefore contains un-eluted CEBPg; the affine fit then
#      subtracted part of the true CEBPg signal. For cycle 6/7 markers the
#      prior set contains both the species-matched secondary-only blank
#      (NC5 for AR/DLL3, NC3 for PDGFa/panCK) and the most recent blank.
#      Negative slopes are dropped and the fit repeated (non-negative model).
#   2. Blank pseudo-marker (false-positive null) is the species/fluor-matched
#      secondary-only control from the run log, corrected against the blanks
#      acquired before it (BLANK_CTRL_OF / prior set). It is used in 03 to
#      anchor the marker threshold.
#   3. Zero padding. COMET writes 0 for FOVs that were not acquired
#      ("Only acquire FOVs with tissue: Yes"). v11 ran Otsu over the padded
#      tile, so the threshold split "acquired vs padding" instead of
#      "nucleus vs background"; on thin cores (PM423: 4,029 cells from 8
#      tiles) nearly every nucleus merged into one over-size object and was
#      discarded. Tissue and nuclear thresholds now ignore unacquired pixels.
#   4. Tile grid now includes the right/bottom partial tiles (v11 dropped up
#      to 2,047 px per edge), and tiles are pre-screened on the pyramid so
#      unacquired / empty tiles are never read at full resolution.
#   5. Per-slide mode for parallel runs: set env PIPE_SLIDES="PM423,PM12_B3"
#      to segment only those slides (writes cells_raw_<sid>.rds only);
#      set env PIPE_MERGE=1 to merge per-slide files and write QC outputs.
#      PIPE_SHARD="i/n" additionally splits one slide's tiles across n
#      processes (cells_raw_<sid>_shard<i>of<n>.rds; merged automatically).
#
# CHANGE IN v11.3:
#   Registration QC now uses robust informative-window cross-cycle DAPI
#   correlation, treating empty or flat windows as unknown rather than failed.
#
# CHANGE IN v11.2:
#   PM0T32 has been re-run on the standard 8-cycle protocol WITH secondary-only
#   negative controls (AR/DLL3 = Cycle 6). It is therefore no longer a protocol
#   variant and is included by default, exactly like every other slide.
#   The INCLUDE_* exclusion switch has been removed.
#
# v11 image-analysis model:
#   1. Nuclear markers:
#        CEBPg, AR measured in eroded nuclear masks.
#
#   2. Membrane / cytoplasmic markers:
#        PDGF, panCK, FOLH1, DLL3 measured in local perinuclear rings.
#        Ring radii: 4, 6, 8 µm.
#        Primary ring radius: 6 µm.
#
#   3. Pixel-level scaled matched-blank subtraction:
#        marker - (blank intercept + blank slope * matched secondary-only
#        negative control), fit on valid non-AF pixels after upper-tail trimming.
#
#   4. AF puncta exclusion:
#        pixels bright in both marker and matched blank are masked,
#        dilated, and excluded before cell aggregation.
#
#   5. Blank false-positive control:
#        matched secondary-only negative control is treated as a pseudo-marker,
#        corrected against fluorophore-matched AF channel.
#
#   6. Radius sensitivity:
#        Primary 6 µm features are named *_ring_*.
#        4 and 8 µm sensitivity features are named *_ring_r4um_* and *_ring_r8um_*.
#
# Main output columns used by 03_classify.R:
#
#   Nuclear marker example:
#        CEBPg_nuc_sub
#        CEBPg_nuc_frac_pos
#        CEBPg_nuc_af_frac
#        CEBPg_nuc_blank_sub
#        CEBPg_nuc_blank_frac_pos
#        CEBPg_nuc_blank_af_frac
#
#   Primary 6 µm ring marker example:
#        PDGF_ring_sub
#        PDGF_ring_frac_pos
#        PDGF_ring_af_frac
#        PDGF_ring_blank_sub
#        PDGF_ring_blank_frac_pos
#        PDGF_ring_blank_af_frac
#
#   Sensitivity ring examples:
#        PDGF_ring_r4um_sub
#        PDGF_ring_r4um_frac_pos
#        PDGF_ring_r8um_sub
#        PDGF_ring_r8um_frac_pos
###############################################################################

###############################################################################
## USER / RUN DEFAULTS
###############################################################################

if (!exists("OUT")) OUT <- "results_ASPC_NEPC"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

if (!exists("ROOT")) ROOT <- Sys.getenv("ED7_ROOT", "/Volumes/T9/Case Studies")

if (!exists("PROCEED_DESPITE_AF")) PROCEED_DESPITE_AF <- TRUE

## These can be overridden before source("02_segment_measure.R").
if (!exists("TILE_SIZE")) TILE_SIZE <- 2048L
if (!exists("MAX_TILES")) MAX_TILES <- Inf
if (!exists("VERBOSE_EVERY")) VERBOSE_EVERY <- 25L

## v11 classifier does not need old raw marker means.
if (!exists("KEEP_RAW_MEANS")) KEEP_RAW_MEANS <- FALSE

###############################################################################
## SLIDE MANIFEST
###############################################################################
## All six slides now share the 8-cycle protocol with secondary-only negative
## controls, so no slide-specific handling is required here.

MANIFEST <- data.frame(
  slide_id = c("PM0T32", "PM677_A8", "PM12_B3",
               "PM1024_B", "PM423", "PM159T5"),
  group    = c("NEPC", "NEPC", "NEPC",
               "ADENO", "ADENO", "ADENO"),
  tiff = file.path(ROOT, c(
    "20260826_155725_1_ybucNp_JJ_full panel_PM0 T32/20260826_155725_1_ybucNp_JJ_full panel_PM0 T32.tiff",
    "20260608_164617_1_ahXvmt_JJ_full panel_PM677 A8/20260608_164617_1_ahXvmt_JJ_full panel_PM677 A8.tiff",
    "20260608_164617_2_9ZdZ8V_JJ_full panel_PM12 B3/20260608_164617_2_9ZdZ8V_JJ_full panel_PM12 B3.tiff",
    "20260608_164617_4_CA2z2F_JJ_full panel_PM1024 B/20260608_164617_4_CA2z2F_JJ_full panel_PM1024 B.tiff",
    "20260625_131559_2_ZxKUqB_JJ_full panel_PM423/20260625_131559_2_ZxKUqB_JJ_full panel_PM423.tiff",
    "20260625_131600_3_VHgSnt_JJ_full panel_PM159T5/20260625_131600_3_VHgSnt_JJ_full panel_PM159T5.tiff")),
  stringsAsFactors = FALSE
)

resolve_tiff <- function(p) {
  if (grepl("\\.ome\\.tiff?$", p, ignore.case = TRUE)) return(p)
  alt <- sub("\\.tiff?$", ".ome.tiff", p, ignore.case = TRUE)
  if (file.exists(alt)) alt else p
}

MANIFEST$tiff_resolved <- vapply(
  MANIFEST$tiff,
  resolve_tiff,
  character(1),
  USE.NAMES = FALSE
)

cat("\nSlide manifest:\n")
for (i in seq_len(nrow(MANIFEST))) {
  ok <- file.exists(MANIFEST$tiff_resolved[i])
  swp <- MANIFEST$tiff_resolved[i] != MANIFEST$tiff[i]

  cat(sprintf(
    "  %-9s %-6s %s%s  %s\n",
    MANIFEST$slide_id[i],
    MANIFEST$group[i],
    if (ok) "[found]" else "[MISSING]",
    if (swp) " [.ome.tiff]" else "",
    basename(MANIFEST$tiff_resolved[i])
  ))
}

if (any(!file.exists(MANIFEST$tiff_resolved))) {
  stop("Missing TIFF(s) -- is /Volumes/T9 mounted?", call. = FALSE)
}

###############################################################################
## LOAD CHANNEL MAP
###############################################################################

channels_file <- file.path(OUT, "channels.rds")

if (!file.exists(channels_file)) {
  stop("Missing channels.rds. Re-run 01_channels.R first.", call. = FALSE)
}

z <- readRDS(channels_file)

slides <- z$slides
info   <- z$info
PX_UM  <- if (!is.null(z$px_um)) z$px_um else 12500 / 44643

unknown <- setdiff(MANIFEST$slide_id, slides$slide_id)
if (length(unknown)) {
  stop(
    "No channel map for: ",
    paste(unknown, collapse = ", "),
    "\nRe-run 01_channels.R (PM0T32 was re-acquired 20260826 and must be re-mapped).",
    call. = FALSE
  )
}

slides <- slides[match(MANIFEST$slide_id, slides$slide_id), , drop = FALSE]
slides$tiff  <- MANIFEST$tiff_resolved
slides$group <- MANIFEST$group
info <- info[MANIFEST$slide_id]

slides$sizeX <- as.integer(unlist(slides$sizeX))
slides$sizeY <- as.integer(unlist(slides$sizeY))

for (cl in c("tier2", "tier3", "has_negctrl")) {
  if (!cl %in% names(slides)) slides[[cl]] <- TRUE
}

for (cl in c("run", "comet_sw")) {
  if (!cl %in% names(slides)) slides[[cl]] <- NA
}

## All slides must now carry negative controls; the correction model requires it.
if (any(!slides$has_negctrl, na.rm = TRUE)) {
  stop(
    "Slide(s) without secondary-only negative controls: ",
    paste(slides$slide_id[!slides$has_negctrl], collapse = ", "),
    "\nv11 matched-blank subtraction requires NC cycles on every slide.",
    call. = FALSE
  )
}

cat("\nGroups:\n")
print(table(slides$group))

###############################################################################
## PARAMETERS
###############################################################################

try(mem.maxVSize(vsize = 64 * 1024), silent = TRUE)

MIN_AREA_UM2 <- 8
MAX_AREA_UM2 <- 400

MIN_AREA_PX <- MIN_AREA_UM2 / PX_UM^2
MAX_AREA_PX <- MAX_AREA_UM2 / PX_UM^2

QC_WIN <- 512L
if (!exists("DAPI_MIN_SD")) DAPI_MIN_SD <- 10
if (!exists("DAPI_MIN_DYNAMIC_RANGE")) DAPI_MIN_DYNAMIC_RANGE <- 30
if (!exists("REG_LOW_COR_THRESHOLD")) REG_LOW_COR_THRESHOLD <- 0.85

## Ring radius sensitivity.
if (!exists("RING_UM_PRIMARY")) RING_UM_PRIMARY <- 6
if (!exists("RING_UM_SET")) RING_UM_SET <- c(4, 6, 8)

.odd <- function(x) {
  x <- as.integer(round(x))
  x + as.integer(x %% 2L == 0L)
}

RING_R_PX_SET <- setNames(
  pmax(3L, .odd(RING_UM_SET / PX_UM)),
  paste0("r", gsub("\\.", "p", as.character(RING_UM_SET)), "um")
)

## Pixel sampling for expensive quantiles/backgrounds.
if (!exists("PIX_SAMPLE_N")) PIX_SAMPLE_N <- 200000L

## AF puncta parameters.
if (!exists("AF_Q")) AF_Q <- 0.995
if (!exists("AF_DILATE_R")) AF_DILATE_R <- 5L
if (!exists("MIN_PIX_MAD")) MIN_PIX_MAD <- 25

PY_ENV <- "r-reticulate"

suppressPackageStartupMessages({
  library(EBImage)
  library(MASS)
  library(dplyr)
  library(purrr)
  library(tibble)
})

###############################################################################
## MARKERS, CONTROLS, COMPARTMENTS
###############################################################################

MARKERS_CORE <- c("CEBPg", "FOLH1", "DLL3", "AR", "panCK", "PDGF")

CTRL_CHANS <- c(
  "AF_Cy5", "AF_TRITC",
  "NC3_Cy5", "NC3_TRITC",
  "NC5_Cy5", "NC5_TRITC",
  "NC7_Cy5", "NC7_TRITC"
)

MARKERS <- c(MARKERS_CORE, CTRL_CHANS)

## From the COMET run log (all 6 imaged slides, identical plane order):
##   cycle 1  AF            (DAPI, TRITC, Cy5)             -> AF_*
##   cycle 2  CEBPg  Cy5    rabbit, anti-Rabbit-647
##   cycle 3  NC            anti-Rabbit-647 / anti-Mouse-555 -> NC3_*
##   cycle 4  FOLH1  Cy5    mouse,  anti-Mouse-647
##   cycle 5  NC            anti-Mouse-647 / anti-Rabbit-555 -> NC5_*
##   cycle 6  AR  Cy5 (mouse, anti-Mouse-647), DLL3 TRITC (rabbit, anti-Rabbit-555)
##   cycle 7  PDGFa Cy5 (rabbit, anti-Rabbit-647), panCK TRITC (mouse, anti-Mouse-555)
##   cycle 8  elution only, re-imaged                        -> NC7_* (terminal)
##
## PRIOR_BLANKS: same-fluorophore blanks acquired before the marker. A blank
## acquired after the marker (e.g. NC3 for CEBPg) contains un-eluted marker
## and must not be regressed out of it.
PRIOR_BLANKS <- list(
  CEBPg = c("AF_Cy5"),
  FOLH1 = c("AF_Cy5", "NC3_Cy5"),
  AR    = c("AF_Cy5", "NC3_Cy5", "NC5_Cy5"),
  DLL3  = c("AF_TRITC", "NC3_TRITC", "NC5_TRITC"),
  PDGF  = c("AF_Cy5", "NC3_Cy5", "NC5_Cy5"),
  panCK = c("AF_TRITC", "NC3_TRITC", "NC5_TRITC")
)

## Species/fluorophore-matched secondary-only control per marker: the
## pseudo-marker used as the false-positive null.
CTRL_OF <- c(
  CEBPg = "NC3_Cy5",
  FOLH1 = "NC5_Cy5",
  AR    = "NC5_Cy5",
  DLL3  = "NC5_TRITC",
  PDGF  = "NC3_Cy5",
  panCK = "NC3_TRITC"
)

MARKER_COMPARTMENT <- c(
  CEBPg = "nuc",
  AR    = "nuc",
  PDGF  = "ring",
  panCK = "ring",
  FOLH1 = "ring",
  DLL3  = "ring"
)

## The other secondary-only control in the same fluorophore. Both NC3 and NC5
## carry anti-Rabbit and anti-Mouse secondaries; blank images can contain
## cycle-specific artifacts (aggregates, debris) absent from the marker image,
## so 03 anchors thresholds to the cleaner of the two nulls.
ALT_CTRL_OF <- c(
  CEBPg = "NC5_Cy5",
  FOLH1 = "NC3_Cy5",
  AR    = "NC3_Cy5",
  DLL3  = "NC3_TRITC",
  PDGF  = "NC5_Cy5",
  panCK = "NC5_TRITC"
)

## For blank false-positive control: the secondary-only control is treated as
## a pseudo-marker and corrected against the blanks acquired before it.
BLANK_REF_OF <- list(
  NC3_Cy5   = c("AF_Cy5"),
  NC5_Cy5   = c("AF_Cy5", "NC3_Cy5"),
  NC3_TRITC = c("AF_TRITC"),
  NC5_TRITC = c("AF_TRITC", "NC3_TRITC")
)

###############################################################################
## CHANNEL MAP GUARD
###############################################################################

bad <- character(0)

for (sid in names(info)) {
  cm <- info[[sid]]$chmap
  da <- as.integer(cm$DAPI_all)

  if (!identical(da, seq_along(da))) {
    bad <- c(bad, paste0(sid, ": DAPI_all = ", paste(da, collapse = ",")))
  }

  if (!is.null(cm$verify) && isFALSE(cm$verify$agrees)) {
    bad <- c(bad, paste0(sid, ": failed empirical DAPI check in 01"))
  }

  for (m in MARKERS) {
    v <- cm[[m]]

    if (is.null(v) || !length(v) || is.na(v[1])) {
      bad <- c(bad, paste0(sid, ": channel '", m, "' is NA"))
    } else if (v[1] %in% da) {
      bad <- c(bad, paste0(sid, ": channel '", m, "' collides with DAPI plane"))
    }
  }

  for (m in names(CTRL_OF)) {
    for (ctl in c(CTRL_OF[[m]], PRIOR_BLANKS[[m]])) {
      if (is.null(cm[[ctl]]) || is.na(cm[[ctl]][1])) {
        bad <- c(bad, paste0(sid, ": missing blank for ", m, " -> ", ctl))
      }
    }
  }

  for (ctl in unique(CTRL_OF)) {
    for (af in BLANK_REF_OF[[ctl]]) {
      if (is.null(af) || is.null(cm[[af]]) || is.na(cm[[af]][1])) {
        bad <- c(bad, paste0(sid, ": missing reference for ", ctl, " -> ", af))
      }
    }
  }
}

if (length(bad)) {
  stop(
    "channels.rds is stale, unverified, or incomplete:\n  ",
    paste(unique(bad), collapse = "\n  "),
    "\nRe-run 01_channels.R before segmenting.",
    call. = FALSE
  )
}

message("chmap guard: OK")

###############################################################################
## IMAGE READER
###############################################################################

READ_MODE <- NA_character_

.setup_python <- function(auto_install = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)

  library(reticulate)

  if (nzchar(PY_ENV)) {
    try(use_virtualenv(PY_ENV, required = TRUE), silent = TRUE)
  }

  cfg <- tryCatch(py_config(), error = function(e) NULL)
  if (is.null(cfg)) return(NULL)

  message("Python: ", cfg$python)

  have <- function(m) {
    tryCatch({
      import(m)
      TRUE
    }, error = function(e) FALSE)
  }

  need <- c("tifffile", "numpy")[!vapply(c("tifffile", "numpy"), have, logical(1))]

  if (length(need) && auto_install) {
    message("Installing Python modules: ", paste(need, collapse = ", "))
    try(py_install(c(need, "imagecodecs")), silent = TRUE)
    need <- need[!vapply(need, have, logical(1))]
  }

  if (length(need)) return(NULL)

  if (!have("zarr") && auto_install) {
    message("Installing zarr ...")
    try(py_install("zarr>=3"), silent = TRUE)
  }

  list(
    tf = import("tifffile", convert = FALSE),
    np = import("numpy", convert = FALSE),
    py = import_builtins(),
    zarr = if (have("zarr")) import("zarr", convert = FALSE) else NULL
  )
}

.PY <- .setup_python()

.pynum <- function(x) as.numeric(unlist(reticulate::py_to_r(x)))

if (!is.null(.PY)) {
  .tf <- .PY$tf
  .np <- .PY$np
  .py <- .PY$py
  .zr <- .PY$zarr
  .cache <- new.env(parent = emptyenv())

  if (!is.null(.zr)) {
    get_store <- function(tiff) {
      k <- basename(tiff)

      if (!is.null(.cache[[k]])) return(.cache[[k]])

      st <- .tf$TiffFile(tiff)$series[0]$aszarr()

      g <- tryCatch(
        .zr$open(store = st, mode = "r"),
        error = function(e) .zr$open(st, mode = "r")
      )

      a <- tryCatch(g$`__getitem__`("0"), error = function(e) g)

      .cache[[k]] <- a
      a
    }

    img_read <- function(tiff, c_index, x0 = 0, y0 = 0, w = NULL, h = NULL) {
      a <- get_store(tiff)

      sl <- function(i, j) .py$slice(as.integer(i), as.integer(j))

      idx <- if (is.null(w)) {
        reticulate::tuple(
          as.integer(c_index - 1L),
          .py$slice(NULL),
          .py$slice(NULL)
        )
      } else {
        reticulate::tuple(
          as.integer(c_index - 1L),
          sl(y0, y0 + h),
          sl(x0, x0 + w)
        )
      }

      m <- reticulate::py_to_r(.np$asarray(a$`__getitem__`(idx)))
      m <- t(m)
      storage.mode(m) <- "double"
      m
    }

    READ_MODE <- "tifffile+zarr"

  } else {
    CACHE_MAX <- 2L

    get_plane <- function(tiff, c_index) {
      k <- paste0(basename(tiff), "#", c_index)

      if (!is.null(.cache[[k]])) return(.cache[[k]])

      pk <- paste0("pages#", basename(tiff))

      if (is.null(.cache[[pk]])) {
        .cache[[pk]] <- .tf$TiffFile(tiff)$series[0]$pages
      }

      pg <- .cache[[pk]]$`__getitem__`(as.integer(c_index - 1L))
      m <- reticulate::py_to_r(.np$asarray(pg$asarray()))

      keys <- setdiff(ls(.cache), grep("^pages#", ls(.cache), value = TRUE))

      if (length(keys) >= CACHE_MAX) {
        rm(list = keys[1], envir = .cache)
        gc(FALSE)
      }

      assign(k, m, envir = .cache)
      m
    }

    img_read <- function(tiff, c_index, x0 = 0, y0 = 0, w = NULL, h = NULL) {
      m <- get_plane(tiff, c_index)

      if (!is.null(w)) {
        m <- m[
          (y0 + 1):min(y0 + h, nrow(m)),
          (x0 + 1):min(x0 + w, ncol(m)),
          drop = FALSE
        ]
      }

      m <- t(m)
      storage.mode(m) <- "double"
      m
    }

    READ_MODE <- "tifffile+pagecache"
  }

} else if (requireNamespace("RBioFormats", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(rJava)
    library(RBioFormats)
  })

  img_read <- function(tiff, c_index, x0 = 0, y0 = 0, w = NULL, h = NULL) {
    sub <- list(C = as.integer(c_index))

    if (!is.null(w)) {
      sub <- c(
        sub,
        list(
          x = as.integer(x0),
          y = as.integer(y0),
          w = as.integer(w),
          h = as.integer(h)
        )
      )
    }

    a <- as.array(RBioFormats::read.image(
      tiff,
      series = 1,
      subset = sub,
      normalize = FALSE,
      filter.metadata = TRUE
    ))

    if (length(dim(a)) > 2) a <- a[, , 1]
    storage.mode(a) <- "double"
    a
  }

  READ_MODE <- "bioformats"

} else {
  stop(
    "No usable image reader. Install reticulate/tifffile/zarr or RBioFormats.",
    call. = FALSE
  )
}

message("READ MODE: ", READ_MODE)
message("PX_UM: ", signif(PX_UM, 4), " µm/px")
message("TILE_SIZE: ", TILE_SIZE)
message("Area filter: ", round(MIN_AREA_PX), "-", round(MAX_AREA_PX), " px")
message("Ring radii: ", paste(RING_UM_SET, collapse = ", "), " µm")
message("Primary ring: ", RING_UM_PRIMARY, " µm")
message("Ring px: ", paste(names(RING_R_PX_SET), RING_R_PX_SET, sep = "=", collapse = ", "))

###############################################################################
## STORE GEOMETRY CHECK
###############################################################################

if (grepl("^tifffile\\+zarr", READ_MODE)) {
  dimtab <- do.call(rbind, lapply(seq_len(nrow(slides)), function(k) {
    a <- get_store(slides$tiff[k])
    d <- .pynum(a$shape)
    ch <- .pynum(a$chunks)

    data.frame(
      slide_id = slides$slide_id[k],
      group = slides$group[k],
      planes = d[1],
      sizeX = d[3],
      sizeY = d[2],
      chunk = paste(ch[2], ch[3], sep = "x"),
      aligned = TILE_SIZE %% ch[2] == 0 && TILE_SIZE %% ch[3] == 0,
      n_dapi = length(info[[slides$slide_id[k]]]$chmap$DAPI_all),
      row.names = NULL
    )
  }))

  cat("\nStore geometry:\n")
  print(dimtab, row.names = FALSE)

  bad_dim <- which(dimtab$sizeX != slides$sizeX | dimtab$sizeY != slides$sizeY)

  if (length(bad_dim)) {
    warning(
      "Dimensions differ from channels.rds for: ",
      paste(dimtab$slide_id[bad_dim], collapse = ", "),
      " -- using on-disk values.",
      call. = FALSE
    )

    slides$sizeX <- as.integer(dimtab$sizeX)
    slides$sizeY <- as.integer(dimtab$sizeY)
  }
}

###############################################################################
## HELPERS
###############################################################################

.qsafe <- function(x, p, default = NA_real_, n_sample = PIX_SAMPLE_N) {
  x <- x[is.finite(x)]

  if (!length(x)) return(default)

  if (length(x) > n_sample) {
    set.seed(1L)
    x <- sample(x, n_sample)
  }

  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
}

.is_informative_dapi <- function(v,
                                 min_sd = DAPI_MIN_SD,
                                 min_range = DAPI_MIN_DYNAMIC_RANGE) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (length(v) < 1000L) return(FALSE)
  s <- stats::sd(v, na.rm = TRUE)
  r <- as.numeric(stats::quantile(v, 0.99, na.rm = TRUE, names = FALSE) -
                    stats::median(v, na.rm = TRUE))
  is.finite(s) && is.finite(r) && s >= min_sd && r >= min_range
}

.safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 1000L) return(NA_real_)
  if (stats::sd(a[ok], na.rm = TRUE) <= 0) return(NA_real_)
  if (stats::sd(b[ok], na.rm = TRUE) <= 0) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok]))
}

.dapi_registration_qc <- function(tiff, cm, x0, y0, w, h) {
  out <- list(
    reg_cor = NA_real_,
    reg_cor_min = NA_real_,
    reg_cor_median = NA_real_,
    reg_cor_p05 = NA_real_,
    reg_cor_n = 0L,
    reg_cor_n_low = NA_integer_,
    reg_cor_low_frac = NA_real_,
    reg_qc_informative = FALSE,
    reg_ref_sd = NA_real_,
    reg_ref_range = NA_real_
  )

  if (length(cm$DAPI_all) <= 1L) return(out)

  qx <- x0 + max(0L, (w - QC_WIN) %/% 2L)
  qy <- y0 + max(0L, (h - QC_WIN) %/% 2L)

  ref <- as.vector(img_read(tiff, cm$DAPI_all[1], qx, qy, QC_WIN, QC_WIN))
  ref_ok <- ref[is.finite(ref)]

  if (length(ref_ok) >= 1000L) {
    out$reg_ref_sd <- stats::sd(ref_ok, na.rm = TRUE)
    out$reg_ref_range <- as.numeric(stats::quantile(ref_ok, 0.99, na.rm = TRUE, names = FALSE) -
                                      stats::median(ref_ok, na.rm = TRUE))
  }

  if (!.is_informative_dapi(ref)) return(out)

  cors <- vapply(cm$DAPI_all[-1], function(i) {
    v <- as.vector(img_read(tiff, i, qx, qy, QC_WIN, QC_WIN))
    if (!.is_informative_dapi(v)) return(NA_real_)
    .safe_cor(ref, v)
  }, numeric(1))

  cors <- cors[is.finite(cors)]
  out$reg_cor_n <- length(cors)
  if (!length(cors)) return(out)

  out$reg_qc_informative <- TRUE
  out$reg_cor_min <- suppressWarnings(min(cors, na.rm = TRUE))
  out$reg_cor_median <- suppressWarnings(stats::median(cors, na.rm = TRUE))
  out$reg_cor_p05 <- suppressWarnings(as.numeric(stats::quantile(cors, 0.05, na.rm = TRUE, names = FALSE)))
  out$reg_cor_n_low <- sum(cors < REG_LOW_COR_THRESHOLD, na.rm = TRUE)
  out$reg_cor_low_frac <- mean(cors < REG_LOW_COR_THRESHOLD, na.rm = TRUE)
  out$reg_cor <- out$reg_cor_median
  out
}

.weighted_sum_by_label <- function(values, labels, nmax) {
  out <- rep(0, nmax)

  if (!length(values)) return(out)

  rs <- rowsum(values, group = labels, reorder = FALSE)
  ii <- suppressWarnings(as.integer(rownames(rs)))
  keep <- is.finite(ii) & ii > 0L & ii <= nmax

  out[ii[keep]] <- as.numeric(rs[keep, 1])
  out
}

.mean_by_label_valid <- function(img, lab, valid, nmax) {
  l <- as.integer(as.vector(lab))
  v <- as.numeric(as.vector(img))

  ok <- as.vector(valid) &
        is.finite(v) &
        l > 0L &
        l <= nmax

  if (!any(ok)) return(rep(NA_real_, nmax))

  ## Important:
  ## base::tabulate() does not support weights=.
  ## Use rowsum() for weighted numerator and base::tabulate() for denominator.
  num <- .weighted_sum_by_label(v[ok], l[ok], nmax)
  den <- base::tabulate(l[ok], nbins = nmax)

  out <- num / pmax(den, 1)
  out[den == 0] <- NA_real_
  out
}

.frac_by_label_valid <- function(mask, lab, valid, nmax) {
  l <- as.integer(as.vector(lab))
  m <- as.logical(as.vector(mask))

  ok <- as.vector(valid) &
        l > 0L &
        l <= nmax

  if (!any(ok)) return(rep(NA_real_, nmax))

  den <- base::tabulate(l[ok], nbins = nmax)

  ok_pos <- ok & m
  if (any(ok_pos)) {
    num <- base::tabulate(l[ok_pos], nbins = nmax)
  } else {
    num <- rep(0, nmax)
  }

  out <- num / pmax(den, 1)
  out[den == 0] <- NA_real_
  out
}

## Otsu threshold on a numeric vector (EBImage::otsu works on whole images,
## which would include unacquired zero padding).
.otsu_vec <- function(v, n_sample = PIX_SAMPLE_N, levels = 256L) {
  v <- v[is.finite(v)]
  if (length(v) < 100L) return(NA_real_)
  if (length(v) > n_sample) {
    set.seed(1L)
    v <- sample(v, n_sample)
  }
  r <- range(v)
  if (r[2] <= r[1]) return(r[1])
  br <- seq(r[1], r[2], length.out = levels + 1L)
  h <- as.numeric(tabulate(findInterval(v, br, rightmost.closed = TRUE), nbins = levels))
  mids <- (br[-1] + br[-length(br)]) / 2
  w0 <- cumsum(h)
  w1 <- sum(h) - w0
  m0 <- cumsum(h * mids) / pmax(w0, 1)
  m1 <- (sum(h * mids) - cumsum(h * mids)) / pmax(w1, 1)
  bcv <- w0 * w1 * (m0 - m1)^2
  mids[which.max(bcv)]
}

## `blank_imgs` is a list of same-fluorophore blank planes (>= 1).
.make_af_mask <- function(marker_img, blank_imgs, valid_pix) {
  qm <- .qsafe(marker_img[valid_pix], AF_Q)

  blank_hi <- Reduce(`|`, lapply(blank_imgs, function(b) {
    qb <- .qsafe(b[valid_pix], AF_Q)
    is.finite(b) & b >= qb
  }))

  af <- is.finite(marker_img) &
        valid_pix &
        marker_img >= qm &
        blank_hi

  af <- EBImage::dilate(
    EBImage::Image(af),
    EBImage::makeBrush(AF_DILATE_R, "disc")
  ) > 0

  as.matrix(af)
}

## Robust non-negative multi-blank fit:
##   marker ~ intercept + sum_k slope_k * blank_k
## Blanks with a negative slope are dropped and the model refit, so a blank
## can only ever be subtracted, never added.
.fit_blank_affine <- function(marker_img, blank_imgs, valid_pix, af) {
  K <- length(blank_imgs)
  out <- list(intercept = 0, slope = rep(0, K), r2 = NA_real_)

  ok <- valid_pix & !af & is.finite(marker_img)
  for (b in blank_imgs) ok <- ok & is.finite(b)
  if (!any(ok)) return(out)

  marker <- as.numeric(marker_img[ok])
  B <- do.call(cbind, lapply(blank_imgs, function(b) as.numeric(b[ok])))

  keep <- marker <= .qsafe(marker, 0.99, default = Inf)
  for (k in seq_len(K)) keep <- keep & B[, k] <= .qsafe(B[, k], 0.99, default = Inf)

  marker <- marker[keep]
  B <- B[keep, , drop = FALSE]

  if (length(marker) > PIX_SAMPLE_N) {
    set.seed(1L)
    idx <- sample.int(length(marker), PIX_SAMPLE_N)
    marker <- marker[idx]
    B <- B[idx, , drop = FALSE]
  }

  if (length(marker) < 10L) return(out)

  active <- which(apply(B, 2, function(z) stats::sd(z) > 0))
  cf <- NULL

  while (length(active)) {
    fit_df <- data.frame(marker = marker, B[, active, drop = FALSE])
    names(fit_df) <- c("marker", paste0("b", active))
    fml <- stats::as.formula(paste("marker ~", paste(names(fit_df)[-1], collapse = " + ")))

    fit <- try(suppressWarnings(MASS::rlm(fml, data = fit_df, maxit = 50)), silent = TRUE)
    if (inherits(fit, "try-error")) fit <- try(stats::lm(fml, data = fit_df), silent = TRUE)
    if (inherits(fit, "try-error")) return(out)

    cf <- stats::coef(fit)
    if (any(!is.finite(cf))) return(out)

    neg <- which(cf[-1] < 0)
    if (!length(neg)) break
    active <- active[-neg]
    cf <- NULL
  }

  if (is.null(cf)) {
    out$intercept <- stats::median(marker)
    out$r2 <- 0
    return(out)
  }

  slope <- rep(0, K)
  slope[active] <- as.numeric(cf[-1])

  fitted <- as.numeric(cf[1] + B %*% slope)
  sst <- sum((marker - mean(marker))^2)
  sse <- sum((marker - fitted)^2)

  out$intercept <- as.numeric(cf[1])
  out$slope <- slope
  if (is.finite(sst) && sst > 0) out$r2 <- max(min(1 - sse / sst, 1), 0)
  out
}

.correct_marker_pixels <- function(marker_img, blank_imgs, valid_pix) {
  if (!is.list(blank_imgs)) blank_imgs <- list(blank_imgs)

  af <- .make_af_mask(marker_img, blank_imgs, valid_pix)
  blank_fit <- .fit_blank_affine(marker_img, blank_imgs, valid_pix, af)

  pred <- blank_fit$intercept
  for (k in seq_along(blank_imgs)) {
    if (blank_fit$slope[k] != 0) pred <- pred + blank_fit$slope[k] * blank_imgs[[k]]
  }
  resid <- marker_img - pred

  bg_pix <- valid_pix &
            !af &
            is.finite(resid)

  vals <- resid[bg_pix]
  vals <- vals[is.finite(vals)]

  if (length(vals) > PIX_SAMPLE_N) {
    set.seed(1L)
    vals <- sample(vals, PIX_SAMPLE_N)
  }

  bg <- median(vals, na.rm = TRUE)
  sc <- mad(vals, constant = 1.4826, na.rm = TRUE)

  if (!is.finite(bg)) bg <- 0
  if (!is.finite(sc) || sc <= 0) sc <- sd(vals, na.rm = TRUE)
  if (!is.finite(sc) || sc <= 0) sc <- MIN_PIX_MAD

  pix_thr <- bg + max(MIN_PIX_MAD, 3 * sc)

  sub <- pmax(resid - bg, 0)
  sub[af] <- NA_real_

  pos_pix <- resid > pix_thr &
             !af &
             valid_pix

  list(
    sub = sub,
    pos_pix = pos_pix,
    af = af,
    bg = bg,
    pix_thr = pix_thr,
    noise = sc,
    blank_intercept = blank_fit$intercept,
    ## summed slope keeps the v11 column meaning (1 = unit gain)
    blank_slope = sum(blank_fit$slope),
    blank_r2 = blank_fit$r2
  )
}

.aggregate_corrected_features <- function(cc, lab, valid_lab, nmax) {
  valid_comp <- valid_lab & !cc$af

  data.frame(
    sub      = .mean_by_label_valid(cc$sub,     lab, valid_comp, nmax),
    frac_pos = .frac_by_label_valid(cc$pos_pix, lab, valid_comp, nmax),
    af_frac  = .frac_by_label_valid(cc$af,      lab, valid_lab,  nmax),
    bg       = rep(cc$bg,      nmax),
    pix_thr  = rep(cc$pix_thr, nmax),
    noise    = rep(cc$noise,   nmax),
    blank_intercept = rep(cc$blank_intercept, nmax),
    blank_slope     = rep(cc$blank_slope,     nmax),
    blank_r2        = rep(cc$blank_r2,        nmax),
    stringsAsFactors = FALSE
  )
}

.col_or_na <- function(df, name) {
  if (name %in% names(df)) df[[name]] else rep(NA_real_, nrow(df))
}

.median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  suppressWarnings(stats::median(x, na.rm = TRUE))
}

.min_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  suppressWarnings(min(x, na.rm = TRUE))
}

.quantile_or_na <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  suppressWarnings(as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE)))
}

.ring_tag <- function(r_um) {
  paste0("r", gsub("\\.", "p", as.character(r_um)), "um")
}

.feature_prefix <- function(marker, compartment, r_um = NULL, primary = RING_UM_PRIMARY) {
  if (compartment == "nuc") {
    return(paste0(marker, "_nuc"))
  }

  if (is.null(r_um) || identical(as.numeric(r_um), as.numeric(primary))) {
    return(paste0(marker, "_ring"))
  }

  paste0(marker, "_ring_", .ring_tag(r_um))
}

.make_ring_labels <- function(nuc, nuc_mask, tis, d_s) {
  max_px <- max(RING_R_PX_SET)

  max_band <- dilate(nuc_mask, makeBrush(max_px, "disc")) &
              !nuc_mask &
              tis

  local_mask <- nuc_mask | max_band

  ## One propagation at largest radius, then subset to 4 / 6 / 8 µm.
  cell_max <- propagate(d_s, seeds = nuc, mask = local_mask, lambda = 1e-4)

  rings <- lapply(seq_along(RING_UM_SET), function(i) {
    r_um <- RING_UM_SET[i]
    r_px <- RING_R_PX_SET[i]

    band <- dilate(nuc_mask, makeBrush(r_px, "disc")) &
            !nuc_mask &
            tis

    ring <- cell_max
    ring[!band] <- 0
    ring
  })

  names(rings) <- vapply(RING_UM_SET, .ring_tag, character(1))

  rm(cell_max, max_band, local_mask)
  rings
}

###############################################################################
## TILE PROCESSING
###############################################################################

USE_CELLPOSE <- FALSE
CP_DIAM <- 36

## v12 nuclear segmentation parameters (see process_tile).
if (!exists("NUC_LOCAL_W")) NUC_LOCAL_W <- 21L          # half-window px (~6 um)
if (!exists("NUC_LOCAL_OFFSET")) NUC_LOCAL_OFFSET <- 0.01
if (!exists("NUC_FLOOR_FRAC")) NUC_FLOOR_FRAC <- 0.5     # floor between glass and Otsu
if (!exists("TISSUE_UM")) TISSUE_UM <- 10

process_tile <- function(tiff, cm, x0, y0, w, h, sid) {

  ## Registration QC.
  ## Empty/flat DAPI windows are unknown, not failures. `reg_cor` is now the
  ## median informative-window cross-cycle DAPI correlation for compatibility.
  reg_qc <- .dapi_registration_qc(tiff, cm, x0, y0, w, h)

  ## DAPI for segmentation.
  dapi_plane <- if (!is.null(cm$DAPI_AR)) cm$DAPI_AR else cm$DAPI_all[1]
  dapi <- img_read(tiff, dapi_plane, x0, y0, w, h)

  mx <- max(dapi, na.rm = TRUE)

  if (!is.finite(mx) || mx <= 0) {
    rm(dapi)
    gc(FALSE)
    return(NULL)
  }

  ## Acquired-FOV mask: COMET writes exact zeros for FOVs it did not image.
  ## Erode slightly so blur halos at FOV borders do not enter thresholds.
  acq <- erode(Image(dapi > 0), makeBrush(5, "disc")) > 0

  if (sum(acq) < 2000) {
    rm(dapi, acq)
    gc(FALSE)
    return(NULL)
  }

  d <- Image(dapi / mx)
  d_s <- gblur(d, sigma = 1)

  ## Nuclear segmentation (v12).
  ## v11 used one global Otsu per tile on linear DAPI inside an Otsu tissue
  ## mask; dim nuclei (about half of them in dense tiles) fell below both and
  ## were lost. v12 thresholds log-DAPI locally (mean filter, ~2 nuclear
  ## diameters) above a permissive floor computed on ACQUIRED pixels only.
  l  <- log1p(pmax(dapi, 0))
  ls <- gblur(Image(l / max(l)), 1.5)
  rm(dapi, l)

  lv <- as.numeric(ls[acq])
  floor_hi <- .otsu_vec(lv)
  floor_lo <- .qsafe(lv, 0.05)
  rm(lv)

  if (!is.finite(floor_hi) || !is.finite(floor_lo)) {
    rm(d, d_s, ls, acq)
    gc(FALSE)
    return(NULL)
  }

  fg_floor <- floor_lo + NUC_FLOOR_FRAC * (floor_hi - floor_lo)

  if (USE_CELLPOSE) {
    cp  <- reticulate::import("cellpose")
    npy <- reticulate::import("numpy")
    mdl <- cp$models$Cellpose(gpu = FALSE, model_type = "nuclei")

    nuc <- Image(reticulate::py_to_r(
      mdl$eval(
        npy$array(as.matrix(d_s)),
        diameter = CP_DIAM,
        channels = c(0L, 0L)
      )[[1]]
    ))

  } else {
    th <- thresh(ls, w = NUC_LOCAL_W, h = NUC_LOCAL_W, offset = NUC_LOCAL_OFFSET) &
          (ls > fg_floor) & acq
    th <- fillHull(opening(th, makeBrush(5, "disc")))

    nuc <- watershed(distmap(th), tolerance = 1.5, ext = 4)
    rm(th)
  }
  rm(ls)

  sh0 <- computeFeatures.shape(nuc)

  if (is.null(sh0) || !nrow(sh0)) {
    rm(d, d_s, nuc, acq)
    gc(FALSE)
    return(NULL)
  }

  drop <- which(
    sh0[, "s.area"] < MIN_AREA_PX |
      sh0[, "s.area"] > MAX_AREA_PX
  )

  if (length(drop)) {
    nuc <- rmObjects(nuc, drop, reenumerate = TRUE)
  }

  nmax <- max(nuc)

  if (nmax < 10) {
    rm(d, d_s, nuc, acq)
    gc(FALSE)
    return(NULL)
  }

  ## Tissue = acquired pixels within TISSUE_UM of a retained nucleus. Rings
  ## (<= 8 um) are clipped to it, so they never extend into empty glass.
  tis <- dilate(nuc > 0, makeBrush(.odd(TISSUE_UM / PX_UM), "disc")) & acq

  ## Compartments:
  ##   - eroded nucleus for CEBPg / AR
  ##   - local perinuclear rings for PDGF / panCK / FOLH1 / DLL3
  nuc_mask <- nuc > 0
  nuc_ero <- nuc * as.integer(erode(nuc_mask, makeBrush(3, "disc")))

  rings <- .make_ring_labels(nuc, nuc_mask, tis, d_s)

  rm(tis, d_s, d, nuc_mask, acq)
  gc(FALSE)

  ## Corrected marker and blank-as-marker features.
  corr_list <- list()

  ## Blank planes are shared by several markers; read each once per tile.
  plane_cache <- new.env(parent = emptyenv())
  .plane <- function(nm) {
    if (is.null(plane_cache[[nm]])) {
      assign(nm, img_read(tiff, cm[[nm]][1], x0, y0, w, h), envir = plane_cache)
    }
    plane_cache[[nm]]
  }

  .all_finite <- function(lst) Reduce(`&`, lapply(lst, is.finite))

  union_ring <- Reduce(`|`, lapply(rings, function(x) x > 0))
  comp_valid <- list(nuc = nuc_ero > 0, ring = union_ring)

  ## Blank pseudo-marker corrections depend only on (blank, compartment), so
  ## they are computed once per tile and shared by all markers.
  blank_cc_cache <- list()
  .blank_cc <- function(ctl, cp) {
    key <- paste(ctl, cp)
    if (is.null(blank_cc_cache[[key]])) {
      blk <- .plane(ctl)
      ref <- lapply(BLANK_REF_OF[[ctl]], .plane)
      valid <- comp_valid[[cp]] & is.finite(blk) & .all_finite(ref)
      blank_cc_cache[[key]] <<- .correct_marker_pixels(blk, ref, valid)
    }
    blank_cc_cache[[key]]
  }

  .emit <- function(cc, lab, prefix, kind) {
    f <- .aggregate_corrected_features(cc = cc, lab = lab, valid_lab = lab > 0, nmax = nmax)
    colnames(f) <- paste0(prefix, kind, colnames(f))
    corr_list[[length(corr_list) + 1]] <<- f
  }

  for (m in MARKERS_CORE) {
    cp  <- MARKER_COMPARTMENT[[m]]

    mk     <- img_read(tiff, cm[[m]][1], x0, y0, w, h)
    priors <- lapply(PRIOR_BLANKS[[m]], .plane)

    valid_marker <- comp_valid[[cp]] & is.finite(mk) & .all_finite(priors)
    cc_marker <- .correct_marker_pixels(mk, priors, valid_marker)

    ## "_blank_"    = species/fluor-matched secondary-only control (CTRL_OF)
    ## "_blankalt_" = the other secondary-only control in the same fluor
    cc_b  <- .blank_cc(CTRL_OF[[m]], cp)
    cc_ba <- .blank_cc(ALT_CTRL_OF[[m]], cp)

    if (cp == "nuc") {
      pre <- .feature_prefix(m, cp)
      .emit(cc_marker, nuc_ero, pre, "_")
      .emit(cc_b,      nuc_ero, pre, "_blank_")
      .emit(cc_ba,     nuc_ero, pre, "_blankalt_")
    } else {
      for (r_um in RING_UM_SET) {
        lab <- rings[[.ring_tag(r_um)]]
        pre <- .feature_prefix(m, cp, r_um)
        .emit(cc_marker, lab, pre, "_")
        .emit(cc_b,      lab, pre, "_blank_")
        .emit(cc_ba,     lab, pre, "_blankalt_")
      }
    }

    rm(mk, priors, cc_marker)
    gc(FALSE)
  }

  rm(plane_cache, blank_cc_cache, union_ring, comp_valid)

  corr_df <- do.call(cbind, corr_list)
  rm(corr_list)
  gc(FALSE)

  ## DAPI mean.
  dapi2 <- img_read(tiff, dapi_plane, x0, y0, w, h)
  dv <- .mean_by_label_valid(dapi2, nuc, nuc > 0, nmax)
  rm(dapi2)

  mo <- computeFeatures.moment(nuc)
  sh <- computeFeatures.shape(nuc)

  out <- data.frame(
    slide_id = sid,
    x_px = mo[, "m.cx"] + x0,
    y_px = mo[, "m.cy"] + y0,
    area = sh[, "s.area"],
    perim = sh[, "s.perimeter"],
    dapi_nuc = dv,
    corr_df,
    reg_cor = reg_qc$reg_cor,
    reg_cor_min = reg_qc$reg_cor_min,
    reg_cor_median = reg_qc$reg_cor_median,
    reg_cor_p05 = reg_qc$reg_cor_p05,
    reg_cor_n = reg_qc$reg_cor_n,
    reg_cor_n_low = reg_qc$reg_cor_n_low,
    reg_cor_low_frac = reg_qc$reg_cor_low_frac,
    reg_qc_informative = reg_qc$reg_qc_informative,
    reg_ref_sd = reg_qc$reg_ref_sd,
    reg_ref_range = reg_qc$reg_ref_range,
    tile_x = x0,
    tile_y = y0,
    stringsAsFactors = FALSE
  )

  rm(nuc, nuc_ero, rings, mo, sh, corr_df)
  gc(FALSE)

  dplyr::filter(out, !is.na(x_px))
}

###############################################################################
## RUN
###############################################################################

if (!exists("RUN_SLIDES")) RUN_SLIDES <- slides$slide_id

## Per-slide / merge modes for parallel runs (see header, change 5).
.env_slides <- Sys.getenv("PIPE_SLIDES", "")
if (nzchar(.env_slides)) RUN_SLIDES <- trimws(strsplit(.env_slides, ",")[[1]])
PIPE_MERGE <- identical(Sys.getenv("PIPE_MERGE", ""), "1")
PER_SLIDE_ONLY <- nzchar(.env_slides) && !PIPE_MERGE

stopifnot(all(RUN_SLIDES %in% slides$slide_id))

if (!PIPE_MERGE && !is.finite(MAX_TILES) && READ_MODE == "tifffile+pagecache") {
  stop(
    "Refusing full run in page-decode mode. Install zarr, restart R.",
    call. = FALSE
  )
}

## Pyramid pre-screen: a tile is processed only if its footprint on a
## down-sampled DAPI level contains acquired (non-zero) pixels above the
## slide-level tissue threshold. Avoids full-resolution reads of padding.
PRESCREEN_LEVEL <- 4L

.prescreen_tiles <- function(tiff, dapi_plane, g, tile) {
  lv <- tryCatch({
    st <- .tf$TiffFile(tiff)$series[0]
    sub <- st$levels[PRESCREEN_LEVEL]
    full_h <- .pynum(st$levels[0L]$shape)[2]
    arr <- reticulate::py_to_r(.np$asarray(
      sub$asarray(key = as.integer(dapi_plane - 1L))
    ))
    list(img = arr, f = full_h / nrow(arr))
  }, error = function(e) NULL)

  if (is.null(lv)) return(rep(TRUE, nrow(g)))

  img <- lv$img
  f <- lv$f
  nz <- img[img > 0]
  if (!length(nz)) return(rep(FALSE, nrow(g)))
  thr <- .otsu_vec(as.numeric(nz))

  vapply(seq_len(nrow(g)), function(i) {
    r0 <- floor(g$y[i] / f) + 1L
    r1 <- min(nrow(img), ceiling((g$y[i] + g$h[i]) / f))
    c0 <- floor(g$x[i] / f) + 1L
    c1 <- min(ncol(img), ceiling((g$x[i] + g$w[i]) / f))
    if (r0 > r1 || c0 > c1) return(FALSE)
    b <- img[r0:r1, c0:c1]
    sum(b > thr) >= 3L
  }, logical(1))
}

.segment_slide <- function(k) {
  sid <- slides$slide_id[k]
  cm  <- info[[sid]]$chmap

  ## Full grid including right/bottom partial tiles.
  xs <- seq(0, slides$sizeX[k] - 1L, by = TILE_SIZE)
  ys <- seq(0, slides$sizeY[k] - 1L, by = TILE_SIZE)
  g <- expand.grid(x = xs, y = ys)
  g$w <- pmin(TILE_SIZE, slides$sizeX[k] - g$x)
  g$h <- pmin(TILE_SIZE, slides$sizeY[k] - g$y)
  g <- g[g$w >= 256 & g$h >= 256, , drop = FALSE]
  n_total <- nrow(g)

  dapi_plane <- if (!is.null(cm$DAPI_AR)) cm$DAPI_AR else cm$DAPI_all[1]
  keep <- .prescreen_tiles(slides$tiff[k], dapi_plane, g, TILE_SIZE)
  g <- g[keep, , drop = FALSE]

  ## Optional sharding of one slide across processes: PIPE_SHARD="i/n"
  ## keeps tiles with (row index %% n) == i and writes a shard file.
  shard_tag <- ""
  shard <- Sys.getenv("PIPE_SHARD", "")
  if (nzchar(shard)) {
    sh <- as.integer(strsplit(shard, "/")[[1]])
    g <- g[(seq_len(nrow(g)) %% sh[2]) == sh[1], , drop = FALSE]
    shard_tag <- sprintf("_shard%dof%d", sh[1], sh[2])
  }

  if (is.finite(MAX_TILES)) {
    mid <- order(
      (g$x - slides$sizeX[k] / 2)^2 +
        (g$y - slides$sizeY[k] / 2)^2
    )
    g <- g[head(mid, MAX_TILES), , drop = FALSE]
  }

  message(
    "=== ", sid, " [", slides$group[k], "] : ",
    nrow(g), " tissue tiles of ", n_total, " grid tiles"
  )

  t0 <- Sys.time()

  out <- map_dfr(seq_len(nrow(g)), function(i) {
    ti <- Sys.time()

    if (i == 1 || i %% VERBOSE_EVERY == 0) {
      message("   starting tile ", i, "/", nrow(g), " at x=", g$x[i], " y=", g$y[i])
    }

    r <- tryCatch(
      process_tile(slides$tiff[k], cm, g$x[i], g$y[i], g$w[i], g$h[i], sid),
      error = function(e) {
        message("   tile ", i, " FAILED: ", conditionMessage(e))
        NULL
      }
    )

    dt <- as.numeric(difftime(Sys.time(), ti, units = "secs"))
    .tile_times <<- c(.tile_times, dt)

    if (i == 1 || i %% VERBOSE_EVERY == 0) {
      message(
        "   finished tile ", i, "/", nrow(g), " [", round(dt, 1), " s]",
        if (is.null(r)) " -> no cells" else paste0(" -> ", nrow(r), " cells")
      )
    }

    gc(FALSE)
    r
  })

  message(
    "   -> ", nrow(out), " cells in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "min")), 1), " min"
  )

  suffix <- paste0(shard_tag, if (is.finite(MAX_TILES)) "_SMOKETEST" else "")
  saveRDS(out, file.path(OUT, paste0("cells_raw_", sid, suffix, ".rds")))
  out
}

.tile_times <- c()

if (PIPE_MERGE) {
  message("\nMerge mode: reading per-slide cells_raw_<sid>.rds")
  cells_raw <- map_dfr(RUN_SLIDES, function(sid) {
    f <- file.path(OUT, paste0("cells_raw_", sid, ".rds"))
    shards <- list.files(OUT, pattern = paste0("^cells_raw_", sid, "_shard[0-9]+of[0-9]+\\.rds$"),
                         full.names = TRUE)
    if (length(shards)) {
      n_exp <- as.integer(sub(".*of([0-9]+)\\.rds$", "\\1", shards[1]))
      if (length(shards) != n_exp) stop(sid, ": found ", length(shards), " of ", n_exp, " shards")
      return(map_dfr(shards, readRDS))
    }
    if (!file.exists(f)) {
      warning("missing ", f, call. = FALSE)
      return(NULL)
    }
    readRDS(f)
  })
} else {
  message(
    "\nRunning slides: ", paste(RUN_SLIDES, collapse = ", "),
    " | max tiles = ", MAX_TILES
  )
  cells_raw <- map_dfr(which(slides$slide_id %in% RUN_SLIDES), .segment_slide)
}

if (PER_SLIDE_ONLY) {
  if (length(.tile_times)) {
    message("Median tile time: ", round(median(.tile_times), 1), " s; total ",
            round(sum(.tile_times) / 3600, 2), " h")
  }
  message("\nPer-slide mode: wrote cells_raw_<sid>.rds for ",
          paste(RUN_SLIDES, collapse = ", "), ". Run with PIPE_MERGE=1 to merge.")
  quit(save = "no", status = 0)
}

if (!nrow(cells_raw)) {
  stop("No cells segmented -- check tissue mask / thresholds.", call. = FALSE)
}

cells_raw <- left_join(
  cells_raw,
  slides[, c(
    "slide_id",
    "group",
    "run",
    "comet_sw",
    "has_negctrl",
    "tier2",
    "tier3"
  )],
  by = "slide_id"
)

out_file <- if (is.finite(MAX_TILES)) "cells_raw_SMOKETEST.rds" else "cells_raw.rds"
saveRDS(cells_raw, file.path(OUT, out_file))

###############################################################################
## QC OUTPUTS
###############################################################################

qc <- cells_raw |>
  group_by(slide_id, group) |>
  summarise(
    n_cells = n(),
    n_tiles = dplyr::n_distinct(paste(tile_x, tile_y)),
    n_reg_qc_informative_tiles = sum(reg_qc_informative %in% TRUE, na.rm = TRUE),
    n_informative_dapi_pairs = sum(reg_cor_n, na.rm = TRUE),
    median_dapi_cor = .median_or_na(reg_cor_median),
    p05_dapi_cor = .quantile_or_na(reg_cor_median, 0.05),
    min_dapi_cor = .min_or_na(reg_cor_min),
    frac_low_dapi_cor = suppressWarnings(sum(reg_cor_n_low, na.rm = TRUE) / pmax(sum(reg_cor_n, na.rm = TRUE), 1)),
    median_ref_sd = .median_or_na(reg_ref_sd),
    median_ref_range = .median_or_na(reg_ref_range),
    median_area_px = .median_or_na(area),
    median_area_um2 = .median_or_na(area) * PX_UM^2,
    .groups = "drop"
  )

cat("\nSegmentation QC:\n")
print(as.data.frame(qc))

write.csv(qc, file.path(OUT, "qc_registration.csv"), row.names = FALSE)

## Corrected primary-feature audit.
corr_audit <- do.call(rbind, lapply(unique(cells_raw$slide_id), function(sid) {
  d <- cells_raw[cells_raw$slide_id == sid, ]

  do.call(rbind, lapply(MARKERS_CORE, function(m) {
    cp <- MARKER_COMPARTMENT[[m]]
    pre <- .feature_prefix(m, cp)

    sub_col  <- paste0(pre, "_sub")
    frac_col <- paste0(pre, "_frac_pos")
    af_col   <- paste0(pre, "_af_frac")

    blank_sub_col  <- paste0(pre, "_blank_sub")
    blank_frac_col <- paste0(pre, "_blank_frac_pos")
    blank_af_col   <- paste0(pre, "_blank_af_frac")
    blank_intercept_col <- paste0(pre, "_blank_intercept")
    blank_slope_col     <- paste0(pre, "_blank_slope")
    blank_r2_col        <- paste0(pre, "_blank_r2")

    blank_ctrl_intercept_col <- paste0(pre, "_blank_blank_intercept")
    blank_ctrl_slope_col     <- paste0(pre, "_blank_blank_slope")
    blank_ctrl_r2_col        <- paste0(pre, "_blank_blank_r2")

    data.frame(
      slide_id = sid,
      marker = m,
      compartment = cp,
      radius_um = ifelse(cp == "ring", RING_UM_PRIMARY, NA),

      median_sub = round(.median_or_na(.col_or_na(d, sub_col)), 3),
      p95_sub = round(.qsafe(.col_or_na(d, sub_col), .95), 3),
      p99_sub = round(.qsafe(.col_or_na(d, sub_col), .99), 3),
      median_frac_pos = round(.median_or_na(.col_or_na(d, frac_col)), 4),
      p95_frac_pos = round(.qsafe(.col_or_na(d, frac_col), .95), 4),
      p99_frac_pos = round(.qsafe(.col_or_na(d, frac_col), .99), 4),
      median_af_frac = round(.median_or_na(.col_or_na(d, af_col)), 4),
      p99_af_frac = round(.qsafe(.col_or_na(d, af_col), .99), 4),
      median_blank_intercept = round(.median_or_na(.col_or_na(d, blank_intercept_col)), 4),
      median_blank_slope = round(.median_or_na(.col_or_na(d, blank_slope_col)), 4),
      median_blank_r2 = round(.median_or_na(.col_or_na(d, blank_r2_col)), 4),

      blank_median_sub = round(.median_or_na(.col_or_na(d, blank_sub_col)), 3),
      blank_p95_sub = round(.qsafe(.col_or_na(d, blank_sub_col), .95), 3),
      blank_p99_sub = round(.qsafe(.col_or_na(d, blank_sub_col), .99), 3),
      blank_median_frac_pos = round(.median_or_na(.col_or_na(d, blank_frac_col)), 4),
      blank_p95_frac_pos = round(.qsafe(.col_or_na(d, blank_frac_col), .95), 4),
      blank_p99_frac_pos = round(.qsafe(.col_or_na(d, blank_frac_col), .99), 4),
      blank_median_af_frac = round(.median_or_na(.col_or_na(d, blank_af_col)), 4),
      blank_p99_af_frac = round(.qsafe(.col_or_na(d, blank_af_col), .99), 4),
      blank_control_median_blank_intercept = round(.median_or_na(.col_or_na(d, blank_ctrl_intercept_col)), 4),
      blank_control_median_blank_slope = round(.median_or_na(.col_or_na(d, blank_ctrl_slope_col)), 4),
      blank_control_median_blank_r2 = round(.median_or_na(.col_or_na(d, blank_ctrl_r2_col)), 4),

      stringsAsFactors = FALSE
    )
  }))
}))

cat("\nCorrected primary marker-feature audit:\n")
print(corr_audit, row.names = FALSE)

write.csv(
  corr_audit,
  file.path(OUT, "qc_corrected_marker_features_primary_6um.csv"),
  row.names = FALSE
)

## Radius sensitivity feature audit before classification.
radius_feature_audit <- do.call(rbind, lapply(unique(cells_raw$slide_id), function(sid) {
  d <- cells_raw[cells_raw$slide_id == sid, ]

  ring_markers <- names(MARKER_COMPARTMENT)[MARKER_COMPARTMENT == "ring"]

  do.call(rbind, lapply(ring_markers, function(m) {
    do.call(rbind, lapply(RING_UM_SET, function(r_um) {
      pre <- .feature_prefix(m, "ring", r_um)

      sub_col  <- paste0(pre, "_sub")
      frac_col <- paste0(pre, "_frac_pos")
      af_col   <- paste0(pre, "_af_frac")

      blank_sub_col  <- paste0(pre, "_blank_sub")
      blank_frac_col <- paste0(pre, "_blank_frac_pos")
      blank_af_col   <- paste0(pre, "_blank_af_frac")
      blank_intercept_col <- paste0(pre, "_blank_intercept")
      blank_slope_col     <- paste0(pre, "_blank_slope")
      blank_r2_col        <- paste0(pre, "_blank_r2")

      blank_ctrl_intercept_col <- paste0(pre, "_blank_blank_intercept")
      blank_ctrl_slope_col     <- paste0(pre, "_blank_blank_slope")
      blank_ctrl_r2_col        <- paste0(pre, "_blank_blank_r2")

      data.frame(
        slide_id = sid,
        marker = m,
        radius_um = r_um,

        median_sub = round(.median_or_na(.col_or_na(d, sub_col)), 3),
        p95_sub = round(.qsafe(.col_or_na(d, sub_col), .95), 3),
        p99_sub = round(.qsafe(.col_or_na(d, sub_col), .99), 3),
        median_frac_pos = round(.median_or_na(.col_or_na(d, frac_col)), 4),
        p95_frac_pos = round(.qsafe(.col_or_na(d, frac_col), .95), 4),
        p99_frac_pos = round(.qsafe(.col_or_na(d, frac_col), .99), 4),
        median_af_frac = round(.median_or_na(.col_or_na(d, af_col)), 4),
        p99_af_frac = round(.qsafe(.col_or_na(d, af_col), .99), 4),
        median_blank_intercept = round(.median_or_na(.col_or_na(d, blank_intercept_col)), 4),
        median_blank_slope = round(.median_or_na(.col_or_na(d, blank_slope_col)), 4),
        median_blank_r2 = round(.median_or_na(.col_or_na(d, blank_r2_col)), 4),

        blank_median_sub = round(.median_or_na(.col_or_na(d, blank_sub_col)), 3),
        blank_p95_sub = round(.qsafe(.col_or_na(d, blank_sub_col), .95), 3),
        blank_p99_sub = round(.qsafe(.col_or_na(d, blank_sub_col), .99), 3),
        blank_median_frac_pos = round(.median_or_na(.col_or_na(d, blank_frac_col)), 4),
        blank_p95_frac_pos = round(.qsafe(.col_or_na(d, blank_frac_col), .95), 4),
        blank_p99_frac_pos = round(.qsafe(.col_or_na(d, blank_frac_col), .99), 4),
        blank_median_af_frac = round(.median_or_na(.col_or_na(d, blank_af_col)), 4),
        blank_p99_af_frac = round(.qsafe(.col_or_na(d, blank_af_col), .99), 4),
        blank_control_median_blank_intercept = round(.median_or_na(.col_or_na(d, blank_ctrl_intercept_col)), 4),
        blank_control_median_blank_slope = round(.median_or_na(.col_or_na(d, blank_ctrl_slope_col)), 4),
        blank_control_median_blank_r2 = round(.median_or_na(.col_or_na(d, blank_ctrl_r2_col)), 4),

        stringsAsFactors = FALSE
      )
    }))
  }))
}))

cat("\nRadius feature audit, ring markers:\n")
print(radius_feature_audit, row.names = FALSE)

write.csv(
  radius_feature_audit,
  file.path(OUT, "qc_radius_feature_summary_4_6_8um.csv"),
  row.names = FALSE
)

bad_reg <- qc$slide_id[
  (is.finite(qc$median_dapi_cor) & qc$median_dapi_cor < REG_LOW_COR_THRESHOLD) |
    (is.finite(qc$frac_low_dapi_cor) & qc$frac_low_dapi_cor > 0.25)
]
if (length(bad_reg)) {
  warning(
    "Low informative-window cross-cycle DAPI registration QC on: ",
    paste(bad_reg, collapse = ", "),
    call. = FALSE
  )
}

no_reg_info <- qc$slide_id[qc$n_reg_qc_informative_tiles == 0]
if (length(no_reg_info)) {
  warning(
    "No DAPI-informative registration QC windows on: ",
    paste(no_reg_info, collapse = ", "),
    " -- registration QC is unknown, not failed.",
    call. = FALSE
  )
}

if (any(qc$median_area_um2 < 15 | qc$median_area_um2 > 150, na.rm = TRUE)) {
  warning(
    "Median nuclear area outside 15-150 um^2 on: ",
    paste(
      qc$slide_id[
        qc$median_area_um2 < 15 |
          qc$median_area_um2 > 150
      ],
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (length(.tile_times)) {
  message("\nMedian tile time: ", round(median(.tile_times), 1), " s")
  message("Total tiles: ", length(.tile_times))
  message("Total tile time: ", round(sum(.tile_times) / 3600, 2), " h")
}

message("\n02_segment_measure.R complete -> ", file.path(OUT, out_file))
