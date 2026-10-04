###############################################################################
# 03_classify.R  v12
#
# Classification using corrected features from 02_segment_measure.R v12.
#
# v12 CHANGES
#   - Restores the file. Commit 7808e70 ("Refactor marker usability and
#     population call logic") replaced 03 with its last ~210 lines only: the
#     input/QC-filter section, gate_marker_one_slide(), qc_gates and the
#     output section (cells_classified.rds) were deleted, so 03 could not run
#     and 04 was reading a stale cells_classified.rds from v11.4.1.
#     Structure below = v11.4.1 header + v11.4.2 population logic (marker
#     usability audit-only) + v11.4.1 outputs.
#   - Marker threshold anchored to the matched blank pseudo-marker
#     (BLANK_FP_TARGET), and blank false-positive rate measured at the
#     marker's own threshold. See gate_marker_one_slide().
#   - Marker usability remains AUDIT ONLY (cell calls are not suppressed),
#     but each slide x marker now carries a confidence flag
#     (enrichment of marker over blank at the same threshold) that 04/05
#     print on every figure.
#   - ASPC_stromal (PDGFa+ panCK-) added as a sensitivity definition; the
#     primary ASPC call is unchanged (PDGFa+).
###############################################################################


if (!exists("OUT")) OUT <- "results_ASPC_NEPC"

suppressPackageStartupMessages({
  library(dplyr)
})

###############################################################################
## INPUT
###############################################################################

infile <- file.path(OUT, "cells_raw.rds")

if (!file.exists(infile)) {
  stop(
    "Missing cells_raw.rds. Re-run revised 02_segment_measure.R first.",
    call. = FALSE
  )
}

cells <- readRDS(infile)

if (!is.data.frame(cells) || !nrow(cells)) {
  stop("cells_raw.rds has zero rows or is not a data.frame.", call. = FALSE)
}

EXPECTED_SLIDES <- c(
  "PM0T32", "PM677_A8", "PM12_B3",
  "PM1024_B", "PM423", "PM159T5"
)

.present <- unique(cells$slide_id)
.absent  <- setdiff(EXPECTED_SLIDES, .present)
.extra   <- setdiff(.present, EXPECTED_SLIDES)

if (length(.absent)) {
  warning(
    "cells_raw.rds is missing expected slide(s): ",
    paste(.absent, collapse = ", "),
    " -- was 02_segment_measure.R run on the full set?",
    call. = FALSE
  )
}

if (length(.extra)) {
  warning(
    "cells_raw.rds contains unexpected slide(s): ",
    paste(.extra, collapse = ", "),
    call. = FALSE
  )
}

if ("has_negctrl" %in% names(cells)) {
  .no_nc <- unique(cells$slide_id[!(cells$has_negctrl %in% TRUE)])
  .no_nc <- .no_nc[!is.na(.no_nc)]

  if (length(.no_nc)) {
    stop(
      "Slide(s) without negative controls in cells_raw.rds: ",
      paste(.no_nc, collapse = ", "),
      "\nThis indicates a stale segmentation. Re-run 01_channels.R and ",
      "02_segment_measure.R against the current acquisitions.",
      call. = FALSE
    )
  }
} else {
  stop(
    "cells_raw.rds is missing has_negctrl. Re-run revised 02_segment_measure.R first.",
    call. = FALSE
  )
}

PX_UM <- 12500 / 44643
if ("area_um2" %in% names(cells)) {
  PX_UM <- sqrt(median(cells$area_um2 / cells$area, na.rm = TRUE))
}

###############################################################################
## QC FILTER
###############################################################################

MIN_AREA_UM2 <- 8
MAX_AREA_UM2 <- 400
MIN_REG_COR  <- 0.85

cells$area_um2 <- cells$area * PX_UM^2

if (!"reg_qc_informative" %in% names(cells)) {
  cells$reg_qc_informative <- is.finite(cells$reg_cor)
}

keep <- is.finite(cells$area_um2) &
        cells$area_um2 >= MIN_AREA_UM2 &
        cells$area_um2 <= MAX_AREA_UM2 &
        (
          !(cells$reg_qc_informative %in% TRUE) |
            is.na(cells$reg_cor) |
            cells$reg_cor >= MIN_REG_COR
        )

cat("QC keep:", sum(keep), "/", nrow(cells), "cells\n")

cells <- cells[keep, , drop = FALSE]

if (!nrow(cells)) {
  stop("No cells remain after QC filtering.", call. = FALSE)
}

###############################################################################
## MARKER ROUTING
###############################################################################

RING_UM_PRIMARY <- 6
RING_UM_SET <- c(4, 6, 8)

MARKER_COMPARTMENT <- c(
  CEBPg = "nuc",
  AR    = "nuc",
  PDGF  = "ring",
  panCK = "ring",
  FOLH1 = "ring",
  DLL3  = "ring"
)

MIN_FRAC <- c(
  CEBPg = 0.08,
  AR    = 0.08,
  PDGF  = 0.06,
  panCK = 0.06,
  FOLH1 = 0.06,
  DLL3  = 0.05
)

HIGH_FRAC <- c(
  CEBPg = NA_real_,
  AR    = NA_real_,
  PDGF  = NA_real_,
  panCK = NA_real_,
  FOLH1 = NA_real_,
  DLL3  = NA_real_
)

MAX_AF_FRAC <- c(
  CEBPg = 0.05,
  AR    = 0.05,
  PDGF  = 0.10,
  panCK = 0.10,
  FOLH1 = 0.10,
  DLL3  = 0.10
)

MIN_MARKER_TO_BLANK_ENRICHMENT <- 2.0
MAX_BLANK_FRAC_POS <- 0.10
MAX_GAIN_R2_FOR_WEAK_ENRICHMENT <- 0.90
GAIN_SLOPE_TOL <- 0.05

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

.required_cols_for <- function(marker, r_um = RING_UM_PRIMARY, blank = FALSE) {
  cp <- MARKER_COMPARTMENT[[marker]]
  pre <- .feature_prefix(marker, cp, r_um)

  if (blank) {
    c(
      paste0(pre, "_blank_sub"),
      paste0(pre, "_blank_frac_pos"),
      paste0(pre, "_blank_af_frac")
    )
  } else {
    c(
      paste0(pre, "_sub"),
      paste0(pre, "_frac_pos"),
      paste0(pre, "_af_frac")
    )
  }
}

all_required <- unique(unlist(lapply(names(MARKER_COMPARTMENT), function(m) {
  if (MARKER_COMPARTMENT[[m]] == "nuc") {
    c(
      .required_cols_for(m, blank = FALSE),
      .required_cols_for(m, blank = TRUE)
    )
  } else {
    unlist(lapply(RING_UM_SET, function(r) {
      c(
        .required_cols_for(m, r, blank = FALSE),
        .required_cols_for(m, r, blank = TRUE)
      )
    }))
  }
})))

missing <- setdiff(all_required, names(cells))

if (length(missing)) {
  stop(
    "Missing v11 corrected feature columns from 02:\n  ",
    paste(missing, collapse = "\n  "),
    "\nRe-run 02_segment_measure.R v11.3+ before 03.",
    call. = FALSE
  )
}

###############################################################################
## HELPERS
###############################################################################

.qsafe <- function(x, p, default = NA_real_) {
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
}

.median_col <- function(d, col, default = NA_real_) {
  if (!col %in% names(d)) return(default)
  x <- d[[col]]
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  median(x, na.rm = TRUE)
}

.robust_thr <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 50) return(Inf)

  med <- median(x, na.rm = TRUE)
  sc  <- mad(x, constant = 1.4826, na.rm = TRUE)

  if (!is.finite(sc) || sc <= 0) {
    sc <- sd(x, na.rm = TRUE)
  }

  if (!is.finite(sc) || sc <= 0) {
    return(.qsafe(x, 0.95, default = Inf))
  }

  med + 3 * sc
}

.get_feature_cols <- function(marker, r_um = RING_UM_PRIMARY, blank = FALSE) {
  cp <- MARKER_COMPARTMENT[[marker]]
  pre <- .feature_prefix(marker, cp, r_um)

  if (blank) {
    list(
      score = paste0(pre, "_blank_sub"),
      frac  = paste0(pre, "_blank_frac_pos"),
      af    = paste0(pre, "_blank_af_frac")
    )
  } else {
    list(
      score = paste0(pre, "_sub"),
      frac  = paste0(pre, "_frac_pos"),
      af    = paste0(pre, "_af_frac")
    )
  }
}


.bool <- function(x) {
  y <- x %in% TRUE
  y[is.na(y)] <- FALSE
  y
}


## v12 gate. A cell is marker-positive when its corrected score exceeds
##   thr = max( robust marker threshold (median + 3 MAD),
##              upper (1 - BLANK_FP_TARGET) quantile of the matched blank
##              pseudo-marker score in the same compartment )
## and its positive-pixel fraction >= MIN_FRAC.
## v11 used the robust marker threshold alone. Because that is relative to
## the slide's own distribution it always calls ~5-15% of cells positive,
## and the matched blank processed the same way was called positive just as
## often (blank/marker enrichment < 1 for 20 of 36 slide x marker pairs).
## Anchoring to the blank caps the blank-channel false-positive rate at
## BLANK_FP_TARGET by construction, and the blank is then evaluated at the
## SAME threshold as the marker, which makes the enrichment ratio honest.
BLANK_FP_TARGET <- 0.01

.gate_eligible <- function(d, marker, r_um, blank) {
  cols <- .get_feature_cols(marker, r_um = r_um, blank = blank)
  score <- d[[cols$score]]
  frac  <- d[[cols$frac]]
  af    <- d[[cols$af]]
  list(
    score = score, frac = frac, af = af,
    eligible = is.finite(score) & is.finite(frac) & is.finite(af) &
               af <= MAX_AF_FRAC[[marker]]
  )
}

.gate_summary <- function(e, pos, thr, marker) {
  list(
    pos = pos,
    thr = thr,
    n_eligible = sum(e$eligible),
    n_pos = sum(pos, na.rm = TRUE),
    frac_pos = mean(pos, na.rm = TRUE),
    median_score = median(e$score, na.rm = TRUE),
    p95_score = .qsafe(e$score, 0.95),
    p99_score = .qsafe(e$score, 0.99),
    median_frac = median(e$frac, na.rm = TRUE),
    p95_frac = .qsafe(e$frac, 0.95),
    p99_frac = .qsafe(e$frac, 0.99),
    median_af_frac = median(e$af, na.rm = TRUE),
    p99_af_frac = .qsafe(e$af, 0.99)
  )
}

## Returns list(marker = <gate summary>, blank = <gate summary at the SAME thr>,
##              thr_robust, thr_blank)
gate_marker_one_slide <- function(d, marker, r_um = RING_UM_PRIMARY, blank = FALSE) {
  em <- .gate_eligible(d, marker, r_um, blank = FALSE)
  eb <- .gate_eligible(d, marker, r_um, blank = TRUE)

  if (sum(em$eligible) < 50) {
    none <- rep(FALSE, nrow(d))
    g <- list(marker = .gate_summary(em, none, Inf, marker),
              blank = .gate_summary(eb, none, Inf, marker),
              thr_robust = Inf, thr_blank = Inf)
    return(if (blank) g$blank else g$marker)
  }

  thr_robust <- .robust_thr(em$score[em$eligible])
  if (!is.finite(thr_robust) || thr_robust <= 0) {
    thr_robust <- .qsafe(em$score[em$eligible], 0.95, default = Inf)
  }

  thr_blank <- .qsafe(eb$score[eb$eligible], 1 - BLANK_FP_TARGET, default = -Inf)

  thr <- max(thr_robust, thr_blank, na.rm = TRUE)

  pos_m <- em$eligible & em$score > thr & em$frac >= MIN_FRAC[[marker]]
  pos_b <- eb$eligible & eb$score > thr & eb$frac >= MIN_FRAC[[marker]]
  pos_m[is.na(pos_m)] <- FALSE
  pos_b[is.na(pos_b)] <- FALSE

  out <- list(
    marker = .gate_summary(em, pos_m, thr, marker),
    blank  = .gate_summary(eb, pos_b, thr, marker),
    thr_robust = thr_robust,
    thr_blank = thr_blank
  )
  out$marker$thr_robust <- thr_robust
  out$marker$thr_blank  <- thr_blank
  out$blank$thr_robust  <- thr_robust
  out$blank$thr_blank   <- thr_blank

  if (blank) out$blank else out$marker
}

###############################################################################
## PRIMARY MARKER GATES, 6 um
###############################################################################

for (m in names(MARKER_COMPARTMENT)) {
  cells[[paste0(m, "_pos")]] <- FALSE
  cells[[paste0(m, "_blank_pos")]] <- FALSE
}

gates <- list()
blank_gates <- list()

.gate_row <- function(sid, m, g, control_type) {
  data.frame(
    slide_id = sid,
    marker = m,
    radius_um = ifelse(MARKER_COMPARTMENT[[m]] == "nuc", NA, RING_UM_PRIMARY),
    compartment = MARKER_COMPARTMENT[[m]],
    control_type = control_type,
    threshold_sub = g$thr,
    threshold_robust = if (is.null(g$thr_robust)) NA_real_ else g$thr_robust,
    threshold_blank_q = if (is.null(g$thr_blank)) NA_real_ else g$thr_blank,
    n_eligible = g$n_eligible,
    n_pos = g$n_pos,
    frac_pos = g$frac_pos,
    median_sub = g$median_score,
    p95_sub = g$p95_score,
    p99_sub = g$p99_score,
    median_frac_pixels_pos = g$median_frac,
    p95_frac_pixels_pos = g$p95_frac,
    p99_frac_pixels_pos = g$p99_frac,
    median_af_frac = g$median_af_frac,
    p99_af_frac = g$p99_af_frac,
    min_frac = MIN_FRAC[[m]],
    max_af_frac = MAX_AF_FRAC[[m]],
    stringsAsFactors = FALSE
  )
}

for (sid in unique(cells$slide_id)) {
  ii <- which(cells$slide_id == sid)
  d  <- cells[ii, , drop = FALSE]

  for (m in names(MARKER_COMPARTMENT)) {
    g  <- gate_marker_one_slide(d, marker = m, r_um = RING_UM_PRIMARY, blank = FALSE)
    bg <- gate_marker_one_slide(d, marker = m, r_um = RING_UM_PRIMARY, blank = TRUE)

    cells[ii, paste0(m, "_pos")] <- g$pos
    cells[ii, paste0(m, "_blank_pos")] <- bg$pos

    gates[[length(gates) + 1]] <- .gate_row(sid, m, g, "marker")
    blank_gates[[length(blank_gates) + 1]] <- .gate_row(sid, m, bg, "blank_as_marker")
  }
}

qc_gates <- bind_rows(gates)
qc_blank_false_positive <- bind_rows(blank_gates)

write.csv(qc_gates, file.path(OUT, "qc_gates_primary_6um.csv"), row.names = FALSE)
write.csv(qc_blank_false_positive, file.path(OUT, "qc_blank_false_positive_primary_6um.csv"),
          row.names = FALSE)

cat("\n=== Primary marker gates, 6 um:\n")
print(qc_gates[, c("slide_id", "marker", "threshold_sub", "threshold_robust",
                   "threshold_blank_q", "n_eligible", "n_pos", "frac_pos")])

cat("\n=== Blank pseudo-marker at the marker threshold, 6 um:\n")
print(qc_blank_false_positive[, c("slide_id", "marker", "threshold_sub",
                                  "n_eligible", "n_pos", "frac_pos")])


###############################################################################
## MARKER USABILITY QC, 6 um -- AUDIT ONLY, NOT USED FOR CELL CALLING
###############################################################################

marker_usability <- do.call(rbind, lapply(unique(cells$slide_id), function(sid) {
  d <- cells[cells$slide_id == sid, , drop = FALSE]

  do.call(rbind, lapply(names(MARKER_COMPARTMENT), function(m) {
    cp <- MARKER_COMPARTMENT[[m]]
    pre <- .feature_prefix(m, cp, RING_UM_PRIMARY)

    marker_row <- qc_gates[
      qc_gates$slide_id == sid &
        qc_gates$marker == m &
        qc_gates$control_type == "marker",
      , drop = FALSE
    ]

    blank_row <- qc_blank_false_positive[
      qc_blank_false_positive$slide_id == sid &
        qc_blank_false_positive$marker == m &
        qc_blank_false_positive$control_type == "blank_as_marker",
      , drop = FALSE
    ]

    marker_frac <- if (nrow(marker_row)) marker_row$frac_pos[1] else NA_real_
    blank_frac  <- if (nrow(blank_row))  blank_row$frac_pos[1]  else NA_real_

    ## Same threshold for both, so this is a like-for-like enrichment.
    ## Floor the blank at 1/n_eligible to avoid Inf for a zero blank count.
    n_el <- if (nrow(blank_row)) max(blank_row$n_eligible[1], 1) else 1
    enrichment <- marker_frac / max(blank_frac, 1 / n_el)

    slope <- .median_col(d, paste0(pre, "_blank_slope"))
    r2    <- .median_col(d, paste0(pre, "_blank_r2"))

    data.frame(
      slide_id = sid,
      marker = m,
      marker_frac_pos = marker_frac,
      blank_frac_pos = blank_frac,
      enrichment_over_blank = enrichment,
      median_blank_slope = slope,
      median_blank_r2 = r2,

      ## Important:
      ## marker_usable is retained only for QC reporting.
      ## It is NOT used to suppress cell calls.
      marker_usable = TRUE,
      confidence = dplyr::case_when(
        !is.finite(enrichment)      ~ "unknown",
        enrichment >= MIN_MARKER_TO_BLANK_ENRICHMENT * 2.5 ~ "high",
        enrichment >= MIN_MARKER_TO_BLANK_ENRICHMENT ~ "moderate",
        TRUE                        ~ "low"
      ),
      fail_reason = "audit_only_not_used_for_calling",

      stringsAsFactors = FALSE
    )
  }))
}))

cat("\n=== Marker usability audit, primary 6 um; NOT used for cell calling:\n")
print(marker_usability, row.names = FALSE)

write.csv(
  marker_usability,
  file.path(OUT, "qc_marker_usability_primary_6um.csv"),
  row.names = FALSE
)

## Add usability columns for compatibility, but keep all markers callable.
for (m in names(MARKER_COMPARTMENT)) {
  usable_col <- paste0(m, "_usable")
  cells[[usable_col]] <- TRUE
}

###############################################################################
## RADIUS SENSITIVITY, 4 / 6 / 8 um
###############################################################################

sens <- list()

ring_markers <- names(MARKER_COMPARTMENT)[MARKER_COMPARTMENT == "ring"]

for (sid in unique(cells$slide_id)) {
  ii <- which(cells$slide_id == sid)
  d  <- cells[ii, , drop = FALSE]

  for (r_um in RING_UM_SET) {
    for (m in ring_markers) {
      g <- gate_marker_one_slide(
        d,
        marker = m,
        r_um = r_um,
        blank = FALSE
      )

      bg <- gate_marker_one_slide(
        d,
        marker = m,
        r_um = r_um,
        blank = TRUE
      )

      sens[[length(sens) + 1]] <- data.frame(
        slide_id = sid,
        marker = m,
        radius_um = r_um,
        marker_frac_pos = g$frac_pos,
        marker_n_pos = g$n_pos,
        marker_n_eligible = g$n_eligible,
        blank_frac_pos = bg$frac_pos,
        ## bg is evaluated at the marker threshold (g$thr) in v12
        blank_n_pos = bg$n_pos,
        blank_n_eligible = bg$n_eligible,
        enrichment_over_blank = ifelse(bg$frac_pos > 0, g$frac_pos / bg$frac_pos, Inf),
        marker_threshold_sub = g$thr,
        blank_threshold_sub = bg$thr,
        marker_p99_sub = g$p99_score,
        blank_p99_sub = bg$p99_score,
        stringsAsFactors = FALSE
      )
    }
  }
}

qc_radius_sensitivity <- bind_rows(sens)

write.csv(
  qc_radius_sensitivity,
  file.path(OUT, "qc_radius_sensitivity_4_6_8um.csv"),
  row.names = FALSE
)

cat("\n=== Radius sensitivity, ring markers:\n")
print(qc_radius_sensitivity)

###############################################################################
## POPULATION CALLS -- PIXEL-CORRECTED CELL CALLS, NO SLIDE-LEVEL BLANK QC GATE
###############################################################################

## Raw marker calls are already based on:
##   - pixel-level matched-blank subtraction from 02
##   - AF puncta masking from 02
##   - corrected per-cell sub and frac_pos gates in 03
##
## Therefore do not suppress cells based on slide-level mean blank positivity.

for (m in names(MARKER_COMPARTMENT)) {
  pos_col <- paste0(m, "_pos")
  qc_col <- paste0(m, "_pos_qc")

  ## Compatibility column.
  ## Since slide-level marker usability QC is audit-only now,
  ## the QC-gated call is identical to the corrected marker call.
  cells[[qc_col]] <- cells[[pos_col]] %in% TRUE
}

cells$CEBPg_neg <- !(cells$CEBPg_pos %in% TRUE)
cells$CEBPg_neg_qc <- !(cells$CEBPg_pos_qc %in% TRUE)

## ASPC:
## Requested definition simplifies to PDGF_pos.
## Final ASPC now uses corrected PDGF cell calls directly.
cells$ASPC_raw <- cells$PDGF_pos %in% TRUE
cells$ASPC_callable <- TRUE
cells$ASPC_not_callable <- FALSE
cells$ASPC <- cells$ASPC_raw

cells$ASPC_CEBPg_pos <- cells$ASPC & (cells$CEBPg_pos %in% TRUE)
cells$ASPC_CEBPg_neg <- cells$ASPC & !(cells$CEBPg_pos %in% TRUE)

## Sensitivity definition: PDGFa ligand is also made by tumour epithelium,
## so a stromal-restricted ASPC call excludes panCK+ cells.
cells$ASPC_stromal <- cells$ASPC & !(cells$panCK_pos %in% TRUE)

cells$ASPC_raw_CEBPg_pos <- cells$ASPC_raw & (cells$CEBPg_pos %in% TRUE)
cells$ASPC_raw_CEBPg_neg <- cells$ASPC_raw & !(cells$CEBPg_pos %in% TRUE)

## Compatibility columns for scripts that expect them.
cells$ASPC_CEBPg_known <- cells$ASPC
cells$ASPC_CEBPg_unknown <- FALSE
cells$ASPC_raw_CEBPg_known <- cells$ASPC_raw
cells$ASPC_raw_CEBPg_unknown <- FALSE

## Epithelial cancer.
cells$EPI_CA_raw <- (cells$panCK_pos %in% TRUE) &
                    ((cells$FOLH1_pos %in% TRUE) | (cells$AR_pos %in% TRUE))

cells$EPI_CA_callable <- TRUE
cells$EPI_CA_not_callable <- FALSE
cells$EPI_CA <- cells$EPI_CA_raw

## Neuroendocrine cancer.
cells$NE_CA_raw <- cells$DLL3_pos %in% TRUE
cells$NE_CA_callable <- TRUE
cells$NE_CA_not_callable <- FALSE
cells$NE_CA <- cells$NE_CA_raw

## Blank-derived false-positive population controls.
## These are retained as QC/audit metrics only. They do not suppress marker calls.
cells$ASPC_blank <- cells$PDGF_blank_pos %in% TRUE

cells$EPI_CA_blank <- (cells$panCK_blank_pos %in% TRUE) &
                      ((cells$FOLH1_blank_pos %in% TRUE) |
                       (cells$AR_blank_pos %in% TRUE))

cells$NE_CA_blank <- cells$DLL3_blank_pos %in% TRUE

## Display label.
cells$pop <- "Other"

cells$pop[cells$ASPC]   <- "ASPC"
cells$pop[cells$EPI_CA] <- "EPI_CA"
cells$pop[cells$NE_CA]  <- "NE_CA"

cells$n_pop_flags <- rowSums(cbind(
  ASPC = cells$ASPC,
  EPI_CA = cells$EPI_CA,
  NE_CA = cells$NE_CA
), na.rm = TRUE)

cells$multi_pop <- cells$n_pop_flags > 1
cells$pop[cells$multi_pop] <- "Multi"


###############################################################################
## OUTPUT HELPERS
###############################################################################

slide_levels <- names(table(cells$slide_id))
n_by_slide <- as.integer(table(cells$slide_id))

.sum_by_slide <- function(x) {
  if (is.character(x)) {
    if (!x %in% names(cells)) {
      return(rep(NA_integer_, length(slide_levels)))
    }
    x <- cells[[x]]
  }

  if (length(x) != nrow(cells)) {
    return(rep(NA_integer_, length(slide_levels)))
  }

  as.integer(tapply(.bool(x), cells$slide_id, sum, na.rm = TRUE)[slide_levels])
}

.mean_by_slide <- function(x) {
  if (is.character(x)) {
    if (!x %in% names(cells)) {
      return(rep(NA_real_, length(slide_levels)))
    }
    x <- cells[[x]]
  }

  if (length(x) != nrow(cells)) {
    return(rep(NA_real_, length(slide_levels)))
  }

  as.numeric(tapply(.bool(x), cells$slide_id, mean, na.rm = TRUE)[slide_levels])
}

.all_by_slide <- function(x) {
  if (is.character(x)) {
    if (!x %in% names(cells)) {
      return(rep(FALSE, length(slide_levels)))
    }
    x <- cells[[x]]
  }

  if (length(x) != nrow(cells)) {
    return(rep(FALSE, length(slide_levels)))
  }

  as.logical(tapply(.bool(x), cells$slide_id, function(v) all(v %in% TRUE))[slide_levels])
}

###############################################################################
## OUTPUTS
###############################################################################

saveRDS(cells, file.path(OUT, "cells_classified.rds"))

cat("\n=== Cells per slide after QC:\n")
print(table(cells$slide_id))

cat("\n=== Single-marker positive fraction by slide, primary 6 um:\n")

marker_frac <- sapply(names(MARKER_COMPARTMENT), function(m) {
  tapply(cells[[paste0(m, "_pos")]], cells$slide_id, mean, na.rm = TRUE)
})

print(round(marker_frac, 4))

write.csv(
  data.frame(slide_id = rownames(marker_frac), marker_frac, row.names = NULL),
  file.path(OUT, "marker_positive_fraction_by_slide_primary_6um.csv"),
  row.names = FALSE
)

cat("\n=== Single-marker QC-gated positive fraction by slide, primary 6 um:\n")

marker_frac_qc <- sapply(names(MARKER_COMPARTMENT), function(m) {
  tapply(cells[[paste0(m, "_pos_qc")]], cells$slide_id, mean, na.rm = TRUE)
})

print(round(marker_frac_qc, 4))

write.csv(
  data.frame(slide_id = rownames(marker_frac_qc), marker_frac_qc, row.names = NULL),
  file.path(OUT, "marker_qc_positive_fraction_by_slide_primary_6um.csv"),
  row.names = FALSE
)

cat("\n=== Single-marker blank false-positive fraction by slide, primary 6 um:\n")

blank_frac <- sapply(names(MARKER_COMPARTMENT), function(m) {
  tapply(cells[[paste0(m, "_blank_pos")]], cells$slide_id, mean, na.rm = TRUE)
})

print(round(blank_frac, 4))

write.csv(
  data.frame(slide_id = rownames(blank_frac), blank_frac, row.names = NULL),
  file.path(OUT, "marker_blank_false_positive_fraction_by_slide_primary_6um.csv"),
  row.names = FALSE
)

population_counts <- data.frame(
  slide_id = slide_levels,
  n_cells = n_by_slide,

  ASPC = .sum_by_slide("ASPC"),
  ASPC_stromal = .sum_by_slide("ASPC_stromal"),
  ASPC_raw = .sum_by_slide("ASPC_raw"),
  ASPC_not_callable = .sum_by_slide("ASPC_not_callable"),
  ASPC_callable = .all_by_slide("ASPC_callable"),

  ASPC_CEBPg_known = .sum_by_slide("ASPC_CEBPg_known"),
  ASPC_CEBPg_pos = .sum_by_slide("ASPC_CEBPg_pos"),
  ASPC_CEBPg_neg = .sum_by_slide("ASPC_CEBPg_neg"),
  ASPC_CEBPg_unknown = .sum_by_slide("ASPC_CEBPg_unknown"),

  ASPC_raw_CEBPg_known = .sum_by_slide("ASPC_raw_CEBPg_known"),
  ASPC_raw_CEBPg_pos = .sum_by_slide("ASPC_raw_CEBPg_pos"),
  ASPC_raw_CEBPg_neg = .sum_by_slide("ASPC_raw_CEBPg_neg"),
  ASPC_raw_CEBPg_unknown = .sum_by_slide("ASPC_raw_CEBPg_unknown"),

  EPI_CA = .sum_by_slide("EPI_CA"),
  EPI_CA_raw = .sum_by_slide("EPI_CA_raw"),
  EPI_CA_not_callable = .sum_by_slide("EPI_CA_not_callable"),
  EPI_CA_callable = .all_by_slide("EPI_CA_callable"),

  NE_CA = .sum_by_slide("NE_CA"),
  NE_CA_raw = .sum_by_slide("NE_CA_raw"),
  NE_CA_not_callable = .sum_by_slide("NE_CA_not_callable"),
  NE_CA_callable = .all_by_slide("NE_CA_callable"),

  Other = .sum_by_slide(cells$pop == "Other"),
  Multi = .sum_by_slide("multi_pop"),

  ASPC_blank = .sum_by_slide("ASPC_blank"),
  EPI_CA_blank = .sum_by_slide("EPI_CA_blank"),
  NE_CA_blank = .sum_by_slide("NE_CA_blank"),

  row.names = NULL
)

cat("\n=== Population counts by slide, primary 6 um:\n")
print(population_counts)

write.csv(
  population_counts,
  file.path(OUT, "population_counts_by_slide_primary_6um.csv"),
  row.names = FALSE
)

population_fractions <- population_counts

for (cl in setdiff(
  names(population_fractions),
  c("slide_id", "n_cells", "ASPC_callable", "EPI_CA_callable", "NE_CA_callable")
)) {
  population_fractions[[cl]] <- population_fractions[[cl]] / population_fractions$n_cells
}

cat("\n=== Population fractions by slide, primary 6 um:\n")
print(population_fractions)

write.csv(
  population_fractions,
  file.path(OUT, "population_fractions_by_slide_primary_6um.csv"),
  row.names = FALSE
)

aspc_callable_summary <- cells |>
  group_by(slide_id) |>
  summarise(
    ASPC_callable = all(.data$ASPC_callable %in% TRUE),
    PDGF_usable = all(.data$PDGF_usable %in% TRUE),
    ASPC = sum(.data$ASPC, na.rm = TRUE),
    ASPC_raw = sum(.data$ASPC_raw, na.rm = TRUE),
    ASPC_not_callable = sum(.data$ASPC_not_callable, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== ASPC callability by slide:\n")
print(as.data.frame(aspc_callable_summary), row.names = FALSE)

write.csv(
  aspc_callable_summary,
  file.path(OUT, "aspc_callability_by_slide.csv"),
  row.names = FALSE
)

population_callability_summary <- cells |>
  group_by(slide_id) |>
  summarise(
    ASPC_callable = all(.data$ASPC_callable %in% TRUE),
    EPI_CA_callable = all(.data$EPI_CA_callable %in% TRUE),
    NE_CA_callable = all(.data$NE_CA_callable %in% TRUE),

    PDGF_usable = all(.data$PDGF_usable %in% TRUE),
    panCK_usable = all(.data$panCK_usable %in% TRUE),
    FOLH1_usable = all(.data$FOLH1_usable %in% TRUE),
    AR_usable = all(.data$AR_usable %in% TRUE),
    DLL3_usable = all(.data$DLL3_usable %in% TRUE),
    CEBPg_usable = all(.data$CEBPg_usable %in% TRUE),

    n_ASPC = sum(.data$ASPC, na.rm = TRUE),
    n_ASPC_raw = sum(.data$ASPC_raw, na.rm = TRUE),
    n_ASPC_not_callable = sum(.data$ASPC_not_callable, na.rm = TRUE),

    n_EPI_CA = sum(.data$EPI_CA, na.rm = TRUE),
    n_EPI_CA_raw = sum(.data$EPI_CA_raw, na.rm = TRUE),
    n_EPI_CA_not_callable = sum(.data$EPI_CA_not_callable, na.rm = TRUE),

    n_NE_CA = sum(.data$NE_CA, na.rm = TRUE),
    n_NE_CA_raw = sum(.data$NE_CA_raw, na.rm = TRUE),
    n_NE_CA_not_callable = sum(.data$NE_CA_not_callable, na.rm = TRUE),

    .groups = "drop"
  )

cat("\n=== Population callability by slide:\n")
print(as.data.frame(population_callability_summary), row.names = FALSE)

write.csv(
  population_callability_summary,
  file.path(OUT, "population_callability_by_slide.csv"),
  row.names = FALSE
)

cat("\n=== Overlapping final callable population flags:\n")
print(table(cells$n_pop_flags))

cat("\n=== Display population labels:\n")
print(table(cells$pop))

message("\n03_classify.R complete -> ", file.path(OUT, "cells_classified.rds"))
