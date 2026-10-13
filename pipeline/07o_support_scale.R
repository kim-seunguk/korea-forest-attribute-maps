# 07o_support_scale.R  [change of support: over how large an area does each map become usable?]
# ============================================================================================
# A 30 m map validated at plot level answers one question; a user summarising a watershed asks
# another. This turns that difference into a curve. The spatial-CV out-of-fold predictions are
# aggregated into progressively coarser square units and R2, RMSE and MBE are recomputed at each.
#
# Aggregation raises R2 on its own, so two controls are carried:
#   (1) RMSE and MBE in the original units alongside it, which averaging cannot inflate the way it
#       inflates a variance ratio.
#   (2) a permutation control — predictions shuffled across the country, then the same aggregation,
#       which should leave R2 at or below zero at every support.
#
# ! Aggregating to units that hold one or two plots lets a handful of sparse units carry the whole
#   statistic. Results are written for each value in MIN_PLOTS_SET, with min_plots as an output column.
#
# In  : outputs/<VERSION>/m4a_pred_points<IN_SUF>.csv  (per-plot spatial-CV out-of-fold obs/pred, direct route)
#       outputs/<VERSION>/cv_folds.csv                 (Plot_ID, UTMK_X, UTMK_Y in EPSG:5179 metres)
# Out : outputs/<VERSION>/m11_support_scale<SUF>.csv + m11_support_scale<SUF>.png
#   FTM_ONLY <- TRUE; source(here::here("pipeline", "07o_support_scale.R"))   # the delivered set: reads the
#       *_ftmonly run of 07a and writes *_ftm, the curve the dataset documentation quotes. Light, seconds.
#   source(...) alone reads the six-set ladder, the unsuffixed m4a_pred_points table, for both predictor sets.
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here::here("R", "paths.R")); source(here::here("R", "run_utils.R")); source(here::here("R", "eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
SEED <- 7214L
# Two ranges are measured together, and they mean different things.
#   Cluster range, 0.1-2 km. The inventory samples in clusters, and the median nearest-neighbour
#     distance between plots is 50 m (measured and printed in step 1b below). A 100-250 m unit
#     therefore holds two to four subplots of the SAME cluster. Averaging them removes reference
#     noise — positioning error, plot-to-cell misregistration, and the sampling variance of a 0.04 ha
#     plot. This range shows how much of the plot-level error belongs to the reference data.
#     Bell et al. (2023, PNW-GTR-1018) state for gradient nearest-neighbour maps that "accuracy
#     assessments are not appropriate at scales finer than the 0.81-ha or 9-pixel footprint", and
#     0.81 ha is 90 m x 90 m. With min_plots = 1 there is no step at 100 m.
#   Regional range, 5-100 km. The change-of-support curve proper, where units are genuine regional
#     means.
SUPPORTS  <- c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 30, 50, 100)   # km; 0 = plot level, no aggregation
# MIN_PLOTS is a set. The cluster range holds two to four plots per unit by
# design, so demanding five empties it; the regional range can let a sparse unit carry the value. All
# values are computed, with min_plots as an output column.
#   1 counts every plot as its own unit, which holds the sample fixed at every support. This is the
#     headline series.
#   2 is the smallest value for which an average exists.
#   4 is a "complete cluster" sensitivity. Five or more is unusable in the cluster range, because a
#     cluster holds four subplots and a unit would have to catch two clusters.
MIN_PLOTS_SET <- c(1L, 2L, 4L)
MIN_PLOTS <- 2L                              # aggregate_to default
NPERM     <- 200L                            # permutation control repeats

# ---- reportability: bootstrap interval width ----------------------------------------------------
# Rule: do not report a support whose 95% bootstrap CI is wider than CI_WIDTH_MAX.
# Basis: 0.15 is just above the gain in biomass R² from plot level to the regional range; a wider
#   interval cannot resolve that gain. One condition then cuts both pathologies: the coarse end,
#   where 100 km leaves a handful of units, and an over-strict guard in the cluster range.
# ! ci_width is always written out, so the threshold can be changed after the run.
CI_WIDTH_MAX <- 0.15
NBOOT     <- 1000L
TARGETS   <- c("coniferfrac_full", "shannon_full", "Bio_a_ha")   # continuous only
PSETS     <- c("clim_s2med", "full")                             # without / with the forest-type map
# FTM_ONLY: the delivered set on its own complete-case sample, read from 07a's FTM_ONLY run
# (*_ftmonly) and written to *_ftm. The default run reads the six-set ladder, whose common subset
# keeps only the plots complete under every set.
FTM    <- exists("FTM_ONLY") && isTRUE(FTM_ONLY)
IN_SUF <- if (FTM) "_ftmonly" else ""; SUF <- if (FTM) "_ftm" else ""
if (FTM) PSETS <- "full"
outd <- here("outputs", VERSION)

tm <- run_timer("m11_support_scale", version = VERSION, budget_min = 10)

# ---- 1. load: per-plot out-of-fold predictions + coordinates ------------------------------------
tm$stage("load")
pp <- fread(require_file(file.path(outd, sprintf("m4a_pred_points%s.csv", IN_SUF))), colClasses = list(character = "Plot_ID"))
fo <- fread(require_file(file.path(outd, "cv_folds.csv")),        colClasses = list(character = "Plot_ID"))
pp <- pp[target %in% TARGETS & predset %in% PSETS]
pp[, `:=`(obs = as.numeric(obs), pred = as.numeric(pred))]        # stored as character: the file also holds classification targets
DT <- merge(pp, fo[, .(Plot_ID, UTMK_X, UTMK_Y)], by = "Plot_ID")
DT <- DT[is.finite(obs) & is.finite(pred) & is.finite(UTMK_X) & is.finite(UTMK_Y)]
cat(sprintf("loaded %d plot-predictions | targets=%s | predsets=%s\n",
            nrow(DT), paste(TARGETS, collapse = ","), paste(PSETS, collapse = ",")))

# ---- 1b. sample layout ---------------------------------------------------------------------------
# Measuring the nearest-neighbour distance leaves the cluster structure in the output: if these
# numbers change, the cluster range in SUPPORTS has to be reconsidered.
tm$stage("spacing")
XY <- unique(DT[, .(Plot_ID, UTMK_X, UTMK_Y)])
CELLM <- 2000
XY[, `:=`(bx = floor(UTMK_X / CELLM), by = floor(UTMK_Y / CELLM))]
setkey(XY, bx, by)
nnd <- vapply(seq_len(nrow(XY)), function(i) {
  nb <- XY[.(rep(XY$bx[i] + (-1:1), each = 3), rep(XY$by[i] + (-1:1), times = 3)), nomatch = 0L]
  d <- sqrt((nb$UTMK_X - XY$UTMK_X[i])^2 + (nb$UTMK_Y - XY$UTMK_Y[i])^2)
  d <- d[d > 0]; if (length(d)) min(d) else NA_real_
}, numeric(1))
qn <- quantile(nnd, c(0, .10, .50, .90, 1), na.rm = TRUE)
cat(sprintf("plot spacing: nearest-neighbour distance min %.0f | p10 %.0f | median %.0f | p90 %.0f | max %.0f m\n",
            qn[1], qn[2], qn[3], qn[4], qn[5]))
cat(sprintf("  share with a neighbour within 100 m: %.1f %%  -> a 100-250 m cell holds SUBPLOTS OF ONE CLUSTER,\n",
            100 * mean(nnd < 100, na.rm = TRUE)))
cat("  so aggregation below ~2 km removes REFERENCE noise, not regional variance.\n")

# ---- 2. aggregation primitive --------------------------------------------------------------------
# An absolute grid: coordinates divided by L, with the origin NOT anchored to the sample's min(). The
# unit definition is then reproducible and identical across predictor sets and targets.
unit_id <- function(x, y, L_km) paste(floor(x / (L_km * 1000)), floor(y / (L_km * 1000)), sep = "_")

aggregate_to <- function(d, L_km, min_plots = MIN_PLOTS) {
  if (L_km == 0) return(data.table(o = d$obs, p = d$pred, n = 1L))
  a <- d[, .(o = mean(obs), p = mean(pred), n = .N), by = .(u = unit_id(UTMK_X, UTMK_Y, L_km))]
  a[n >= min_plots, .(o, p, n)]
}

mbe <- function(o, p) mean(p - o)

# ---- 3. support curve + permutation control ------------------------------------------------------
tm$stage("curve")
rows <- list()
for (ps in PSETS) for (tg in TARGETS) {
  d <- DT[predset == ps & target == tg]
  if (!nrow(d)) next
  for (L in SUPPORTS) for (mp in MIN_PLOTS_SET) {
    # The plot level is independent of min_plots, and the row is emitted for every series so that a
    # reader filtering on min_plots still has the ladder's anchor. The values are identical by
    # construction.
    a <- aggregate_to(d, L, mp)
    if (nrow(a) < 3L) next
    # bci() fixes its own seed and restores the RNG state, and the permutation control below is
    # seeded here, so the admissible flag is reproducible.
    ci <- bci(a$o, a$p, r2, NBOOT)
    # Control: shuffle the predictions nationally, then run the same aggregation. This is the null
    # distribution of the R2 that aggregation alone produces.
    # ! Run it for the two series that are reporting candidates and name them explicitly. Deriving
    #   the set from the first element of MIN_PLOTS_SET would silently move the control whenever that
    #   set changes. min_plots = 4 is a sensitivity row and needs no null distribution.
    CTL_SET <- c(1L, 2L)
    ctl <- rep(NA_real_, NPERM)
    if (mp %in% CTL_SET) {
      set.seed(SEED)
      ctl <- vapply(seq_len(NPERM), function(i) {
        dd <- copy(d); dd[, pred := sample(pred)]
        aa <- aggregate_to(dd, L, mp)
        if (nrow(aa) < 3L) NA_real_ else r2(aa$o, aa$p)
      }, numeric(1))
    }
    rows[[length(rows) + 1L]] <- data.table(
      version = VERSION, target = tg, predset = ps, support_km = L,
      n_units = nrow(a), n_plots = sum(a$n), mean_plots_per_unit = round(mean(a$n), 1),
      r2 = round(r2(a$o, a$p), 4), ci_lo = round(ci[1], 4), ci_hi = round(ci[2], 4),
      ci_width = round(ci[2] - ci[1], 4),
      admissible = as.integer((ci[2] - ci[1]) <= CI_WIDTH_MAX),   # reportability, per the rule above
      rmse = round(rmse(a$o, a$p), 4), mbe = round(mbe(a$o, a$p), 4),
      ctl_mean = round(mean(ctl, na.rm = TRUE), 4),
      ctl_lo = round(quantile(ctl, .025, na.rm = TRUE), 4),
      ctl_hi = round(quantile(ctl, .975, na.rm = TRUE), 4),
      min_plots = mp, nperm = NPERM)
  }
}
out <- rbindlist(rows)
fwrite(out, file.path(outd, sprintf("m11_support_scale%s.csv", SUF)))

# ---- 4. consistency gate: L=0 must equal the plot-level point in m4a_ladder.csv -------------------
tm$stage("gate")
lad <- fread(require_file(file.path(outd, sprintf("m4a_ladder%s.csv", IN_SUF))))
ref <- lad[scheme == "spatial" & predset %in% PSETS & target %in% TARGETS, .(target, predset, point)]
chk <- merge(out[support_km == 0, .(target, predset, r2)], ref, by = c("target", "predset"))
chk[, delta := round(r2 - point, 4)]
print(chk, row.names = FALSE)
bad <- chk[abs(delta) > 0.002]
if (nrow(bad)) stop("GATE FAILED - L=0 R2 does not match the m4a_ladder point; the join is wrong. See delta above:\n",
                    paste(capture.output(print(bad)), collapse = "\n"), call. = FALSE)
cat("gate OK — plot-level R² reproduces m4a_ladder.csv within 0.002\n")

# ---- 5. figure -----------------------------------------------------------------------------------
tm$stage("figure")
png(file.path(outd, sprintf("m11_support_scale%s.png", SUF)), width = 1500, height = 760, res = 150)
op <- par(mfrow = c(1, length(PSETS)), mar = c(4.6, 4.2, 3.4, 0.8))
# The figure draws one series only, so that the cluster range and the regional range form a single
# continuous curve. The others stay in the CSV under the min_plots column.
FIG <- out[min_plots == MIN_PLOTS_SET[1]]
xs <- SUPPORTS; xs[xs == 0] <- 0.03                      # a slot for the plot level on a log axis
cols <- c(coniferfrac_full = "#1b7837", shannon_full = "#762a83", Bio_a_ha = "#b35806")
for (ps in PSETS) {
  plot(NA, xlim = range(xs), ylim = c(-0.55, 1), log = "x", xaxt = "n",
       xlab = "Aggregation support (km)", ylab = expression(R^2),
       main = ifelse(ps == "clim_s2med", "without the forest-type map", "the delivered predictor set"))
  axis(1, at = xs, labels = ifelse(SUPPORTS == 0, "plot\n(30 m)", SUPPORTS), cex.axis = .8)
  abline(h = 0, col = "grey70", lty = 2)
  # The boundary between the two ranges. Left of it, averaging subplots of one cluster removes
  # reference noise; right of it, units are genuine regional means.
  abline(v = 3, col = "grey55", lty = 4)
  mtext("| within-cluster (reference noise)      between-cluster (regional means) |",
        side = 3, line = 0.2, cex = 0.62, col = "grey35")
  ctl <- FIG[predset == ps, .(lo = min(ctl_lo), hi = max(ctl_hi)), by = support_km][order(support_km)]
  polygon(c(xs, rev(xs)), c(ctl$lo, rev(ctl$hi)), col = "#00000012", border = NA)
  for (tg in TARGETS) {
    d <- FIG[predset == ps & target == tg][order(support_km)]
    if (!nrow(d)) next
    xx <- xs[SUPPORTS %in% d$support_km]
    # Bootstrap intervals.
    segments(xx, d$ci_lo, xx, d$ci_hi, col = cols[[tg]], lwd = 1)
    segments(xx * 0.94, d$ci_lo, xx * 1.06, d$ci_lo, col = cols[[tg]], lwd = 1)
    segments(xx * 0.94, d$ci_hi, xx * 1.06, d$ci_hi, col = cols[[tg]], lwd = 1)
    # Solid line and filled points only where the value is reportable. Beyond that the curve is drawn
    # as a faint dotted line with open points.
    ok <- d$admissible == 1L
    lines(xx, d$r2, col = paste0(cols[[tg]], "55"), lwd = 1.5, lty = 3)
    if (any(ok)) lines(xx[ok], d$r2[ok], col = cols[[tg]], lwd = 2)
    points(xx[ok],  d$r2[ok],  col = cols[[tg]], pch = 16)
    points(xx[!ok], d$r2[!ok], col = cols[[tg]], pch = 1, cex = .9)
  }
  if (ps == PSETS[1]) {
    legend("bottomright", bty = "n", lwd = 2, col = unname(cols[TARGETS]),
      legend = c("coniferous fraction", "tree-species diversity", "aboveground biomass"))
    legend("bottom", bty = "n", cex = .72, pch = c(16, 1), col = "grey25",
      legend = c(sprintf("filled = 95%% CI width <= %.2f (reportable)", CI_WIDTH_MAX),
                 "open = CI too wide to resolve the effect"))
  }
  mtext("grey band = permutation control", side = 3, line = -1.2, cex = 0.7, col = "grey40")
}
par(op); dev.off()

cat(sprintf("\n===== HEADLINE: R² by aggregation level (min_plots = %d) =====\n", MIN_PLOTS_SET[1]))
print(dcast(out[min_plots == MIN_PLOTS_SET[1]], target + predset ~ support_km, value.var = "r2"),
      row.names = FALSE)
cat("\nRMSE in original units (cannot be inflated by aggregation the way R² can):\n")
print(dcast(out[min_plots == MIN_PLOTS_SET[1]], target + predset ~ support_km, value.var = "rmse"),
      row.names = FALSE)
cat(sprintf("\nn_units (min_plots = %d, predset %s):\n", MIN_PLOTS_SET[1], PSETS[1]))
print(dcast(out[min_plots == MIN_PLOTS_SET[1] & predset == PSETS[1]], target ~ support_km,
            value.var = "n_units"), row.names = FALSE)
cat(sprintf("\n95%% bootstrap CI width (report only where <= %.2f):\n", CI_WIDTH_MAX))
print(dcast(out[min_plots == MIN_PLOTS_SET[1]], target + predset ~ support_km, value.var = "ci_width"),
      row.names = FALSE)
cat(sprintf("\nreportable support range per target (min_plots = %d, predset full):\n", MIN_PLOTS_SET[1]))
for (tg in TARGETS) {
  d <- out[min_plots == MIN_PLOTS_SET[1] & predset == "full" & target == tg & admissible == 1L]
  if (!nrow(d)) { cat(sprintf("  %-18s none\n", tg)); next }
  b <- d[which.max(r2)]
  cat(sprintf("  %-18s up to %g km | best R2 %.3f at %g km (CI width %.3f) | plot-level %.3f\n",
              tg, max(d$support_km), b$r2, b$support_km, b$ci_width,
              out[min_plots == MIN_PLOTS_SET[1] & predset == "full" & target == tg &
                    support_km == 0, r2]))
}
cat("\nsensitivity — same curve requiring a complete cluster (>=4 plots per unit):\n")
print(dcast(out[min_plots == 4L], target + predset ~ support_km, value.var = "r2"), row.names = FALSE)
cat(sprintf("\nwrote m11_support_scale%s.csv (%d rows) + m11_support_scale%s.png\n", SUF, nrow(out), SUF))
tm$done(n = nrow(DT), params = sprintf("supports=%s min_plots=%s nperm=%d nboot=%d",
                                       paste(SUPPORTS, collapse = "/"),
                                       paste(MIN_PLOTS_SET, collapse = "+"), NPERM, NBOOT))
