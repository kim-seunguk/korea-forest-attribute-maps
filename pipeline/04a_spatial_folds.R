# 04a_spatial_folds.R — freeze the validation design. Run this before fitting anything.
# ============================================================================================
# Freezes the one split every later stage reuses:
#   - 50 km blocks. The variogram is computed and reported, but it does not choose the size. At
#     4 km plot spacing the semivariogram is nugget-dominated and its range unreliable, so the
#     block was set where measured leakage stops falling. Both numbers are printed.
#   - an independent test set of whole held-out blocks, spatially disjoint from the training data,
#     plus K spatial folds, and K random folds as the leaky contrast.
# Engine = R/eval_utils.R; instrumentation = R/run_utils.R (stage timing + run_timings.csv).
#
# Inputs : data/interim/nfi8/rd_points.csv  (Plot_ID, UTMK_X/Y, Bio_a_ha, BAPH, dom_ba_share)
# Outputs: outputs/nfi8/cv_folds.csv         (Plot_ID, block_id, role{test,train}, fold_spatial, fold_random, block_km)
#          outputs/nfi8/variogram_range.csv  (target, range_km, sill) + chosen block_km
#          outputs/nfi8/variogram_range.png  (justification plot)
#   source(here::here("pipeline","04a_spatial_folds.R"))
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here::here("R", "paths.R")); source(here::here("R", "run_utils.R")); source(here::here("R", "eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
K <- 10L; TEST_FRAC <- 0.2; SEED <- 7214L
VG_TARGETS <- c("Bio_a_ha", "BAPH", "dom_ba_share")        # continuous proxies for spatial structure (model-free)
BLOCK_CLAMP <- c(15, 80)                                    # km: floor keeps enough plots/block; ceiling keeps enough blocks for K
# The variogram range is biased low by the nugget for plot attributes, so it is not used directly.
# The block size is instead set from where measured leakage stops falling: spatial-CV accuracy for
# biomass reaches its floor well before 50 km, and 50 km still leaves enough blocks for a stable
# K = 10. Set this to NULL to fall back on the variogram's own pick, which is printed either way.
BLOCK_KM_OVERRIDE <- 50L
outd <- here("outputs", VERSION); dir.create(outd, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("\n===== spatial folds [%s]  K=%d  test_frac=%.2f  seed=%d =====\n", VERSION, K, TEST_FRAC, SEED))
tm <- run_timer("m1_eval_scaffold", version = VERSION, budget_min = 15)

# ---- 1. load ------------------------------------------------------------------------------
tm$stage("load")
rd <- fread(require_file(here("data", "interim", VERSION, "rd_points.csv")),
            encoding = "UTF-8", colClasses = list(character = "Plot_ID"))
stopifnot(all(c("UTMK_X", "UTMK_Y", "Plot_ID") %in% names(rd)))
n <- nrow(rd); cat(sprintf("loaded %d plots\n", n))

# ---- 2. variogram block-sizing -------------------------------------------------------------
# The variogram runs on the raw, undetrended response, because that is the spatial structure a
# coordinate-aware model can exploit and so the leakage-relevant scale. The detrended range is
# reported beside it as context: it is the intrinsic residual range, appropriate only for a model
# given no coordinates. Neither of them sets the block; BLOCK_KM_OVERRIDE above does. What the
# variogram would have picked is computed anyway, and printed, so the two can be compared.
tm$stage("variogram")
vg_rows <- list(); vg_curves <- list()
for (tg in VG_TARGETS) {
  if (!tg %in% names(rd)) { cat(sprintf("  skip %s (absent)\n", tg)); next }
  e_raw <- estimate_block_size(rd$UTMK_X, rd$UTMK_Y, rd[[tg]], detrend = FALSE, seed = SEED)
  e_dt  <- estimate_block_size(rd$UTMK_X, rd$UTMK_Y, rd[[tg]], detrend = TRUE,  seed = SEED)
  vg_rows[[tg]]   <- data.frame(target = tg, range_raw_km = round(e_raw$range_km, 1),
                                range_detrended_km = round(e_dt$range_km, 1), sill = signif(e_raw$sill, 4))
  vg_curves[[tg]] <- transform(e_raw$variogram, target = tg)        # plot the RAW curves (the exploitable structure)
  cat(sprintf("  [%-12s] raw range = %.1f km   (detrended %.1f km)\n", tg, e_raw$range_km, e_dt$range_km))
}
vg_tab   <- do.call(rbind, vg_rows)
ranges   <- vg_tab$range_raw_km[is.finite(vg_tab$range_raw_km)]
raw_blk  <- if (length(ranges)) max(ranges) else 50          # >= largest RAW range (NUGGET-BIASED low for plot attrs)
BLOCK_VG <- max(BLOCK_CLAMP[1], min(BLOCK_CLAMP[2], round(raw_blk / 5) * 5))
BLOCK_KM <- if (!is.null(BLOCK_KM_OVERRIDE)) as.integer(BLOCK_KM_OVERRIDE) else BLOCK_VG
cat(sprintf("variogram block = %d km (max RAW range %.1f km; nugget-biased)  ->  CHOSEN BLOCK_KM = %d km%s\n",
            BLOCK_VG, raw_blk, BLOCK_KM, if (!is.null(BLOCK_KM_OVERRIDE)) "  [set by BLOCK_KM_OVERRIDE]" else ""))
vg_tab <- rbind(vg_tab,
                data.frame(target = "block_km_variogram", range_raw_km = BLOCK_VG, range_detrended_km = NA, sill = NA),
                data.frame(target = "CHOSEN_block_km",     range_raw_km = BLOCK_KM, range_detrended_km = NA, sill = NA))
fwrite(vg_tab, file.path(outd, "variogram_range.csv"))

# justification plot (RAW semivariogram curves + chosen block size)
try({
  vgc <- do.call(rbind, vg_curves)
  png(file.path(outd, "variogram_range.png"), width = 1400, height = 1000, res = 180)
  cols <- seq_along(VG_TARGETS)
  plot(NA, xlim = c(0, max(vgc$lag_km, na.rm = TRUE)), ylim = c(0, max(vgc$gamma, na.rm = TRUE) * 1.05),
       xlab = "lag (km)", ylab = "semivariance (raw, centred)", main = sprintf("%s — raw response variogram", VERSION))
  for (i in seq_along(VG_TARGETS)) { tg <- VG_TARGETS[i]; d <- vgc[vgc$target == tg, ]
    if (nrow(d)) lines(d$lag_km, d$gamma, col = i, lwd = 2) }
  abline(v = BLOCK_KM, lty = 2, col = "grey40")
  legend("bottomright", c(VG_TARGETS, sprintf("block = %d km", BLOCK_KM)),
         col = c(cols, "grey40"), lwd = c(rep(2, length(VG_TARGETS)), 1), lty = c(rep(1, length(VG_TARGETS)), 2), bty = "n")
  dev.off()
}, silent = TRUE)

# ---- 3. blocks + test/train + spatial/random folds ----------------------------------------
tm$stage("folds")
block_id <- make_spatial_blocks(rd$UTMK_X, rd$UTMK_Y, BLOCK_KM)
folds    <- assign_folds(block_id, rd$Plot_ID, K = K, test_frac = TEST_FRAC, seed = SEED)
folds$UTMK_X <- rd$UTMK_X; folds$UTMK_Y <- rd$UTMK_Y; folds$block_km <- BLOCK_KM

n_blk   <- length(unique(block_id))
n_test  <- sum(folds$role == "test"); n_train <- sum(folds$role == "train")
n_tblk  <- length(unique(folds$block_id[folds$role == "train"]))
cat(sprintf("blocks=%d (train-blocks=%d)  TEST=%d plots (%.1f%%)  TRAIN=%d plots\n",
            n_blk, n_tblk, n_test, 100 * n_test / n, n_train))
if (n_tblk < K) cat(sprintf("WARN only %d train blocks for K=%d spatial folds — folds will be coarse.\n", n_tblk, K))
cat("spatial fold sizes (train):", paste(table(folds$fold_spatial), collapse = " "), "\n")

# ---- 4. write the frozen split -------------------------------------------------------------
tm$stage("write")
setcolorder(setDT(folds), c("Plot_ID", "UTMK_X", "UTMK_Y", "block_id", "block_km", "role", "fold_spatial", "fold_random"))
fwrite(folds, file.path(outd, "cv_folds.csv"))
cat(sprintf("wrote %s/cv_folds.csv  (the frozen split every later stage reads via load_cv_folds('%s'))\n", outd, VERSION))

tm$done(n = n, params = sprintf("BLOCK_KM=%d K=%d test_frac=%.2f", BLOCK_KM, K, TEST_FRAC))
cat("NEXT: stage 04b builds the targets on these frozen folds (load_cv_folds).\n")
