# R/eval_utils.R — the validation scaffold every stage scores against.
#
# Design choices:
#   (1) spatial-block CV, plus an independent test set of whole held-out blocks that is spatially
#       disjoint from the training data, plus random folds as the leaky contrast.
#   (2) block size is sized against the spatial-autocorrelation range of the raw response, which is
#       what a coordinate-aware model can exploit. 04a may override the variogram's own pick; it
#       prints both.
#   (3) blocks, random folds and bootstrap resamples are drawn by inventory cluster, not by plot. The
#       four subplots of a cluster sit about 50 m apart and are not independent of one another.
#
# This file is the single home for: variogram block-sizing, the test/train + spatial/random fold split, and the
# metric primitives (nir / bootstrap CI / r2 / rmse). Source AFTER R/paths.R (+ R/run_utils.R).
# Base R only. The runnable freezer that calls these and writes cv_folds.csv is pipeline/04a_spatial_folds.R.

# ---- 1. block size from the spatial-autocorrelation range ----------------------------------
# The block must be at least the autocorrelation range the model can exploit. The predictor set
# includes the plot coordinates and smoothly varying climate, so the model can interpolate spatially,
# and the leakage-relevant quantity is the raw response autocorrelation range. That is the default.
# detrend = TRUE removes the 2-D spatial trend first and reports the intrinsic residual range, which
# is smaller and appropriate only for a model given no coordinates; it is kept as context, not to
# choose the block. Range read off a sampled empirical (semi)variogram.
estimate_block_size <- function(x, y, values, detrend = FALSE, n_sample = 3000L, n_bins = 40L,
                                max_dist_frac = 0.5, sill_frac = 0.95, seed = 7214L) {
  ok <- is.finite(values) & is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; v <- values[ok]
  if (detrend) {
    fit <- try(lm(v ~ poly(x, 2) + poly(y, 2)), silent = TRUE)      # intrinsic residual range (context only)
    res <- if (inherits(fit, "try-error")) v - mean(v) else as.numeric(residuals(fit))
  } else res <- v - mean(v)                                         # RAW range = what a coord-aware model exploits
  set.seed(seed)
  idx <- if (length(res) > n_sample) sort(sample.int(length(res), n_sample)) else seq_along(res)
  d_km  <- as.numeric(dist(cbind(x[idx], y[idx]))) / 1000           # pairwise distance (UTMK metres -> km)
  gamma <- 0.5 * as.numeric(dist(res[idx]))^2                        # per-pair semivariance 0.5*(ri-rj)^2
  max_d <- as.numeric(quantile(d_km, max_dist_frac))
  keep  <- d_km <= max_d
  br    <- seq(0, max_d, length.out = n_bins + 1L)
  lag   <- (head(br, -1) + tail(br, -1)) / 2
  vg    <- tapply(gamma[keep], cut(d_km[keep], br, include.lowest = TRUE), mean)
  vg_df <- data.frame(lag_km = lag, gamma = as.numeric(vg))
  sill  <- var(res)                                                 # total residual variance ~ sill
  reach <- which(vg_df$gamma >= sill_frac * sill & is.finite(vg_df$gamma))
  range_km <- if (length(reach)) vg_df$lag_km[reach[1]] else max_d  # practical range
  list(range_km = range_km, sill = sill, variogram = vg_df)
}

# ---- 2. spatial blocks + test/train + spatial/random folds ---------------------------------
# cluster: placed each by its own coordinate, the subplots of a cluster can straddle a block edge and
# so end up on both sides of the train/test split. With `cluster` given, every subplot takes its
# cluster's mean position, so a cluster lies in exactly one block. The grid origin is still taken
# from the plots themselves.
make_spatial_blocks <- function(ux, uy, block_km, cluster = NULL) {
  ox <- min(ux); oy <- min(uy)
  if (!is.null(cluster)) { ux <- ave(ux, cluster); uy <- ave(uy, cluster) }
  bx <- floor((ux - ox) / (block_km * 1000)); by <- floor((uy - oy) / (block_km * 1000))
  paste(bx, by, sep = "_")
}

# Hold out WHOLE blocks as the independent TEST set (spatially disjoint from train, ~test_frac of plots);
# partition the remaining (train) blocks into K SPATIAL folds; assign train plots to K RANDOM folds (the
# leaky contrast that exposes the optimism). Test rows get NA folds (never used in CV).
# cluster: the random folds are assigned per cluster, so that subplots 50 m apart never sit on both
# sides of a random-CV split. Drawn per plot, the random-vs-spatial gap would carry that sibling
# leakage. NULL keeps the per-plot assignment.
assign_folds <- function(block_id, plot_id, K = 10L, test_frac = 0.2, seed = 7214L, cluster = NULL) {
  ublk  <- unique(block_id)
  blk_n <- tapply(seq_along(block_id), block_id, length)            # plots per block
  set.seed(seed); blk_order <- sample(names(blk_n))
  cumf  <- cumsum(blk_n[blk_order]) / length(block_id)
  test_blocks  <- blk_order[cumf <= test_frac]
  if (!length(test_blocks)) test_blocks <- blk_order[1]
  role  <- ifelse(block_id %in% test_blocks, "test", "train")
  train_blocks <- setdiff(ublk, test_blocks)
  set.seed(seed + 1L); blk_fold <- setNames(sample(rep_len(1:K, length(train_blocks))), train_blocks)
  fold_spatial <- as.integer(ifelse(role == "train", blk_fold[block_id], NA_integer_))
  tr <- which(role == "train"); fold_random <- rep(NA_integer_, length(plot_id))
  if (is.null(cluster)) {
    set.seed(seed + 2L); fold_random[tr] <- sample(rep_len(1:K, length(tr)))
  } else {
    uc <- unique(cluster[tr]); set.seed(seed + 2L); cf <- setNames(sample(rep_len(1:K, length(uc))), uc)
    fold_random[tr] <- cf[cluster[tr]]
  }
  data.frame(Plot_ID = plot_id, block_id = block_id, role = role,
             fold_spatial = fold_spatial, fold_random = as.integer(fold_random),
             stringsAsFactors = FALSE)
}

# ---- 3. metric primitives — one copy, so every stage scores the same way -------------------
r2   <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
rmse <- function(o, p) sqrt(mean((o - p)^2))
nir  <- function(o) max(table(o)) / length(o)                       # no-information rate (majority class)
# Bootstrap CI of a metric. cluster: whole clusters are resampled rather than plots, because the
# residuals of the subplots in one cluster are correlated and plot resampling understates the
# interval. A fixed seed makes the interval reproducible; the caller's RNG state is restored on exit.
bci  <- function(o, p, fn, nboot = 1000L, cluster = NULL, seed = 7214L) {
  if (exists(".Random.seed", envir = globalenv())) {
    old <- get(".Random.seed", envir = globalenv()); on.exit(assign(".Random.seed", old, envir = globalenv()))
  }
  set.seed(seed)
  if (is.null(cluster)) {
    v <- replicate(nboot, { i <- sample.int(length(o), replace = TRUE); fn(o[i], p[i]) })
  } else {
    idx <- split(seq_along(o), cluster)
    v <- replicate(nboot, { i <- unlist(idx[sample.int(length(idx), replace = TRUE)], use.names = FALSE); fn(o[i], p[i]) })
  }
  stats::quantile(v, c(.025, .975), na.rm = TRUE)
}
# The cluster a plot belongs to: digits 1-12 of the 13-digit 신표본점번호, the last digit being the subplot.
cluster_of <- function(plot_id) substr(as.character(plot_id), 1L, 12L)

# Read the frozen split, so every stage shares identical folds. Call it once, at the top of a stage.
load_cv_folds <- function(version) {
  f <- here::here("outputs", version, "cv_folds.csv")
  if (!file.exists(f)) stop("cv_folds.csv missing — run pipeline/04a_spatial_folds.R first: ", f, call. = FALSE)
  utils::read.csv(f, colClasses = c(Plot_ID = "character"), stringsAsFactors = FALSE)
}
