# 07r_block_size_test.R — find the block size at which spatial cross-validation stops leaking.
# ============================================================================================
# Plot-level forest attributes are nugget-dominated, so the variogram range reflects local noise
# rather than the scale over which a model can interpolate from its neighbours. The operational test
# of "is the block big enough?" is to run a baseline at increasing block sizes and find where the
# spatial-CV score stops falling: the plateau is where the leakage is severed. Two predictor sets
# locate the channel as well, climate and terrain with the coordinates, and the same set without
# them. If the first keeps falling while the second plateaus early, the coordinates are the leak.
# The chosen block is where the curves plateau; 04a_spatial_folds.R holds it as BLOCK_KM_OVERRIDE.
#
# Engine: spatial-block folds from R/eval_utils.R, drawn by inventory cluster as 04a draws them, and
# a random forest per fold. Diagnostic only, no test holdout. Target = biomass, the starkest leakage
# case. The random-fold baseline does not depend on the block size and is fitted once.
#
# Inputs : data/interim/nfi8/rd_points.csv
# Outputs: outputs/nfi8/block_size_test<SUF>.csv  (predset, block_km, n_blocks, spatial_r2, random_r2, gap, n) + .png
#   source(here::here("pipeline", "07r_block_size_test.R"))
#   PILOT <- TRUE first for a quick shakedown on fewer trees; a pilot writes its own file.
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table); library(randomForest); library(foreach); library(doParallel) })
source(here::here("R", "paths.R")); source(here::here("R", "run_utils.R")); source(here::here("R", "eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
if (!exists("PILOT"))   PILOT   <- FALSE
K <- 10L; SEED <- 7214L; TARGET <- "Bio_a_ha"
NTREE  <- if (PILOT) 200L else 500L
BLOCKS <- c(15, 25, 40, 50, 65, 80, 100)                   # km: from the variogram's own range past the chosen block
preds_clim    <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
preds_nolonlat<- setdiff(preds_clim, c("Lon","Lat","sin_y","cos_x"))   # the same no-coordinate rung as 07a
PSETS <- list(clim = preds_clim, clim_nolonlat = preds_nolonlat)
outd  <- here("outputs", VERSION); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
SUF   <- if (PILOT) "_PILOT" else ""                        # a pilot writes its own file

cat(sprintf("\n===== block-size test [%s]  target=%s  K=%d  ntree=%d  blocks={%s} =====\n",
            VERSION, TARGET, K, NTREE, paste(BLOCKS, collapse = ",")))
tm <- run_timer("block_size_test", version = VERSION, budget_min = if (PILOT) 25 else 90)

# ---- load + common-complete mask -----------------------------------------------------------
tm$stage("load")
rd <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding = "UTF-8",
            colClasses = list(character = "Plot_ID"))
needed <- unique(c(preds_clim, TARGET, "UTMK_X","UTMK_Y"))
rd <- rd[complete.cases(rd[, ..needed])]; n <- nrow(rd); cat(sprintf("n=%d (common-complete)\n", n))
# CLUSTER_FOLDS (default TRUE): random folds and blocks by inventory cluster, as 04a draws them, so
# the subplots of one cluster never sit on both sides of a split. FALSE draws them per plot.
if (!exists("CLUSTER_FOLDS")) CLUSTER_FOLDS <- TRUE
cl_id <- if (CLUSTER_FOLDS) cluster_of(rd$Plot_ID) else NULL
if (CLUSTER_FOLDS) {
  uc <- unique(cl_id); set.seed(SEED + 99L); cf <- setNames(sample(rep_len(1:K, length(uc))), uc)
  fold_random <- as.integer(cf[cl_id])                               # random folds: independent of block size
} else { set.seed(SEED + 99L); fold_random <- sample(rep_len(1:K, n)) }

cl <- register_cluster(K)
cv_r2 <- function(fold, preds) {                                     # spatial/random CV R2 for one (fold-vector, predset)
  oof <- foreach(f = 1:K, .combine = rbind, .packages = c("randomForest","data.table"),
                 .export = c("rd","NTREE","SEED","TARGET")) %dopar% {
    tr <- which(fold != f & !is.na(fold)); te <- which(fold == f); if (!length(te)) return(NULL)
    set.seed(SEED + f)
    rf <- randomForest(x = as.data.frame(rd[tr, preds, with = FALSE]), y = rd[[TARGET]][tr], ntree = NTREE)
    data.frame(pos = te, pred = as.numeric(predict(rf, as.data.frame(rd[te, preds, with = FALSE]))))
  }
  oof <- oof[order(oof$pos), ]; r2(rd[[TARGET]][oof$pos], oof$pred)
}

# ---- random baseline once (block-independent) ----------------------------------------------
tm$stage("random")
rnd <- sapply(names(PSETS), function(ps) cv_r2(fold_random, PSETS[[ps]]))
cat(sprintf("random-CV R2:  %s\n", paste(sprintf("%s=%.3f", names(rnd), rnd), collapse = "  ")))

# ---- spatial CV across block sizes ---------------------------------------------------------
rows <- list(); t_blk0 <- NULL
for (bi in seq_along(BLOCKS)) {
  bkm <- BLOCKS[bi]; tm$stage(sprintf("block_%dkm", bkm))
  blk <- make_spatial_blocks(rd$UTMK_X, rd$UTMK_Y, bkm, cluster = cl_id)
  set.seed(SEED + 1L); bf <- setNames(sample(rep_len(1:K, length(unique(blk)))), unique(blk))
  fold_sp <- as.integer(bf[blk]); nblk <- length(unique(blk))
  for (ps in names(PSETS)) {
    sp <- cv_r2(fold_sp, PSETS[[ps]])
    rows[[length(rows) + 1]] <- data.frame(predset = ps, block_km = bkm, n_blocks = nblk,
      spatial_r2 = round(sp, 4), random_r2 = round(rnd[[ps]], 4), gap = round(rnd[[ps]] - sp, 4), n = n)
    progress(bi, length(BLOCKS), sprintf("%-13s block=%3dkm (%3d blk)  spatial=%.3f  gap=%.3f", ps, bkm, nblk, sp, rnd[[ps]] - sp))
  }
  if (bi == 1) { t_blk0 <- as.numeric(difftime(Sys.time(), tm$t_run0, units = "secs"))
                 estimate_total(t_blk0, length(BLOCKS), budget_min = if (PILOT) 25 else 90, label = "block sweep") }
}
stopCluster(cl)

out <- rbindlist(rows); fwrite(out, file.path(outd, sprintf("block_size_test%s.csv", SUF)))
cat("\n===== block-size sweep (spatial_r2 should DROP then PLATEAU; pick the plateau block) =====\n")
print(out[, .(predset, block_km, n_blocks, spatial_r2, random_r2, gap)], row.names = FALSE)

# plateau plot
try({
  png(file.path(outd, sprintf("block_size_test%s.png", SUF)), width = 1500, height = 1000, res = 180)
  plot(NA, xlim = range(BLOCKS), ylim = range(c(out$spatial_r2, out$random_r2)),
       xlab = "spatial block size (km)", ylab = sprintf("CV R2 (%s)", TARGET),
       main = sprintf("%s: leakage against block size", VERSION))
  for (i in seq_along(PSETS)) { ps <- names(PSETS)[i]; d <- out[predset == ps]
    lines(d$block_km, d$spatial_r2, col = i, lwd = 2, type = "b", pch = 19)
    abline(h = d$random_r2[1], col = i, lty = 3) }
  legend("topright", c(paste0(names(PSETS), " (spatial)"), "random (dotted)"),
         col = c(seq_along(PSETS), "grey40"), lwd = 2, lty = c(rep(1, length(PSETS)), 3), bty = "n")
  dev.off()
}, silent = TRUE)

tm$done(n = n, params = sprintf("target=%s K=%d ntree=%d blocks=%d", TARGET, K, NTREE, length(BLOCKS)))
cat("READ: where each curve flattens the leakage is severed; that is the block to set as BLOCK_KM_OVERRIDE in 04a_spatial_folds.R.\n")
cat("      if clim keeps dropping but clim_nolonlat plateaus early, the coordinates are the leak channel.\n")
