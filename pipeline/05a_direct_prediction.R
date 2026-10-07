# 05a_direct_prediction.R — fit one random forest per attribute and predict wall to wall.
# ============================================================================================
# The first of the two product families. One random forest per attribute, fitted on every plot, then
# predicted over the whole forest extent, with a per-cell uncertainty layer beside each map.
#
# The model is the one stage 07a scores, refit on every plot instead of on a training fold. Engine,
# predictors and parameters are the same in both places.
#
# One pass over the grid produces all four attributes and their uncertainty from a single read of
# each chunk. The chunking, streaming and core-sizing are shared with 05b.
#
# Model
#   ntree = 500, seed 7214, randomForest defaults: regression mtry = p/3 and nodesize = 5,
#   classification mtry = sqrt(p) and nodesize = 1.
#   Predictors: 14 climate and terrain + 42 Sentinel-2 medians + forest type, diameter class,
#               density class and stand origin from the forest-type map = 60.
#   Targets: coniferous fraction, Shannon diversity and biomass as regression; dominant species as
#            classification over the legend in data/species_codes.csv, with the pine races merged
#            and any legend species below MIN_SP plots pooled into the "other" class (see below).
#
# Uncertainty
#   regression     -> the standard deviation across the ensemble, from predict.all = TRUE
#   classification -> one minus the vote margin between the top two classes, so a cell where two
#                     species are nearly tied reads as uncertain
#
# Masking
#   The maps are bounded by the forest-type map: its attributes are predictors, and they are missing
#   outside its surveyed, stocked extent. A random forest cannot accept a missing value, so every
#   predicted cell needs a complete predictor vector.
#
# Inputs
#   data/interim/covar_5179/<preds_clim>.tif                     (14 climate/topo/coord single-band)
#   data/interim/covar_5179_s2/s2_bands/s2band_%02d.tif          (42 single-band S2 medians)
#   data/interim/fmap_30m/<vintage>/{V_VG,V_DC,V_DN,V_FROR}.tif  (forest-type map)
#   outputs/nfi8/{cv_folds.csv, m2_targets.csv}  · data/interim/nfi8/rd_points.csv  · data/interim/nfi8/plot_features_s2.csv
#   data/species_codes.csv                                       (code to species name, the full legend)
# OUTPUT: outputs/nfi8/m8_maps_30m/ftm_{conifer_fraction,dominant_species,biomass,shannon}_30m.tif
#         outputs/nfi8/m8_uncertainty_30m/ftm_{conifer_fraction,biomass,shannon}_sd_30m.tif
#         outputs/nfi8/m8_uncertainty_30m/ftm_dominant_species_uncertainty_30m.tif  (1 - vote margin)
#         These are the names the dataset is published under (R/paths.R::product_file).
#         outputs/nfi8/m8_maps_30m/m8_domsp_legend_used.csv   (the code<->species rows actually trained)
#
# ! These ARE the delivered files, and they are written in place. A re-run replaces them with
#   unfiltered predictions, so 06b_apply_landcover_filter.R has to be run again afterwards.
#
# HEAVY (prediction over the whole forest extent, block-processed) — run in RStudio:
#   source(here::here("pipeline","05a_direct_prediction.R"))
#   (optional levers, same names/behaviour as 05b: IMPUTE_FREE_GB, IMPUTE_CHUNK, IMPUTE_CORES, IMPUTE_CORES_MAX)
# AFTER: pipeline/06b_apply_landcover_filter.R — re-apply the built/barren/water filter.
# ============================================================================================
suppressPackageStartupMessages({
  library(here); library(data.table); library(raster); library(randomForest)
  library(foreach); library(doSNOW)
})
source(here("R", "paths.R"))
NTREE <- 500L; SEED <- 7214L; MIN_SP <- 20L                        # the ntree, seed and rare-species rule 07a scores with

# The single delivered predictor set. "full" is the run parameter R/paths.R keys the published file
# names by (FILE_TIER), and the token stage 07 uses for the same set in its validation tables.
TIER <- "full"
log_time <- function(m) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), m)); flush.console() }
log_time(sprintf("=== direct-RF production (ntree=%d) ===", NTREE))

outd    <- here("outputs", VERSION)
map_dir <- file.path(outd, "m8_maps_30m")                          # the delivered maps
unc_dir <- file.path(outd, "m8_uncertainty_30m")                   # their uncertainty layers
for (d in c(map_dir, unc_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
tname   <- unname(FILE_TIER[TIER])                                 # delivered file prefix; see R/paths.R

# ===== 0. training data — IDENTICAL loading to 07a_predictor_ladder.R / 05b_plot_index_imputation.R ===========
folds <- fread(file.path(outd, "cv_folds.csv"),  colClasses = list(character = "Plot_ID"))
m2    <- fread(file.path(outd, "m2_targets.csv"), encoding = "UTF-8", colClasses = list(character = "Plot_ID"))
rd    <- fread(require_file(here("data", "interim", VERSION, "rd_points.csv")), encoding = "UTF-8",
               colClasses = list(character = "Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses = list(character = "Plot_ID"))
dup <- intersect(setdiff(names(m2), "Plot_ID"), names(rd)); if (length(dup)) rd <- rd[, setdiff(names(rd), dup), with = FALSE]
DT <- Reduce(function(a, b) merge(a, b, by = "Plot_ID", all.x = TRUE),
             list(folds[, .(Plot_ID, role, fold_spatial)], rd, m2, feat))

# ---- predictor set = 07a's "full" (clim + S2 medians + V*/FROR) ; SAME sort() pin as 05b -----
preds_clim <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
s2_med <- sort(intersect(names(DT), S2_BAND_NAMES))    # the 42 declared bands; sort pins the order, see 05b
vstar  <- V_PREDICTORS                                                                   # forest-type map predictors; stand age is not among them
PRED   <- c(preds_clim, s2_med, vstar, "FROR")
stopifnot(length(s2_med) == 42, all(PRED %in% names(DT)))

DT[, V_VG := factor(as.character(V_VG))]; DT[, FROR := factor(as.character(FROR))]       # factor predictors
for (v in c("V_DC","V_DN")) DT[[v]] <- as.numeric(DT[[v]])

# ---- targets: 3 regression + domsp (pine-merged; coded to the published legend) ------------
REG <- c("coniferfrac_full", "shannon_full", "Bio_a_ha")
leg <- fread(here("data","species_codes.csv"), encoding = "UTF-8")
sp2code <- setNames(as.integer(leg$code), leg$species); SP_SET <- unique(leg$species)    # species -> published integer code
DT[, domsp := enc2utf8(as.character(domsp_full))]
DT[domsp %in% c("강원지방소나무","중부지방소나무"), domsp := "소나무"]                    # pine race -> species (idempotent if merged at 04b)
DT[!(domsp %in% SP_SET) & !is.na(domsp) & domsp != "", domsp := "기타"]                  # collapse to the legend's class set
# The rare-species rule the validation applies (07a, 07b, 07d, 07g: fewer than MIN_SP plots over the whole
# table -> 기타) decides the class set here too. A legend species below MIN_SP joins 기타 and its code is left
# unused; codes are never renumbered, so the legend file keeps every row and the trained classes are the
# subset written to m8_domsp_legend_used.csv.
sp_n   <- DT[!is.na(domsp) & domsp != "", .N, by = domsp]
stopifnot(all(sp_n[N >= MIN_SP & domsp != "기타", domsp] %in% SP_SET))                 # a new common species needs a legend code first
drop_sp <- setdiff(SP_SET[SP_SET != "기타"], sp_n[N >= MIN_SP, domsp])
if (length(drop_sp)) { DT[domsp %in% drop_sp, domsp := "기타"]; SP_SET <- setdiff(SP_SET, drop_sp)
  log_time(sprintf("domsp: %d legend class(es) below MIN_SP=%d pooled to 기타: %s", length(drop_sp), MIN_SP,
                   paste(sprintf("%s (code %d)", drop_sp, sp2code[drop_sp]), collapse = ", "))) }
DT[, domsp := factor(domsp, levels = SP_SET)]                                            # train on the legend's class set/order
fwrite(leg[species %in% SP_SET], file.path(map_dir, "m8_domsp_legend_used.csv"), bom = TRUE)   # the classes actually trained
log_time(sprintf("domsp legend = %d classes (published codes, pine-merged)", length(SP_SET)))

# reference = ALL plots (train+test = 100% DELIVERY) with a COMPLETE predictor+target vector (per target).
trDT <- DT[role %in% c("train", "test")]
mask_x <- complete.cases(trDT[, PRED, with = FALSE])
log_time(sprintf("reference pool = %d plots (role in train/test) | predictors = %d",
                 nrow(trDT), length(PRED)))

# ===== 1. train one randomForest per target on ALL plots (100% delivery) ============================
# Each target uses the plots complete in BOTH predictors and that target (the per-target mask, matching
# 07a's per-predset complete.cases). Regression -> numeric y ; domsp -> factor y (legend levels).
fit_one <- function(tg, is_class) {
  ok <- mask_x & complete.cases(trDT[, tg, with = FALSE])
  if (is_class) ok <- ok & !is.na(trDT$domsp)
  x  <- as.data.frame(trDT[ok, PRED, with = FALSE])
  y  <- if (is_class) droplevels(trDT$domsp[ok]) else trDT[[tg]][ok]
  set.seed(SEED)                                                   # delivery refit (no fold) -> single seed, matches 07a's SEED base
  rf <- randomForest(x = x, y = y, ntree = NTREE)
  log_time(sprintf("  trained %-16s n=%5d %s", tg, sum(ok), if (is_class) sprintf("(%d classes)", nlevels(y)) else ""))
  rf
}
log_time("training delivery models (ntree=500, all plots)...")
models <- setNames(lapply(REG, fit_one, is_class = FALSE), REG)
models[["domsp"]] <- fit_one("domsp", is_class = TRUE)
vg_levels   <- levels(trDT$V_VG[mask_x])
fror_levels <- levels(trDT$FROR[mask_x])

# ===== 2. wall-to-wall predictor stack on the covar_5179 grid — IDENTICAL assembly to 05b ============
covar_f <- setNames(file.path(here("data","interim","covar_5179"), paste0(preds_clim, ".tif")), preds_clim)
s2_band_dir <- here("data", "interim", "covar_5179_s2", "s2_bands")                       # single-band S2 — see 03c; one band per file is what makes this readable
s2_f <- if (dir.exists(s2_band_dir)) file.path(s2_band_dir, sprintf("s2band_%02d.tif", seq_len(42))) else {
          s2_f_tiled <- here("data", "interim", "covar_5179_s2", "s2med42_5179_tiled.tif")
          if (file.exists(s2_f_tiled)) s2_f_tiled else here("data", "interim", "covar_5179_s2", "s2med42_5179.tif") }
for (f in c(covar_f, s2_f)) require_file(f)
covar_r <- raster::stack(unname(covar_f)); names(covar_r) <- preds_clim
s2_r <- raster::stack(s2_f); stopifnot(raster::nlayers(s2_r) == 42L)
names(s2_r) <- S2_BAND_NAMES                                # the file's physical band order; see R/paths.R
stopifnot(setequal(names(s2_r), s2_med))
vg_dir  <- here("data", "interim", "fmap_30m", FTM_VINTAGE)
vstar_f <- setNames(file.path(vg_dir, paste0(vstar, ".tif")), vstar)
fror_f  <- file.path(vg_dir, "V_FROR.tif")
for (f in c(vstar_f, fror_f)) require_file(f)
e  <- raster::extent(covar_r); ce <- terra::ext(e[1], e[2], e[3], e[4])                 # EXTEND the 4 small FTM layers up to covar (see 05b)
vg_ext <- file.path(tempdir(), sprintf("vstar_ext_%s.tif", FTM_VINTAGE))
fr_ext <- file.path(tempdir(), sprintf("fror_ext_%s.tif",  FTM_VINTAGE))
terra::extend(terra::rast(unname(vstar_f)), ce, filename = vg_ext, overwrite = TRUE)
terra::extend(terra::rast(fror_f),          ce, filename = fr_ext, overwrite = TRUE)
vstar_r <- raster::stack(vg_ext);  names(vstar_r) <- vstar
fror_r  <- raster::raster(fr_ext); names(fror_r)  <- "FROR"
TD <- raster::stack(covar_r, s2_r, vstar_r, fror_r)[[PRED]]
stopifnot(identical(names(TD), PRED))
nr <- nrow(TD); nc <- ncol(TD)
log_time(sprintf("TD stack ready (%d layers, %d x %d, covar_5179 grid)", raster::nlayers(TD), nr, nc))

# ---- PRE-FLIGHT: verify the raster predictors agree with the training features before the run (see 05b) ---------
if (max(abs(rd$Lon), na.rm = TRUE) > 180) {
  set.seed(1L); chk <- sample(which(mask_x), min(300L, sum(mask_x)))
  xy <- as.matrix(rd[match(trDT$Plot_ID[chk], rd$Plot_ID), .(Lon, Lat)])
  e  <- raster::extract(TD[[s2_med]], xy)
  Xs <- as.data.frame(trDT[chk, s2_med, with = FALSE])
  cors <- vapply(s2_med, function(cn) suppressWarnings(stats::cor(e[, cn], Xs[[cn]], use = "complete.obs")), numeric(1))
  med_cor <- stats::median(cors, na.rm = TRUE)
  log_time(sprintf("PRE-FLIGHT S2 alignment: median band cor = %.3f over %d plots (want > 0.7)", med_cor, length(chk)))
  if (is.na(med_cor) || med_cor < 0.5) stop("S2 raster MISALIGNED with features (cor ", round(med_cor,3), ") — fix before the run.", call. = FALSE)
  # Every predictor, not only S2: the training table and the wall-to-wall stack must hold the same numbers
  # at the plots, so a stale climate or terrain column cannot slip through. Factors are compared as codes.
  ea <- raster::extract(TD, xy); Xa <- as.data.frame(trDT[chk, PRED, with = FALSE])
  dif <- vapply(PRED, function(cn) { a <- suppressWarnings(as.numeric(as.character(ea[, cn])))
                                     b <- suppressWarnings(as.numeric(as.character(Xa[[cn]])))
                                     k <- is.finite(a) & is.finite(b); if (!any(k)) NA_real_ else max(abs(a[k] - b[k])) }, numeric(1))
  bad <- names(dif)[is.finite(dif) & dif > 1e-3 * pmax(1, abs(sapply(Xa, function(x) stats::median(suppressWarnings(as.numeric(as.character(x))), na.rm = TRUE))))]
  log_time(sprintf("PRE-FLIGHT all %d predictors: max |raster - table| per predictor, worst %s = %.4g", length(PRED),
                   names(which.max(dif)), max(dif, na.rm = TRUE)))
  if (length(bad)) stop("training table and raster stack disagree at the plots for: ", paste(bad, collapse = ", "), call. = FALSE)
} else log_time("PRE-FLIGHT skipped (plot coords not in 5179 metres).")

# ===== 3. wall-to-wall prediction — STREAMING (bounded RAM), same harness as 05b ====================
# Per-worker peak ~ 4 RF models + the predict.all spike (chunk forest cells x NTREE x 8 bytes for the 3 REG
# targets). predict.all is the memory driver here (05b's was newtargets); a SMALL CHUNK bounds it. Levers
# (IMPUTE_*) behave exactly as in 05b. Prediction is cheaper per cell than k = 1 matching.
CHUNK <- if (exists("IMPUTE_CHUNK")) as.integer(IMPUTE_CHUNK) else 4L
mdl_gb <- as.numeric(utils::object.size(models)) / 1e9
FOREST_FR <- 0.75; SAFE <- 0.78
phys <- max(1L, parallel::detectCores(logical = FALSE))
free_gb <- if (exists("IMPUTE_FREE_GB")) as.numeric(IMPUTE_FREE_GB) else NA_real_
per_worker_gb <- mdl_gb + CHUNK * nc * FOREST_FR * (length(REG) * NTREE + length(PRED)) * 8 / 1e9   # predict.all spike (REG x NTREE) dominates
CORES_CAP <- if (exists("IMPUTE_CORES_MAX")) as.integer(IMPUTE_CORES_MAX) else 16L
cores <- if (exists("IMPUTE_CORES")) as.integer(IMPUTE_CORES) else
         if (!is.na(free_gb)) max(1L, min(phys, CORES_CAP, as.integer(floor(SAFE * free_gb / per_worker_gb)))) else
         min(CORES_CAP, max(1L, parallel::detectCores() - 1L))
BATCH <- 40L
starts <- seq(1L, nr, by = CHUNK); nbatch <- ceiling(length(starts) / BATCH)
log_time(sprintf("ADAPTIVE plan: cores=%d chunk=%d | models=%.2fGB per-worker≈%.2fGB peak≈%.0fGB / free=%s | %d chunks / %d batches",
                 cores, CHUNK, mdl_gb, per_worker_gb, cores * per_worker_gb,
                 if (is.na(free_gb)) "?(set IMPUTE_FREE_GB to auto-size)" else sprintf("%.0fGB", free_gb), length(starts), nbatch))
log_time("  -> if memory runs short, lower IMPUTE_CHUNK; holding all trees per cell is what costs.")

# outputs: 4 maps + 3 regression SD layers + 1 dominant-species margin (8 rasters)
#
# m2 names a target by its stratum, coniferfrac_full, while the published file carries no stratum,
# because only the whole-stand version is delivered. R/paths.R keys the published names, so the
# translation belongs here rather than in product_file, which must not quietly accept a stratum it
# does not publish.
PUBLISHED <- c(coniferfrac_full = "coniferfrac", shannon_full = "shannon",
               Bio_a_ha = "Bio_a_ha", domsp = "domsp")
stopifnot(all(c(REG, "domsp") %in% names(PUBLISHED)))
OUT_MAP <- setNames(file.path(map_dir, vapply(c(REG, "domsp"),
                    function(v) product_file(TIER, PUBLISHED[[v]], "map"), "")), c(REG, "domsp"))
OUT_UNC <- c(setNames(file.path(unc_dir, vapply(REG,
                    function(v) product_file(TIER, PUBLISHED[[v]], "sd"), "")), REG),
             domsp = file.path(unc_dir, product_file(TIER, "domsp", "uncertainty")))
new_writer <- function(path, dtype = "FLT4S", na = NA) {
  w <- raster::raster(TD); raster::crs(w) <- "EPSG:5179"           # stamp 5179 (covar grid holds params but no authority code)
  raster::writeStart(w, filename = path, overwrite = TRUE, datatype = dtype, NAflag = na,
                     options = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER")) }
# The dominant-species map is a class raster, so it is written as uint8 with nodata 255, the delivered
# convention; 06d then gives it its colour table and MODE overviews. Every other layer is float32/NaN.
W_MAP <- lapply(names(OUT_MAP), function(k) if (k == "domsp") new_writer(OUT_MAP[[k]], "INT1U", 255) else new_writer(OUT_MAP[[k]]))
names(W_MAP) <- names(OUT_MAP)
W_UNC <- lapply(OUT_UNC, new_writer)

predict_chunk <- function(r0) {                                   # worker: models/TD from clusterExport
  n <- min(CHUNK, nr - r0 + 1L)
  m <- raster::getValues(TD, row = r0, nrows = n)                 # (n*nc) x length(PRED), columns = PRED
  N <- nrow(m)
  res <- list(coniferfrac_full = rep(NA_real_, N), shannon_full = rep(NA_real_, N), Bio_a_ha = rep(NA_real_, N),
              domsp = rep(NA_real_, N), coniferfrac_full_sd = rep(NA_real_, N), shannon_full_sd = rep(NA_real_, N),
              Bio_a_ha_sd = rep(NA_real_, N), domsp_uncertainty = rep(NA_real_, N))
  ok <- stats::complete.cases(m)
  if (any(ok)) {
    d <- as.data.frame(m[ok, , drop = FALSE])
    d$V_VG <- factor(as.character(d$V_VG), levels = vg_levels)
    d$FROR <- factor(as.character(d$FROR), levels = fror_levels)
    ok2 <- stats::complete.cases(d)                              # drop cells whose V_VG/FROR code is unseen in training
    if (any(ok2)) {
      dd  <- d[ok2, , drop = FALSE]; pos <- which(ok)[ok2]
      for (tg in c("coniferfrac_full","shannon_full","Bio_a_ha")) {   # REGRESSION: value + tree-to-tree SD
        pa <- predict(models[[tg]], dd, predict.all = TRUE)
        res[[tg]][pos]           <- as.numeric(pa$aggregate)
        ind <- pa$individual; M <- ncol(ind)                          # per-cell tree-to-tree SD (sample, base R; no matrixStats dep)
        res[[paste0(tg,"_sd")]][pos] <- sqrt(pmax(0, (rowMeans(ind^2) - rowMeans(ind)^2)) * (M / (M - 1)))
      }
      pb <- predict(models[["domsp"]], dd, type = "prob")        # CLASSIFICATION: class code + vote margin
      m1 <- max.col(pb, ties.method = "first"); rr <- seq_len(nrow(pb))
      top1 <- pb[cbind(rr, m1)]; pb[cbind(rr, m1)] <- -Inf
      m2i <- max.col(pb, ties.method = "first"); top2 <- pb[cbind(rr, m2i)]
      res[["domsp"]][pos]             <- SP2CODE[colnames(pb)[m1]]  # label -> published integer code
      res[["domsp_uncertainty"]][pos] <- 1 - (top1 - top2)          # 1 - vote margin (higher = more uncertain, matching the SD layers and the deposit)
    }
  }
  res
}

log_time(sprintf("creating SOCK cluster (%d workers, sequential launch)...", cores))
cl <- parallel::makeCluster(cores, type = "SOCK", outfile = "", setup_strategy = "sequential", setup_timeout = 180)
registerDoSNOW(cl)
parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({ library(raster); library(randomForest) })
  raster::rasterOptions(maxmemory = 1e9); NULL })
SP2CODE <- sp2code
parallel::clusterExport(cl, c("models","TD","vg_levels","fror_levels","CHUNK","nr","SP2CODE","REG"))
parallel::clusterExport(cl, "predict_chunk")

pb <- txtProgressBar(min = 0, max = length(starts), style = 3)
done <- 0L; n_forest <- 0.0; perf <- vector("list", nbatch); t0_all <- Sys.time()
for (b in seq_len(nbatch)) {
  tb0 <- Sys.time()
  bi <- starts[((b - 1L) * BATCH + 1L):min(b * BATCH, length(starts))]
  resb <- foreach(r0 = bi, .noexport = c("models","TD","predict_chunk"),
                  .options.snow = list(progress = function(k) setTxtProgressBar(pb, done + k)),
                  .packages = c("raster","randomForest")) %dopar% predict_chunk(r0)
  nf_b <- 0.0
  for (i in seq_along(bi)) {
    r <- resb[[i]]
    for (tg in c(REG, "domsp")) W_MAP[[tg]] <- raster::writeValues(W_MAP[[tg]], r[[tg]], bi[i])
    for (u in names(OUT_UNC))  W_UNC[[u]]  <- raster::writeValues(W_UNC[[u]],  r[[if (u=="domsp") "domsp_uncertainty" else paste0(u,"_sd")]], bi[i])
    nf <- sum(!is.na(r[["coniferfrac_full"]])); if (nf) { n_forest <- n_forest + nf; nf_b <- nf_b + nf }
  }
  done <- done + length(bi); setTxtProgressBar(pb, done); rm(resb); gc(FALSE)
  dt_b <- as.numeric(difftime(Sys.time(), tb0, units = "secs")); cum_s <- as.numeric(difftime(Sys.time(), t0_all, units = "secs"))
  eta_m <- (cum_s / b) * (nbatch - b) / 60; cps <- nf_b / max(dt_b, 1e-6)
  perf[[b]] <- data.frame(batch = b, chunks = length(bi), forest_cells = nf_b, wall_sec = round(dt_b, 1),
                          cells_per_sec = round(cps), cores = cores, chunk = CHUNK, cum_min = round(cum_s/60, 2), eta_min = round(eta_m, 1))
  cat(sprintf("\n[%s] batch %d/%d: %.0fs | %s forest cells | %s cells/s | cum %.1f min | ETA ~%.0f min\n",
              format(Sys.time(), "%H:%M:%S"), b, nbatch, dt_b, format(nf_b, big.mark = ","),
              format(round(cps), big.mark = ","), cum_s/60, eta_m)); flush.console()
}
close(pb); parallel::stopCluster(cl)
for (w in W_MAP) raster::writeStop(w); for (w in W_UNC) raster::writeStop(w)
total_s <- as.numeric(difftime(Sys.time(), t0_all, units = "secs"))
data.table::fwrite(data.table::rbindlist(perf), file.path(outd, sprintf("predict_perf_%s.csv", tname)), bom = TRUE)
log_time(sprintf("DONE: %d chunks, %d cores, %.1f min | forest cells = %s (%.2f Mha)",
                 length(starts), cores, total_s/60, format(n_forest, big.mark = ","), n_forest * 0.09 / 1e6))
cat(sprintf("\nmaps -> %s\n  %s\nNext: source pipeline/06b_apply_landcover_filter.R (re-apply built/barren/water filter), then pipeline/06d_finalize_rasters.py.\n",
            map_dir, paste(basename(unlist(OUT_MAP)), collapse = ", ")))
