# 05b_plot_index_imputation.R — assign every forest cell one matched inventory plot.
# ================================================================================================
# The second product family. yaImpute with method = "randomForest" and k = 1 gives each cell the
# index of the single most similar plot, so the attributes decoded from it in 05c are a combination
# that was actually observed together in the field.
#
# What the matching runs on
#   y = the plot's field forest structure: forest type, diameter class, density class
#   x = 14 climate, terrain and position variables + 42 Sentinel-2 seasonal medians,
#       plus the forest-type map attributes and stand origin in the FTM-added tier
#
# The four published attributes are deliberately not part of y. Matching on them would tune the
# selection to the quantities the accuracy is then reported for. Structure is directly observed and
# discriminates well, whereas diversity and biomass are weak derived targets that would split the
# forest almost at random and pull ecologically dissimilar plots together. Keeping y field-measured
# while x comes from the map also keeps the two tiers non-circular. Stand age is not matched on
# either; it is largely implied by diameter class and travels with the matched plot.
#
# Why the delivered map uses every plot
#   k = 1 imputation is donor-based, so a plot held out of the library can never appear on the map.
#   Holding out 20 % would permanently erase the rare forest types and species that fell in that
#   20 %. The library is therefore all plots, and accuracy is measured separately in 07b under
#   spatial cross-validation and a disjoint test set, which is a conservative bound on this map.
#
#   ! Sampling this map at plot locations does not validate it. Every plot is a donor, and most
#     cells containing a plot match it back to itself, so most of the comparison is a plot against
#     its own record. 07c does sample it, and reports how often that happened alongside the
#     metrics.
#
# Tiers (IMPUTE_TIER, default "full")
#   full = 60 predictors, bounded by the forest-type map   -> index/imputation_ftm.tif
#   base = 56 predictors, no map attributes                -> index/imputation_base.tif
#   The base tier needs the northward-extended covariate grid and a forest-masked Sentinel-2
#   composite, which together let it reach the far north-east where the map does not survey.
#
# In : data/interim/covar_5179_s2/s2med42_5179.tif (03a -> 03b), data/interim/nfi8/rd_points.csv,
#      data/interim/nfi8/plot_features_s2.csv
# Out: outputs/nfi8/index/imputation{_base,_ftm}.tif
#      outputs/nfi8/models/train_reference{_base,_ftm}.csv and the fitted yaImpute object
#
# Then: 05c decodes the index into attribute rasters, and 06b must be re-run because the
#       land-cover filter is lost whenever these rasters are rewritten.
#
# Heavy — k = 1 over every forest cell in the country, processed in blocks. Run in RStudio:
#   IMPUTE_TIER <- "base"                                   # omit for the default "full"
#   source(here::here("pipeline", "05b_plot_index_imputation.R"))
# ============================================================================================
suppressPackageStartupMessages({
  library(here); library(data.table); library(raster); library(yaImpute); library(randomForest)
  library(foreach); library(doSNOW)      # parallel over row-chunks (all cores); gower loaded on workers
})
source(here("R", "paths.R"))
NTREE <- 500L; KNN <- 1L; SEED <- 7214L                     # match 07b FULL (k=1 TreeMap, ntree=500)

# ---- TIER SWITCH (the same two tiers 07b_compare_families.R scores) --------------------------
#   full (default) = 14 climate, terrain and position + 42 Sentinel-2 medians + forest type,
#                    diameter class, density class and stand origin (60) -> index/imputation_ftm.tif
#   base           = the same 56 without the map attributes              -> index/imputation_base.tif
#
#   Each tier also writes index/nearest_reference_distance{,_base}.tif: the random-forest distance
#   between a cell and the plot it matched, one minus the proportion of trees in which the two share
#   a terminal node, so a higher value means a less similar match. yaImpute returns it as $neiDstTrgs.
#   It is a distance, not a proximity, and it is not an area-of-applicability layer — no convex-hull
#   or threshold test is run, so neither of those names applies. See the dataset documentation.
#
#   base predicts wherever the covariates and the forest-masked Sentinel-2 composite are both present,
#   which reaches the far north-east where the forest-type map does not survey, and it needs no
#   national map, so it is the portable configuration. It requires the northward-extended covariate
#   grid from 02c and a forest-masked composite.
#   Set IMPUTE_TIER <- "base" before sourcing; the default keeps the full path.
TIER <- if (exists("IMPUTE_TIER")) IMPUTE_TIER else "full"; stopifnot(TIER %in% c("base", "full"))
log_time <- function(m) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), m)); flush.console() }  # flush: a long non-yielding run otherwise buffers its output and silence looks like a hang
log_time(sprintf("=== imputation TIER = '%s' ===", TIER))

outd    <- here("outputs", VERSION)
mdl_dir <- file.path(outd, "models"); dir.create(mdl_dir, showWarnings = FALSE, recursive = TRUE)
idx_dir <- file.path(outd, "index");  dir.create(idx_dir, showWarnings = FALSE, recursive = TRUE)
tsuf    <- tier_suffix(TIER)                                           # "_base" | "_ftm"; see R/paths.R
out_tif <- file.path(idx_dir, sprintf("imputation%s.tif", tsuf))
dst_tif <- file.path(idx_dir, sprintf("nearest_reference_distance%s.tif", tsuf))   # RF distance to the matched plot

# ===== 0. training data — IDENTICAL loading to pipeline/07b_compare_families.R ======================
folds <- fread(file.path(outd, "cv_folds.csv"),  colClasses = list(character = "Plot_ID"))
m2    <- fread(file.path(outd, "m2_targets.csv"), colClasses = list(character = "Plot_ID"))
rd    <- fread(require_file(here("data", "interim", VERSION, "rd_points.csv")), encoding = "UTF-8",
               colClasses = list(character = "Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses = list(character = "Plot_ID"))
dup <- intersect(setdiff(names(m2), "Plot_ID"), names(rd)); if (length(dup)) rd <- rd[, setdiff(names(rd), dup), with = FALSE]
DT <- Reduce(function(a, b) merge(a, b, by = "Plot_ID", all.x = TRUE),
             list(folds[, .(Plot_ID, role, fold_spatial)], rd, m2, feat))

# ---- predictor set = predictor tier ('full' or 'base') ; y = FIELD structure {FT,DC,DN} (matched-on) --------------
preds_clim <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
s2_med <- sort(intersect(names(DT), S2_BAND_NAMES))                       # the 42 declared bands (R/paths.R)
# sort() PINS the predictor column order. yaImpute's k=1 randomForest samples mtry columns in a
# seed-driven order, so the same predictors in a different order fit a different forest and change
# a large share of the matched plots. Left unsorted, the order follows the input CSV's column order,
# which is regenerated per machine and not byte-identical -> the deposited product would differ by
# machine. Must match the
# identical sort in pipeline/07b_compare_families.R so the product and its validation agree.
vstar  <- V_PREDICTORS                                                    # FTM raster-X (full tier); V_AC DROPPED
PRED   <- if (TIER == "full") c(preds_clim, s2_med, vstar, "FROR") else   # full = 60 vars (07b 'full')
          c(preds_clim, s2_med)                                           # base = 56 vars (07b 'clim_s2med')
# y = FIELD forest structure (NOT the deliverables). fld_V_VG = field FT (~15 groups) via EF-weighted IV (01b);
# fld_V_DC/fld_V_DN = field diameter/crown-density classes. Deliverables ride as passengers in train_reference.
TARGETS <- c("fld_V_VG","fld_V_DC","fld_V_DN")
stopifnot(length(s2_med) == 42, all(PRED %in% names(DT)), all(TARGETS %in% names(DT)))

DT[, V_VG := factor(as.character(V_VG))]                        # FTM raster-X predictors (full tier)
DT[, FROR := factor(as.character(FROR))]
for (v in c("V_DC","V_DN")) DT[[v]] <- as.numeric(DT[[v]])      # raster-X DC/DN (full predictors)
# y = field structure: FT (~15) as factor + DC/DN numeric. No domsp lumping — domsp is a passenger now, not a target.
DT[, fld_V_VG := factor(as.character(fld_V_VG))]
for (v in c("fld_V_DC","fld_V_DN")) DT[[v]] <- as.numeric(DT[[v]])

# reference library = ALL plots (train+test = 100% DELIVERY) with a COMPLETE predictor+target vector.
# Accuracy is estimated separately by 07b (80/20 + spatial K-fold CV); see header. yaImpute needs x & y complete.
trDT <- DT[role %in% c("train", "test")]
ok   <- complete.cases(trDT[, c(PRED, TARGETS), with = FALSE])
tr   <- trDT[ok]
log_time(sprintf("reference (ALL plots 100%%, complete) = %d plots | predictors = %d | targets = %d",
                 nrow(tr), length(PRED), length(TARGETS)))

X.df <- as.data.frame(tr[, PRED, with = FALSE])
Y.df <- as.data.frame(tr[, TARGETS, with = FALSE]); Y.df$fld_V_VG <- droplevels(factor(Y.df$fld_V_VG))
rownames(X.df) <- rownames(Y.df) <- as.character(seq_len(nrow(X.df)))   # plot_idx = 1..N = index raster values
vg_levels   <- if (TIER == "full") levels(X.df$V_VG)  else NULL         # base has no V_VG/FROR predictors
fror_levels <- if (TIER == "full") levels(X.df$FROR)  else NULL

# ---- train_reference for 05c (plot_idx -> the matched plot's attributes); tier-tagged filename ----
# Stamp = the matched plot's FIELD structure (the Y) + PASSENGER deliverables. Column names V_VG/V_DC/V_DN
# are kept (05c's decode expects them), but the VALUES are the FIELD fld_* — NOT the raster-X predictors.
ref <- tr[, .(plot_idx = seq_len(.N), Plot_ID,
              V_VG = as.character(fld_V_VG), V_DC = fld_V_DC, V_DN = fld_V_DN,   # matched field FT/DC/DN (the Y)
              dom_species = domsp_full, dom_ba_share = dom_ba_share_full,        # 04b domsp (BA, 소나무 병합) = direct-map + 07b axis (NOT 01b dom_species)
              Bio_a_ha, Carbon_a_ha, BAPH,                                        # 04b EF-weighted structure
              shannon_full, coniferfrac_full)]                                    # deliverable passengers (validated, not matched)
fwrite(ref, file.path(mdl_dir, sprintf("train_reference%s.csv", tsuf)), bom = TRUE)
log_time(sprintf("train_reference%s.csv saved (%d plots -> field FT/DC/DN + dom_species(04b) + biomass/carbon/BAPH + shannon/coniferfrac)", tsuf, nrow(ref)))

# STAMP_ONLY rewrites train_reference only, skipping training and the wall-to-wall pass. The index
# raster depends on the matching axis and the predictors, not on the passengers, and `tr` is rebuilt
# deterministically from the same inputs, so plot_idx ordering is unchanged and the existing index
# raster still aligns with the refreshed stamp.
#   IMPUTE_STAMP_ONLY <- TRUE; IMPUTE_TIER <- "<tier>"; source(...05b...)   then 05c
if (isTRUE(get0("IMPUTE_STAMP_ONLY", ifnotfound = FALSE))) {
  log_time("STAMP_ONLY=TRUE -> train_reference rewritten; reusing the existing index raster. Run 05c next.")
} else {

# ===== 1. train yaImpute (random-forest matching supervised by field structure {FT,DC,DN}; k = 1) =====
set.seed(SEED + NTREE)
opt_md <- yai(x = X.df, y = Y.df, method = "randomForest", k = KNN, ntree = NTREE)
saveRDS(opt_md, file.path(mdl_dir, sprintf("yai_s2%s_ntree_%d.rds", tsuf, NTREE)))
log_time(sprintf("yaImpute trained (tier=%s, k=%d, ntree=%d, %d predictors)", TIER, KNN, NTREE, length(PRED)))

# ===== 2. wall-to-wall predictor stack on the covar_5179 grid (all layers share this grid) ===
covar_f <- setNames(file.path(here("data","interim","covar_5179"), paste0(preds_clim, ".tif")), preds_clim)
# Reading order matters more than it looks. `raster` is slow at pulling one band out of the 42-band
# composite: the bands are interleaved, so a few-row window ends up reading close to a whole band.
# `terra` handles the same file well, but a terra SpatRaster cannot be sent to SOCK workers, and the
# parallel design needs that. The way out is to keep `raster` and give it single-band files, which it
# reads quickly — hence the split in 03c. The fallbacks below work but are far too slow for a full
# run.
s2_band_dir <- here("data", "interim", "covar_5179_s2", "s2_bands")
s2_f <- if (dir.exists(s2_band_dir)) file.path(s2_band_dir, sprintf("s2band_%02d.tif", seq_len(42))) else {
          s2_f_tiled <- here("data", "interim", "covar_5179_s2", "s2med42_5179_tiled.tif")
          if (file.exists(s2_f_tiled)) s2_f_tiled else here("data", "interim", "covar_5179_s2", "s2med42_5179.tif")
        }
for (f in c(covar_f, s2_f)) require_file(f)

covar_r <- raster::stack(unname(covar_f)); names(covar_r) <- preds_clim
# Label the stack by the file's physical band order, which is what R/paths.R's S2_BAND_NAMES holds,
# and not by s2_med's order. s2_med comes from the feature table, which Earth Engine exported with
# its properties sorted alphabetically; assigning that order here mislabels every physical band and
# scrambles the S2 predictors, which collapsed the imputation once. It was a real bug.
s2_r <- raster::stack(s2_f); stopifnot(raster::nlayers(s2_r) == 42L)
names(s2_r) <- S2_BAND_NAMES
stopifnot(setequal(names(s2_r), s2_med))                              # same 42 semantic bands; TD[[PRED]] reselects by name

if (TIER == "full") {
  vg_dir  <- here("data", "interim", "fmap_30m", FTM_VINTAGE)
  vstar_f <- setNames(file.path(vg_dir, paste0(vstar, ".tif")), vstar)
  fror_f  <- file.path(vg_dir, "V_FROR.tif")                          # file V_FROR.tif -> layer renamed 'FROR'
  for (f in c(vstar_f, fror_f)) require_file(f)
  # Align the two grids by EXTENDING the four small map layers up to the covariate grid, not by
  # cropping the 56 large covariate and Sentinel-2 layers down. Cropping rewrites every one of those
  # layers through per-block temporary files; extending four layers writes through GDAL once. The
  # added north and edge cells get NA for the map attributes and are dropped by complete.cases in the
  # worker, which is correct — the full tier is bounded by the map, and there is nothing to impute
  # where the map does not reach. Extending to the covariate extent also puts both tiers on one grid,
  # so their outputs align cell for cell.
  e  <- raster::extent(covar_r); ce <- terra::ext(e[1], e[2], e[3], e[4])
  vg_ext <- file.path(tempdir(), sprintf("vstar_ext_%s.tif", FTM_VINTAGE))
  fr_ext <- file.path(tempdir(), sprintf("fror_ext_%s.tif",  FTM_VINTAGE))
  terra::extend(terra::rast(unname(vstar_f)), ce, filename = vg_ext, overwrite = TRUE)   # pad V_* (3 layers) with NA
  terra::extend(terra::rast(fror_f),          ce, filename = fr_ext, overwrite = TRUE)   # pad FROR (1 layer) with NA
  vstar_r <- raster::stack(vg_ext);  names(vstar_r) <- vstar
  fror_r  <- raster::raster(fr_ext); names(fror_r)  <- "FROR"
  TD <- raster::stack(covar_r, s2_r, vstar_r, fror_r)[[PRED]]         # 60 layers on the covar grid, in training-X order
} else {
  TD <- raster::stack(covar_r, s2_r)[[PRED]]                          # 56 layers (base = covar + S2, extended grid)
}
stopifnot(identical(names(TD), PRED))
nr <- nrow(TD); nc <- ncol(TD)
log_time(sprintf("TD stack ready (%d layers, %d x %d, covar_5179 grid)", raster::nlayers(TD), nr, nc))

# ---- PRE-FLIGHT: verify the raster predictors ALIGN with the training features, before the run ----
# Extract TD's 42 S2 bands at ~300 train plots, correlate with X.df's S2 by band. A scrambled/misaligned
# stack (the bug above) gives cor≈0; a correct stack gives cor high (temporal-pooling shift keeps it ≳0.7).
if (max(abs(rd$Lon), na.rm = TRUE) > 180) {                          # rd Lon/Lat are EPSG:5179 metres = TD grid
  set.seed(1L); chk <- sample(seq_len(nrow(tr)), min(300L, nrow(tr)))
  xy <- as.matrix(rd[match(tr$Plot_ID[chk], rd$Plot_ID), .(Lon, Lat)])
  e  <- raster::extract(TD[[s2_med]], xy)
  cors <- vapply(s2_med, function(cn) suppressWarnings(stats::cor(e[, cn], X.df[chk, cn], use = "complete.obs")), numeric(1))
  med_cor <- stats::median(cors, na.rm = TRUE)
  log_time(sprintf("PRE-FLIGHT S2 alignment: median band cor(raster, features) = %.3f over %d plots (want > 0.7)", med_cor, length(chk)))
  if (is.na(med_cor) || med_cor < 0.5)
    stop("S2 raster MISALIGNED with training features (median cor ", round(med_cor, 3),
         ") — fix align_s2 band order/names before the expensive run.", call. = FALSE)
} else log_time("PRE-FLIGHT skipped (plot coords not in 5179 metres).")

# ===== 3+4. wall-to-wall k = 1 imputation — STREAMING (bounded RAM) ==========================
# Memory-safe by construction — collecting every chunk before writing would exhaust memory:
#   (a) opt_md/TD are sent to workers ONCE via clusterExport (NOT re-serialised per task by .export);
#   (b) the grid is processed in super-batches of BATCH chunks — each batch is written to disk and freed
#       immediately (rm+gc), so master RAM never accumulates the full result;
#   (c) `cores` is CAPPED: every worker holds a full opt_md copy, so peak worker RAM ~ cores*(opt_md +
#       CHUNK read buffer). The script PRINTS the projected peak below — LOWER cores/CHUNK if it is near
#       your free RAM. `raster` (not terra) backs the stack: a RasterStack serialises cleanly to SOCK.
# Adaptive sizing. The binding constraint is not total memory but the transient spike inside each
# worker: roughly the forest cells in a chunk times the size of the reference library times eight
# bytes, plus one copy of the model per worker. Because that spike scales with the library, a core
# count that fits one configuration can exhaust memory in another. Rather than fix it, the worker
# count is fitted to the free memory from the measured model size, the library size and the grid
# width. A smaller chunk shrinks the spike, which lets more cores run at once, and CPU is rarely the
# limit here. Levers: set IMPUTE_FREE_GB to the free memory in GB to enable the automatic fit, then
# lower IMPUTE_CHUNK to raise cores; or pin IMPUTE_CORES / IMPUTE_CHUNK directly. Per-batch throughput+ETA are
# logged live and saved to models/impute_perf*.csv so you can tune the NEXT run empirically.
CHUNK <- if (exists("IMPUTE_CHUNK")) as.integer(IMPUTE_CHUNK) else 6L     # smaller default -> smaller spike -> more cores
md_gb  <- as.numeric(utils::object.size(opt_md)) / 1e9
n_ref  <- nrow(tr); FOREST_FR <- 0.75; SAFE <- 0.78                       # dense-chunk forest fraction; RAM safety margin
phys   <- max(1L, parallel::detectCores(logical = FALSE))
free_gb <- if (exists("IMPUTE_FREE_GB")) as.numeric(IMPUTE_FREE_GB) else NA_real_   # set it (GB) to enable auto-sizing
per_worker_gb <- md_gb + CHUNK * nc * (n_ref * FOREST_FR + length(PRED)) * 8 / 1e9  # opt_md + newtargets spike + read buf
# CORES_CAP is a launch-safe ceiling, and the binding constraint is worker startup rather than
# memory: a SOCK cluster asked for more than roughly 16 workers can stall while they launch, since
# each one is a separate R process contending for sockets. Memory would often allow more. Override
# with IMPUTE_CORES if your machine launches them cleanly.
CORES_CAP <- if (exists("IMPUTE_CORES_MAX")) as.integer(IMPUTE_CORES_MAX) else 16L
cores <- if (exists("IMPUTE_CORES")) as.integer(IMPUTE_CORES) else
         if (!is.na(free_gb)) max(1L, min(phys, CORES_CAP, as.integer(floor(SAFE * free_gb / per_worker_gb)))) else
         min(CORES_CAP, max(1L, parallel::detectCores() - 1L))           # no IMPUTE_FREE_GB -> spike-safe default
BATCH <- 40L                                                             # chunks written per flush (bounds master RAM)
starts <- seq(1L, nr, by = CHUNK); nbatch <- ceiling(length(starts) / BATCH)
log_time(sprintf("ADAPTIVE plan: cores=%d chunk=%d | opt_md=%.2fGB refs=%d per-worker≈%.2fGB peak≈%.0fGB / free=%s | %d chunks / %d batches",
                 cores, CHUNK, md_gb, n_ref, per_worker_gb, cores * per_worker_gb,
                 if (is.na(free_gb)) "?(set IMPUTE_FREE_GB to auto-size)" else sprintf("%.0fGB", free_gb),
                 length(starts), nbatch))
log_time("  -> if memory runs short, lower IMPUTE_CHUNK; a smaller chunk also lets more cores run.")

log_time(sprintf("creating SOCK cluster (%d workers, sequential launch)... errors out after 180s instead of hanging forever; if it fails, set a smaller IMPUTE_CORES.", cores))
# setup_strategy = "sequential" starts the workers one at a time. The default starts them all at
# once, which can hang before any work begins. setup_timeout turns a stuck launch into an error
# rather than an unattended run that never returns, and outfile = "" lets worker output reach the
# console.
cl <- parallel::makeCluster(cores, type = "SOCK", outfile = "",
                            setup_strategy = "sequential", setup_timeout = 180)
registerDoSNOW(cl)
log_time("SOCK cluster up; exporting model + stack to workers...")
parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({
  library(raster); library(yaImpute); library(randomForest); library(gower) })
  raster::rasterOptions(maxmemory = 1e9); NULL })
parallel::clusterExport(cl, c("opt_md","TD","vg_levels","fror_levels","CHUNK","nr","TIER"))   # ONCE — not per task

impute_chunk <- function(r0) {                                  # runs on worker; opt_md/TD from clusterExport
  n <- min(CHUNK, nr - r0 + 1L)
  m <- raster::getValues(TD, row = r0, nrows = n)               # (n*nc) x length(PRED) matrix, columns = PRED
  out <- rep(NA_real_, nrow(m)); dst <- rep(NA_real_, nrow(m))  # dst = nearest-reference distance to the matched plot
  ok <- stats::complete.cases(m)
  if (any(ok)) {
    d <- as.data.frame(m[ok, , drop = FALSE])
    if (TIER == "full") {                                      # base has no V_VG/FROR predictors -> skip factor recode
      d$V_VG <- factor(as.character(d$V_VG), levels = vg_levels)
      d$FROR <- factor(as.character(d$FROR), levels = fror_levels)
    }
    ok2 <- stats::complete.cases(d)                            # (full) drop cells whose V_VG/FROR code is unseen in training
    if (any(ok2)) {
      dd <- d[ok2, , drop = FALSE]; rownames(dd) <- paste0("T", seq_len(nrow(dd)))
      nt <- yaImpute::newtargets(opt_md, newdata = dd)          # $neiIdsTrgs = matched plot, $neiDstTrgs = its distance
      out[which(ok)[ok2]] <- as.numeric(nt$neiIdsTrgs[, 1])
      dst[which(ok)[ok2]] <- as.numeric(nt$neiDstTrgs[, 1])     # RF distance to the matched plot (1 - terminal-node proximity)
    }
  }
  list(idx = out, dst = dst)                                   # each length n*nc, row-major
}
parallel::clusterExport(cl, "impute_chunk")

out_r <- raster::raster(TD)
raster::crs(out_r) <- "EPSG:5179"                               # stamp EPSG code (covar grid holds 5179 params but no authority code -> would ship as untagged)
out_r <- raster::writeStart(out_r, filename = out_tif, overwrite = TRUE, datatype = "FLT4S",
                            options = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
dst_r <- raster::raster(TD)                                     # parallel nearest-reference distance raster, same grid
raster::crs(dst_r) <- "EPSG:5179"
dst_r <- raster::writeStart(dst_r, filename = dst_tif, overwrite = TRUE, datatype = "FLT4S",
                            options = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
pb <- txtProgressBar(min = 0, max = length(starts), style = 3)
done <- 0L; n_forest <- 0.0; idx_min <- Inf; idx_max <- -Inf
perf <- vector("list", nbatch); t0_all <- Sys.time()            # per-batch throughput/ETA report -> impute_perf*.csv
for (b in seq_len(nbatch)) {
  tb0 <- Sys.time()
  bi <- starts[((b - 1L) * BATCH + 1L):min(b * BATCH, length(starts))]
  resb <- foreach(r0 = bi, .noexport = c("opt_md","TD","impute_chunk"),   # already on workers via clusterExport
                  .options.snow = list(progress = function(k) setTxtProgressBar(pb, done + k)),
                  .packages = c("raster","yaImpute","randomForest","gower")) %dopar% impute_chunk(r0)
  nf_b <- 0.0
  for (i in seq_along(bi)) {
    v <- resb[[i]]$idx; out_r <- raster::writeValues(out_r, v, bi[i])
    dst_r <- raster::writeValues(dst_r, resb[[i]]$dst, bi[i])
    nf <- sum(!is.na(v))
    if (nf) { n_forest <- n_forest + nf; nf_b <- nf_b + nf
              idx_min <- min(idx_min, min(v, na.rm = TRUE)); idx_max <- max(idx_max, max(v, na.rm = TRUE)) }
  }
  done <- done + length(bi); setTxtProgressBar(pb, done)
  rm(resb); gc(FALSE)                                          # <-- release this batch before the next
  dt_b  <- as.numeric(difftime(Sys.time(), tb0,     units = "secs"))
  cum_s <- as.numeric(difftime(Sys.time(), t0_all,  units = "secs"))
  eta_m <- (cum_s / b) * (nbatch - b) / 60                     # linear ETA from mean batch time so far
  cps   <- nf_b / max(dt_b, 1e-6)
  perf[[b]] <- data.frame(batch = b, chunks = length(bi), forest_cells = nf_b, wall_sec = round(dt_b, 1),
                          cells_per_sec = round(cps), cores = cores, chunk = CHUNK,
                          cum_min = round(cum_s / 60, 2), eta_min = round(eta_m, 1))
  cat(sprintf("\n[%s] batch %d/%d: %.0fs | %s forest cells | %s cells/s | cum %.1f min | ETA ~%.0f min\n",
              format(Sys.time(), "%H:%M:%S"), b, nbatch, dt_b, format(nf_b, big.mark = ","),
              format(round(cps), big.mark = ","), cum_s / 60, eta_m))
  flush.console()                                               # show per-batch ETA live (RStudio buffers %dopar% output otherwise)
}
close(pb); parallel::stopCluster(cl)
out_r <- raster::writeStop(out_r)
dst_r <- raster::writeStop(dst_r)
log_time(sprintf("  nearest-reference distance raster -> %s", dst_tif))
total_s <- as.numeric(difftime(Sys.time(), t0_all, units = "secs"))
perf_f  <- file.path(mdl_dir, sprintf("impute_perf%s.csv", tsuf))
data.table::fwrite(data.table::rbindlist(perf), perf_f, bom = TRUE)          # per-batch throughput/ETA log
log_time(sprintf("imputation done (%d chunks x %d rows, %d cores, %d batches, %.1f min) | perf -> %s",
                 length(starts), CHUNK, cores, nbatch, total_s / 60, basename(perf_f)))
log_time(sprintf("  throughput: %s forest cells/s overall (%s /s/core)",
                 format(round(n_forest / max(total_s, 1e-6)), big.mark = ","),
                 format(round(n_forest / max(total_s, 1e-6) / cores), big.mark = ",")))

cat(sprintf("\nDONE (tier=%s) -> %s\n  imputed forest cells = %s (%.2f Mha) | index range = %g..%g (expect 1..%d)\n",
            TIER, out_tif, format(n_forest, big.mark = ","), n_forest * 0.09 / 1e6, idx_min, idx_max, nrow(ref)))
cat(sprintf("Next (keep IMPUTE_TIER<-\"%s\" set):\n", TIER),
    " 1) source pipeline/05c_index_to_attributes.R          (index -> attribute rasters + dom_species legend)\n",
    " 2) source pipeline/06b_apply_landcover_filter.R          (RE-APPLY built/barren/water filter)\n",
    " 3) validate (pipeline/07b_compare_families.R): imputed FT vs field FT (+ vs FTM FT) + passenger biomass/shannon/domsp.\n")

}   # end else (non-STAMP_ONLY: training + wall-to-wall)
