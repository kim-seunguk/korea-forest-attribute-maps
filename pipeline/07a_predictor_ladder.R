# 07a_predictor_ladder.R — score direct prediction for every target and predictor set.
# ================================================================================================
# One random forest per target, fitted independently. This is the first of the two product families
# and the reference the imputation in 07b is compared against.
#
# Each target is scored three ways, so the optimism that a random split would introduce is measured
# on the frozen spatial folds, on random folds over the same training plots, and on the spatially
# disjoint test set. The gap between spatial and random folds is the optimism that random
# partitioning would have introduced. A no-information rate and a bootstrap confidence interval
# accompany every score.
#
# Targets are scored on canopy-only strata as well as on the whole stand, since a satellite sees the
# canopy and not what grows beneath it.
#
# Fitting runs in parallel over the flattened predictor-set, scheme and fold grid, so a worker takes
# one (predictor set, scheme, fold) at a time rather than a whole target.
#
# In : outputs/nfi8/{cv_folds.csv, m2_targets.csv} · data/interim/nfi8/rd_points.csv
#      data/interim/nfi8/plot_features_s2.csv
# Out: outputs/nfi8/m4a_ladder.csv · outputs/nfi8/m4a_pred_points.csv
#
# Heavy. Runs the full pass, which is what the reported numbers come from. For a quick shakedown
# on fewer trees and fewer bootstrap resamples, set PILOT first:
#   source(here::here("pipeline", "07a_predictor_ladder.R"))
#   PILOT <- TRUE; source(here::here("pipeline", "07a_predictor_ladder.R"))
# ================================================================================================
suppressPackageStartupMessages({ library(here); library(data.table); library(randomForest); library(foreach); library(doParallel) })
source(here::here("R","paths.R")); source(here::here("R","run_utils.R")); source(here::here("R","eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
if (!exists("PILOT"))   PILOT   <- FALSE
SEED <- 7214L; MIN_SP <- 20L
cfg   <- run_sizing(VERSION, "m4a_ladder", pilot = PILOT, budget_min = if (PILOT) 30 else 180)
K <- cfg$K; NTREE <- cfg$ntree; NBOOT <- cfg$nboot
outd  <- here("outputs", VERSION); dir.create(outd, showWarnings = FALSE, recursive = TRUE)

# ---- 1. load + merge (frozen folds + the targets + features) --------------------------------
tm <- run_timer("m4a_ladder", version = VERSION, budget_min = if (PILOT) 30 else 180)
tm$stage("load")
folds <- fread(file.path(outd, "cv_folds.csv"), colClasses = list(character = "Plot_ID"))
m2    <- fread(file.path(outd, "m2_targets.csv"), colClasses = list(character = "Plot_ID"))
rd    <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding = "UTF-8",
               colClasses = list(character = "Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses = list(character = "Plot_ID"))
dup <- intersect(setdiff(names(m2), "Plot_ID"), names(rd))          # Bio_a_ha/Carbon_a_ha/BAPH live in BOTH (cor=1.0)
if (length(dup)) rd <- rd[, setdiff(names(rd), dup), with = FALSE]  # keep the 04b versions -> avoid .x/.y suffixing
DT <- Reduce(function(a,b) merge(a, b, by = "Plot_ID", all.x = TRUE),
             list(folds[, .(Plot_ID, role, fold_spatial, fold_random, block_id)], rd, m2, feat))
cat(sprintf("merged %d plots | train=%d test=%d | dropped dup rd cols: %s\n",
            nrow(DT), sum(DT$role=="train"), sum(DT$role=="test"), if (length(dup)) paste(dup, collapse=",") else "none"))

# ---- 2. predictor sets + target typing -----------------------------------------------------
# sin_y and cos_x are the sine and cosine of terrain ASPECT, not of a coordinate — the names are
# historical. Lon and Lat are the coordinates, so they alone come out for the no-coordinate rung.
preds_clim    <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
preds_nolonlat<- setdiff(preds_clim, c("Lon","Lat"))
s2_med  <- intersect(names(DT), S2_BAND_NAMES)          # the 42 declared bands, in the table's own order
vstar   <- intersect(V_PREDICTORS, names(DT))            # forest-type map predictors; stand age is not among them
has_fror<- "FROR" %in% names(DT)
DT[, V_VG := factor(as.character(V_VG))]; for (v in c("V_DC","V_DN")) if (v %in% names(DT)) DT[[v]] <- as.numeric(DT[[v]])
if (has_fror) DT[, FROR := factor(as.character(FROR))]
# clim_s2med is the base predictor set the dataset is delivered under; full is FTM-added. 07d and
# 07b name the same two base and full, and 07e calls the no-map tier none. The token in each output
# is that script's own.
PSETS <- list(
  clim          = preds_clim,
  clim_nolonlat = preds_nolonlat,
  clim_s2med    = c(preds_clim, s2_med),
  clim_vstar    = c(preds_clim, vstar),
  full          = unique(c(preds_clim, s2_med, vstar, if (has_fror) "FROR")))
cat(sprintf("predsets: s2_med=%d vstar=%d fror=%s\n", length(s2_med), length(vstar), has_fror))

REG <- c("coniferfrac_full","coniferfrac_crown","coniferfrac_ba",
         "shannon_full","shannon_crown","shannon_ba","Bio_a_ha")    # regression ladder (composition + diversity + structure)
CLS <- c("domsp_full")                                              # classification (dominant species)
TOL <- c(coniferfrac_full=0.1, coniferfrac_crown=0.1, coniferfrac_ba=0.1,
         shannon_full=0.25, shannon_crown=0.25, shannon_ba=0.25)    # value-tier band (NA -> skip value tier)
for (clt in CLS) {                                                  # rare dominant species -> 기타 (tractable, comparable)
  DT[[clt]] <- as.character(DT[[clt]]); rare <- names(which(table(DT[[clt]]) < MIN_SP))
  DT[DT[[clt]] %in% rare, (clt) := "기타"]; DT[[clt]] <- factor(DT[[clt]])
}
# BASE_ONLY restricts the run to the base predictor set. The completeness mask is built from the
# union of every predictor in play, so dropping the FTM-added set widens the sample from the plots
# complete under both sets (about 12,032) to every plot complete under base alone (about 12,926).
# That is the sample the headline figures are reported on, and it matches 07b, 07e and 07g.
#
# Leave it FALSE when comparing base against FTM-added, where the two must share plots to be
# comparable. The masking line itself is unchanged either way; it simply sees fewer predictors.
if (exists("BASE_ONLY") && isTRUE(BASE_ONLY)) {
  PSETS <- PSETS[c("clim", "clim_s2med")]                            # both are subsets of the base predictors
  REG   <- c("coniferfrac_full", "shannon_full", "Bio_a_ha")         # whole-stand targets only
  cat("BASE_ONLY: predsets =", paste(names(PSETS), collapse = ","), "| targets =", length(REG) + length(CLS), "\n")
}
all_pred <- unique(unlist(PSETS))
run_banner("m4a_ladder", PILOT, cfg, n = nrow(DT), extra = sprintf("targets=%d predsets=%d", length(REG)+length(CLS), length(PSETS)))

# ---- 3. cluster: jobs parallelize over predset x scheme x fold (DT sent via .export per target) ---
tm$stage("cluster")
cl <- register_cluster(16L)                                        # workers inherit master libPaths (fixes .doSnowGlobals)

# ---- 4. run the ladder ---------------------------------------------------------------------
rows <- list()
# Three columns describe accuracy, and they answer different questions.
#   point  exactness at the plot: R2 for a continuous target, top-1 for a class.
#   value  usefulness within a tolerance: within TOL of the truth, or in the top three classes.
#   window aggregated to the 50 km block before comparing, which removes plot-to-cell mismatch and
#          asks whether the map is right about the region even where it is wrong about the plot.
# The paper reports point. The other two are diagnostics and are not comparable across targets.
metric_row <- function(target, stratum, ps, scheme, o, p, in3, cell, is_class, tol) {
  if (is_class) { pt <- mean(p == o); vl <- mean(in3); wn <- NA_real_; nr <- nir(o)
                  ci <- bci(o, p, function(a,b) mean(a == b), NBOOT) }
  else { pt <- r2(o, p); vl <- if (is.na(tol)) NA_real_ else mean(abs(o - p) <= tol)
         ao <- tapply(o, cell, mean); ap <- tapply(p, cell, mean)[names(ao)]
         wn <- if (length(ao) > 2) r2(as.numeric(ao), as.numeric(ap)) else NA_real_
         nr <- NA_real_; ci <- bci(o, p, r2, NBOOT) }
  rows[[length(rows) + 1L]] <<- data.frame(target = target, stratum = stratum, predset = ps, scheme = scheme,
    point = round(pt, 4), value = round(vl, 4), window = round(wn, 4),
    rmse = if (is_class) NA else round(rmse(o, p), 4), nir = round(nr, 4),
    ci_lo = round(ci[1], 4), ci_hi = round(ci[2], 4), n = length(o), testbed = VERSION)
}
run_target <- function(tg, is_class) {
  mask <- complete.cases(DT[, c(all_pred, tg), with = FALSE])
  tr <- which(DT$role == "train" & mask); te <- which(DT$role == "test" & mask)
  fsp <- DT$fold_spatial[tr]; frd <- DT$fold_random[tr]
  stratum <- if (grepl("_crown$", tg)) "crown" else if (grepl("_ba$", tg)) "ba" else if (grepl("_full$", tg)) "full" else "-"
  tol <- if (tg %in% names(TOL)) TOL[[tg]] else NA_real_
  jobs <- list()                                                   # (predset x scheme x fold) + (predset x test)
  for (ps in names(PSETS)) { for (sc in c("spatial","random")) for (f in 1:K) jobs[[length(jobs)+1L]] <- list(ps=ps, sc=sc, f=f)
                             jobs[[length(jobs)+1L]] <- list(ps=ps, sc="test", f=NA_integer_) }
  res <- foreach(j = jobs, .combine = rbind, .packages = c("randomForest","data.table"),
                 .export = c("DT","NTREE","SEED","PSETS")) %dopar% {   # locals (tr/te/fsp/frd/tg/is_class) auto-detected
    preds <- PSETS[[j$ps]]
    if (j$sc == "test") { trr <- tr; pidx <- te } else {
      fo <- if (j$sc == "spatial") fsp else frd; trr <- tr[fo != j$f]; pidx <- tr[fo == j$f] }
    set.seed(SEED + (if (is.na(j$f)) 0L else j$f))
    y <- DT[[tg]][trr]; if (is_class) y <- droplevels(factor(y))
    rf <- randomForest(x = as.data.frame(DT[trr, preds, with = FALSE]), y = y, ntree = NTREE)
    if (is_class) {
      pb <- predict(rf, as.data.frame(DT[pidx, preds, with = FALSE]), type = "prob")
      t3 <- apply(pb, 1, function(p) names(sort(p, decreasing = TRUE))[1:3])
      data.frame(ps = j$ps, sc = j$sc, pos = pidx, pred = colnames(pb)[max.col(pb)],
                 in3 = vapply(seq_along(pidx), function(i) as.character(DT[[tg]][pidx][i]) %in% t3[, i], logical(1)),
                 stringsAsFactors = FALSE)
    } else data.frame(ps = j$ps, sc = j$sc, pos = pidx,
                      pred = as.numeric(predict(rf, as.data.frame(DT[pidx, preds, with = FALSE]))), in3 = NA)
  }
  for (ps in names(PSETS)) for (sc in c("spatial","random","test")) {
    d <- res[res$ps == ps & res$sc == sc, ]; if (!nrow(d)) next; d <- d[order(d$pos), ]
    metric_row(tg, stratum, ps, sc, DT[[tg]][d$pos], d$pred, d$in3, DT$block_id[d$pos], is_class, tol)
  }
  for (ps in FIG_PSETS) {                                            # spatial OOF per plot -> global accumulator (feeds the accuracy figure)
    d <- res[res$ps == ps & res$sc == "spatial", ]; if (!nrow(d)) next
    pred_pts[[length(pred_pts) + 1L]] <<- data.table(
      Plot_ID = DT$Plot_ID[d$pos], predset = ps, target = tg,
      obs = as.character(DT[[tg]][d$pos]), pred = as.character(d$pred))
  }
}
# --- per-plot spatial-CV OOF predictions for the direct-vs-imputation figure (mirrors 07b section 4b) ---
# the predicted-against-observed panels must be the SAME run as the accuracy table, so emit 07a's own out-of-fold predictions here (base + FTM tiers) instead
# of relying on the separate 07g uncertainty run (which gave slightly different R² -> figure/table mismatch).
FIG_PSETS <- intersect(c("clim_s2med", "full"), names(PSETS))
pred_pts  <- list()
allT <- c(REG, CLS)
for (i in seq_along(allT)) {
  tg <- allT[i]; tm$stage(sprintf("t:%s", tg)); progress(i, length(allT), tg)
  run_target(tg, is_class = tg %in% CLS)
  if (i == 1L) estimate_total(as.numeric(difftime(Sys.time(), tm$t_run0, units = "secs")), length(allT),
                              budget_min = if (PILOT) 30 else 180, label = "ladder")
}
stopCluster(cl)

# ---- 5. write + headline summary -----------------------------------------------------------
out <- rbindlist(rows); fwrite(out, file.path(outd, "m4a_ladder.csv"))
if (length(pred_pts)) {                                              # per-plot OOF (route=direct) for fig05_accuracy.R
  pp <- rbindlist(pred_pts)
  fwrite(pp, file.path(outd, "m4a_pred_points.csv"), bom = TRUE)
  cat(sprintf("wrote m4a_pred_points.csv (%d plot-predictions; route=direct, spatial-CV, predset in {%s})\n",
              nrow(pp), paste(FIG_PSETS, collapse = ", ")))
}
cat("\n===== HEADLINE: composition vs Shannon, full vs canopy (spatial folds) =====\n")
hd <- out[scheme == "spatial" & predset %in% c("clim_s2med","full") &
            grepl("coniferfrac|shannon", target), .(target, predset, point, window, ci_lo, ci_hi, n)]
print(hd[order(target, predset)], row.names = FALSE)
cat("\nDoes composition map better than diversity? Compare coniferfrac point vs shannon point.\n")
cat("Does restricting to the canopy help? Compare *_crown / *_ba against *_full.\n")
cat("Optimism: per (target,predset) compare scheme spatial vs random `point` in m4a_ladder.csv.\n")
tm$done(n = nrow(DT), params = sprintf("K=%d ntree=%d nboot=%d targets=%d psets=%d", K, NTREE, NBOOT, length(allT), length(PSETS)))
