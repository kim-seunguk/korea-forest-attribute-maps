# 07d_extra_metrics.R — mean bias and the dominant-species confusion matrix.
# ============================================================================================
# Adds the metrics 07a does not emit: mean bias error for the continuous targets, and the
# dominant-species confusion matrix with per-class producer's and user's accuracy and F1. R² and
# RMSE are repeated as a cross-check against m4a_ladder.csv; %RMSE = 100 * RMSE / mean(observed).
#
# It reuses 07a's data build — the frozen folds, the targets and the stage 03 features — over the
# delivered predictor set and the same set without the forest-type map, on the same spatial-block CV
# and the same independent test set. No bootstrap, so it is light. Set NCORE before sourcing to
# override the core count.
#   source(here::here("pipeline","07d_extra_metrics.R"))
#   FTM_ONLY <- TRUE; source(...)    # the delivered set alone, on its own complete-case sample
#
# Outputs: outputs/nfi8/m4f_continuous_metrics<SUF>.csv, m4f_domsp_confusion_<predset>_<scheme><SUF>.csv,
#          m4f_domsp_perclass<SUF>.csv
#          SUF names the SAMPLE, not the predictor set: "" for the common subset of plots complete
#          under both sets; "_wholeftm" when FTM_ONLY is TRUE, every plot with a complete record under
#          the delivered set, which is what the published per-class table is built from (07k);
#          "_wholebase" when BASE_ONLY is TRUE, every plot complete under the no-map set.
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table); library(randomForest); library(foreach); library(doParallel) })
source(here::here("R","paths.R")); source(here::here("R","eval_utils.R")); source(here::here("R","run_utils.R"))
if (!exists("VERSION")) VERSION <- "nfi8"
if (!exists("NCORE"))   NCORE   <- max(1L, parallel::detectCores() - 1L)   # jobs cap the useful count anyway
SEED <- 7214L; MIN_SP <- 20L; K <- 10L; NTREE <- 500L
outd <- here("outputs", VERSION)

# ---- 1. load + merge (identical to 07a) ----------------------------------------------------
folds <- fread(file.path(outd, "cv_folds.csv"), colClasses = list(character = "Plot_ID"))
m2    <- fread(file.path(outd, "m2_targets.csv"), encoding = "UTF-8", colClasses = list(character = "Plot_ID"))
rd    <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding = "UTF-8",
               colClasses = list(character = "Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses = list(character = "Plot_ID"))
dup <- intersect(setdiff(names(m2), "Plot_ID"), names(rd)); if (length(dup)) rd <- rd[, setdiff(names(rd), dup), with = FALSE]
DT <- Reduce(function(a,b) merge(a, b, by = "Plot_ID", all.x = TRUE),
             list(folds[, .(Plot_ID, role, fold_spatial, block_id)], rd, m2, feat))

# ---- 2. predictor sets + typing (identical to 07a) -----------------------------------------
preds_clim <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
s2_med <- intersect(names(DT), S2_BAND_NAMES)                                                    # the 42 declared bands
vstar <- intersect(V_PREDICTORS, names(DT)); has_fror <- "FROR" %in% names(DT)   # stand age is not among them
DT[, V_VG := factor(as.character(V_VG))]; for (v in c("V_DC","V_DN")) if (v %in% names(DT)) DT[[v]] <- as.numeric(DT[[v]])
if (has_fror) DT[, FROR := factor(as.character(FROR))]
# full is the delivered predictor set, which 07a and 07b also call full; base is the same without the
# forest-type map, which they call clim_s2med.
PSETS <- list(base = c(preds_clim, s2_med), full = unique(c(preds_clim, s2_med, vstar, if (has_fror) "FROR")))
REG <- c("coniferfrac_full","shannon_full","Bio_a_ha")
CLS <- "domsp_full"
DT[[CLS]] <- as.character(DT[[CLS]]); rare <- names(which(table(DT[[CLS]]) < MIN_SP))
DT[DT[[CLS]] %in% rare, (CLS) := "기타"]; DT[[CLS]] <- factor(DT[[CLS]])
# The completeness mask below is built from the union of the sets in play. With both sets it is the
# common subset of plots complete under both; dropping one set widens the sample to every plot
# complete under the other, which is a different sample and so writes its own files (see the header).
SUF <- ""
if (exists("FTM_ONLY") && isTRUE(FTM_ONLY)) { PSETS <- PSETS["full"]; SUF <- "_wholeftm"; cat("FTM_ONLY: predsets = full\n") }
if (exists("BASE_ONLY") && isTRUE(BASE_ONLY)) { PSETS <- PSETS["base"]; SUF <- "_wholebase"; cat("BASE_ONLY: predsets = base\n") }
all_pred <- unique(unlist(PSETS)); allT <- c(REG, CLS)

# ---- 3. flatten jobs (target x predset x {fold 1..K | test}) + run in parallel --------------
jobs <- list()
for (tg in allT) for (ps in names(PSETS)) {
  for (f in 1:K) jobs[[length(jobs)+1L]] <- list(tg = tg, ps = ps, sc = "spatial", f = f)
  jobs[[length(jobs)+1L]] <- list(tg = tg, ps = ps, sc = "test", f = NA_integer_)
}
cat(sprintf("07d: %d jobs (%d targets x %d predsets x (K=%d + test)) on %d cores\n",
            length(jobs), length(allT), length(PSETS), K, NCORE))
cl <- register_cluster(NCORE)
t0 <- Sys.time()
res <- foreach(j = jobs, .combine = rbind, .packages = c("randomForest","data.table"),
               .export = c("DT","NTREE","SEED","PSETS","all_pred","K","CLS")) %dopar% {
  tg <- j$tg; preds <- PSETS[[j$ps]]; is_class <- tg %in% CLS
  mask <- complete.cases(DT[, c(all_pred, tg), with = FALSE])
  tr <- which(DT$role == "train" & mask); te <- which(DT$role == "test" & mask); fsp <- DT$fold_spatial[tr]
  if (j$sc == "test") { trr <- tr; pidx <- te; seed <- SEED } else { trr <- tr[fsp != j$f]; pidx <- tr[fsp == j$f]; seed <- SEED + j$f }
  if (!length(pidx)) return(NULL)
  set.seed(seed); y <- DT[[tg]][trr]; if (is_class) y <- droplevels(factor(y))
  rf <- randomForest(x = as.data.frame(DT[trr, preds, with = FALSE]), y = y, ntree = NTREE)
  p  <- predict(rf, as.data.frame(DT[pidx, preds, with = FALSE]))
  data.frame(tg = tg, ps = j$ps, sc = j$sc, pos = pidx, pred = as.character(p), stringsAsFactors = FALSE)
}
stopCluster(cl)
setDT(res)
cat(sprintf("predictions done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- 4. continuous: R2 / RMSE / %RMSE / MBE (spatial pooled OOF + test) ---------------------
crows <- list()
for (tgt in REG) for (pst in names(PSETS)) for (scm in c("spatial","test")) {
  d <- res[tg == tgt & ps == pst & sc == scm]; if (!nrow(d)) next
  o <- as.numeric(DT[[tgt]][d$pos]); p <- as.numeric(d$pred)
  crows[[length(crows)+1L]] <- data.table(target = tgt, predset = pst, scheme = scm,
    r2 = round(r2(o,p),4), rmse = round(rmse(o,p),4), pct_rmse = round(100*rmse(o,p)/mean(o),1),
    mbe = round(mean(p - o),4), n = length(o))
}
fwrite(rbindlist(crows), file.path(outd, sprintf("m4f_continuous_metrics%s.csv", SUF)))

# ---- 5. dominant species: confusion matrix + per-class PA/UA/F1 -----------------------------
prows <- list()
for (pst in names(PSETS)) for (scm in c("spatial","test")) {
  d <- res[tg == CLS & ps == pst & sc == scm]; if (!nrow(d)) next
  obs <- as.character(DT[[CLS]][d$pos]); prd <- d$pred; lv <- sort(union(unique(obs), unique(prd)))
  cm <- table(obs = factor(obs, lv), pred = factor(prd, lv))
  fwrite(as.data.table(as.data.frame.matrix(cm), keep.rownames = "obs"),
         file.path(outd, sprintf("m4f_domsp_confusion_%s_%s%s.csv", pst, scm, SUF)))
  for (cc in lv) { tp <- cm[cc,cc]; pa <- tp / sum(cm[cc,]); ua <- tp / sum(cm[,cc])
    # F1 = 2TP / (2TP + FP + FN): 0 whenever the class was observed but never hit, even if it was never
    # predicted (user's accuracy undefined). NA only for a class with no observations.
    f1 <- if (sum(cm[cc,]) == 0) NA_real_ else if (tp == 0) 0 else 2*pa*ua/(pa+ua)
    prows[[length(prows)+1L]] <- data.table(predset = pst, scheme = scm, species = cc,
      n_obs = sum(cm[cc,]), producers_acc = round(pa,4), users_acc = round(ua,4), f1 = round(f1,4)) }
  cat(sprintf("domsp %s / %s: OA=%.4f\n", pst, scm, sum(diag(cm))/sum(cm)))
}
fwrite(rbindlist(prows), file.path(outd, sprintf("m4f_domsp_perclass%s.csv", SUF)))
cat(sprintf("\nwrote m4f_continuous_metrics%s.csv ; m4f_domsp_confusion_*_*%s.csv ; m4f_domsp_perclass%s.csv to %s\n",
            SUF, SUF, SUF, outd))
