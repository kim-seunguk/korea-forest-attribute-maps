# 07g_uncertainty.R — build the uncertainty layers and test whether they track real error.
# ============================================================================================
# For the direct-prediction model, on the frozen spatial folds, this attaches a per-prediction
# uncertainty and then tests whether it is calibrated: does a higher predicted uncertainty actually
# correspond to a larger error? A layer that fails this test can still rank locations, but it cannot
# be read as an error bar.
#   - regression (3 continuous targets: conifer fraction, AG biomass, Shannon): spatial-OOF pred/obs + uncertainty
#     = SD across RF trees (predict.all); calibration cor(SD, |error|) per target, binned by decile.
#   - classification (dom species):  confidence = vote-probability margin (top1-top2); accuracy by margin decile.
#   Per-plot pred/obs for all three continuous targets feed the predicted-against-observed panels.
#
# FTM_ONLY <- TRUE runs the delivered predictor set, 05a's (climate and terrain, Sentinel-2, the
# forest-type map attributes and stand origin), on its own complete-case sample and writes *_ftm
# files. The delivered uncertainty layers are that model's, so the calibration figure, 07h and 07i
# read the *_ftm files. Without the switch the same set minus the forest-type map is scored, for
# comparison, into the unsuffixed files.
#
# Inputs : outputs/nfi8/{cv_folds.csv,m2_targets.csv} · data/interim/nfi8/rd_points.csv · data/interim/nfi8/plot_features_s2.csv
# Output : outputs/nfi8/m7_uncertainty_points<SUF>.csv (per-plot) + m7_calibration<SUF>.csv (binned)
#   FTM_ONLY <- TRUE; source(here::here("pipeline", "07g_uncertainty.R"))   # the delivered set
#   source(...) alone scores the no-map comparison; PILOT <- TRUE first for a quick shakedown.
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table); library(randomForest); library(foreach); library(doParallel) })
source(here::here("R","paths.R")); source(here::here("R","run_utils.R")); source(here::here("R","eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
if (!exists("PILOT"))   PILOT   <- FALSE
SEED <- 7214L; MIN_SP <- 20L; NBIN <- 10L
cfg <- run_sizing(VERSION, "m7_uncertainty", pilot=PILOT, budget_min=if (PILOT) 20 else 120)
K <- cfg$K; NTREE <- cfg$ntree
outd <- here("outputs", VERSION)
SUF  <- if (PILOT) "_PILOT" else ""                                # a pilot writes its own files

tm <- run_timer("m7_uncertainty", version=VERSION, budget_min=if (PILOT) 20 else 120)
tm$stage("load")
folds <- fread(file.path(outd,"cv_folds.csv"), colClasses=list(character="Plot_ID"))
m2    <- fread(file.path(outd,"m2_targets.csv"), colClasses=list(character="Plot_ID"))
rd    <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding="UTF-8", colClasses=list(character="Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses=list(character="Plot_ID"))
dupc  <- intersect(setdiff(names(m2),"Plot_ID"), names(rd)); if (length(dupc)) rd <- rd[, setdiff(names(rd),dupc), with=FALSE]
DT <- Reduce(function(a,b) merge(a,b,by="Plot_ID",all.x=TRUE),
             list(folds[,.(Plot_ID,role,fold_spatial)], rd, m2, feat))

preds_clim <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
s2_med <- intersect(names(DT), S2_BAND_NAMES)                      # the 42 declared bands (R/paths.R)
if (sum(complete.cases(DT[, ..s2_med])) < 10000L)                  # guard: stale/mid-copy inputs -> abort loudly, not silent small-n
  stop("merged S2 completeness < 10k — inputs look stale/mid-update; Restart R + re-read files.", call. = FALSE)
PREDS <- c(preds_clim, s2_med)                                    # without the forest-type map: the comparison
if (exists("FTM_ONLY") && isTRUE(FTM_ONLY)) {                     # the delivered set, typed as in 05a
  DT[, V_VG := factor(as.character(V_VG))]; DT[, FROR := factor(as.character(FROR))]
  for (v in c("V_DC","V_DN")) DT[[v]] <- as.numeric(DT[[v]])
  PREDS <- c(preds_clim, s2_med, V_PREDICTORS, "FROR"); SUF <- paste0(SUF, "_ftm")
  cat("FTM_ONLY: predictors = climate/terrain + Sentinel-2 + V_VG/V_DC/V_DN/FROR\n")
}
REG_TS <- c("coniferfrac_full", "Bio_a_ha", "shannon_full")       # coniferous fraction, biomass, Shannon H'
CLS_T  <- "domsp_full"
DT[[CLS_T]] <- as.character(DT[[CLS_T]]); rare <- names(which(table(DT[[CLS_T]]) < MIN_SP))
DT[DT[[CLS_T]] %in% rare, (CLS_T) := "기타"]; DT[[CLS_T]] <- factor(DT[[CLS_T]])
fold_sp <- DT$fold_spatial
run_banner("m7_uncertainty", PILOT, cfg, n=nrow(DT),
           extra=sprintf("reg=%s cls=%s", paste(REG_TS, collapse=","), CLS_T))

tm$stage("cluster"); cl <- register_cluster(16L)
# ---- regression: spatial-OOF pred/obs + per-tree SD, for every continuous target --------------
# The three continuous targets shown as predicted against observed in the figure, all on the same
# predictor set so they share identical predictors. The reliability curves are computed for each of
# them.
tm$stage("reg")
reg_list <- list(); reg_cal_list <- list(); reg_cor <- setNames(numeric(0), character(0))
for (RT in REG_TS) {
  mask <- complete.cases(DT[, c(PREDS, RT), with=FALSE]); tr <- which(DT$role=="train" & mask)
  rf_oof <- foreach(f=1:K, .combine=rbind, .packages=c("randomForest","data.table"), .export=c("DT","NTREE","SEED","PREDS","RT","fold_sp")) %dopar% {
    trr <- tr[fold_sp[tr]!=f]; te <- tr[fold_sp[tr]==f]; set.seed(SEED+f)
    rf <- randomForest(x=as.data.frame(DT[trr, PREDS, with=FALSE]), y=DT[[RT]][trr], ntree=NTREE)
    pa <- predict(rf, as.data.frame(DT[te, PREDS, with=FALSE]), predict.all=TRUE)$individual    # n_te x ntree
    data.frame(pos=te, pred=rowMeans(pa), tree_sd=apply(pa,1,sd)) }
  setDT(rf_oof)
  rf_oof[, `:=`(obs=DT[[RT]][pos], Plot_ID=DT$Plot_ID[pos], target=RT)][, abs_err := abs(obs-pred)]
  reg_list[[RT]] <- rf_oof
  reg_cor[RT] <- cor(rf_oof$tree_sd, rf_oof$abs_err)             # calibration for every continuous target
  rb <- cut(rf_oof$tree_sd, unique(quantile(rf_oof$tree_sd, 0:NBIN/NBIN)), include.lowest=TRUE)
  reg_cal_list[[RT]] <- data.table(target=RT, metric="tree_sd", bin=seq_along(levels(rb)),
                        mean_uncert=round(as.numeric(tapply(rf_oof$tree_sd, rb, mean)),4),
                        mean_abs_err=round(as.numeric(tapply(rf_oof$abs_err, rb, mean)),4), n=as.integer(table(rb)))
}
reg     <- rbindlist(reg_list)
reg_cal <- rbindlist(reg_cal_list)

# ---- classification: vote-margin as confidence (spatial OOF) -------------------------------
tm$stage("cls")
mask <- complete.cases(DT[, c(PREDS, CLS_T), with=FALSE]); tr <- which(DT$role=="train" & mask)
cls <- foreach(f=1:K, .combine=rbind, .packages=c("randomForest","data.table"), .export=c("DT","NTREE","SEED","PREDS","CLS_T","fold_sp")) %dopar% {
  trr <- tr[fold_sp[tr]!=f]; te <- tr[fold_sp[tr]==f]; set.seed(SEED+f)
  rf <- randomForest(x=as.data.frame(DT[trr, PREDS, with=FALSE]), y=droplevels(factor(DT[[CLS_T]][trr])), ntree=NTREE)
  pb <- predict(rf, as.data.frame(DT[te, PREDS, with=FALSE]), type="prob")
  m  <- t(apply(pb, 1, function(p){ s<-sort(p,decreasing=TRUE); c(s[1], s[1]-s[2]) }))
  data.frame(pos=te, pred=colnames(pb)[max.col(pb)], top1prob=m[,1], margin=m[,2], stringsAsFactors=FALSE) }
cls$obs <- as.character(DT[[CLS_T]][cls$pos]); cls$correct <- as.integer(cls$pred==cls$obs); cls$Plot_ID <- DT$Plot_ID[cls$pos]
cls_cor <- cor(cls$margin, cls$correct)
cb <- cut(cls$margin, unique(quantile(cls$margin, 0:NBIN/NBIN)), include.lowest=TRUE)
cls_cal <- data.table(target=CLS_T, metric="vote_margin", bin=seq_along(levels(cb)),
                      mean_uncert=round(as.numeric(tapply(cls$margin, cb, mean)),4),
                      accuracy=round(as.numeric(tapply(cls$correct, cb, mean)),4), n=as.integer(table(cb)))
stopCluster(cl)
setDT(reg); setDT(cls)                                            # foreach .combine=rbind returns data.frame; .() needs data.table

# ---- write + calibration verdict -----------------------------------------------------------
fwrite(rbind(
  reg[, .(target, Plot_ID, pred=round(pred,4), obs=round(obs,4), uncertainty=round(tree_sd,4), err=round(abs_err,4))],
  cls[, .(target=CLS_T, Plot_ID, pred, obs, uncertainty=round(margin,4), err=round(1-correct,4))]),
  file.path(outd, sprintf("m7_uncertainty_points%s.csv", SUF)))
fwrite(rbind(reg_cal, cls_cal, fill=TRUE), file.path(outd, sprintf("m7_calibration%s.csv", SUF)))
cat("\nREG cor(tree_SD, |err|) per target (positive => uncertainty tracks error = calibrated):\n")
print(round(reg_cor, 3))
print(reg_cal[, .(target, bin, mean_uncert, mean_abs_err, n)], row.names=FALSE)
cat(sprintf("\nCLS (%s): cor(vote_margin, correct) = %+.3f  (positive => higher margin = more accurate)\n", CLS_T, cls_cor))
print(cls_cal[, .(bin, mean_uncert, accuracy, n)], row.names=FALSE)
cat("\nCalibrated (monotone error-vs-uncertainty) => ship the uncertainty layer with the product.\n")
tm$done(n=nrow(DT), params=sprintf("K=%d ntree=%d", K, NTREE))
