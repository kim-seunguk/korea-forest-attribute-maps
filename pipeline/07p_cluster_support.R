# 07p_cluster_support.R — what changes if the map is built at the inventory cluster's support.
# ============================================================================================
# The delivered product is 30 m. 07o aggregates PREDICTIONS into windows after the fact. This stage
# asks what a product rebuilt at the cluster's support would look like. Rebuilding means
# retraining, which changes three things at once (listed below).
#
# Cluster support follows the sampling design: a cluster holds up to four subplots about 50 m apart,
# so it covers a footprint of roughly 100 m, larger than the 30 m cell. Building at cluster support
# requires no new rasters: averaging the subplots' predictor vectors is a spatial average of the
# predictors over the footprint, sampled at the subplot positions inside it.
#
# Configs 1 and 3 alone would confound two separate effects, so config 2 is added to separate them:
#
#   1. train plot,    evaluate plot      the delivered product
#   2. train plot,    evaluate cluster   keeps every training plot, averages the subplots'
#                                        out-of-fold predictions inside each cluster
#   3. train cluster, evaluate cluster   one training unit and one prediction per cluster
#
# Config 2 isolates the reference-noise gain with no loss of sample, and it is what a product built
# from the same plots at the cluster's support would actually be.
#
# ! Config 2 costs the imputation its defining property. Four subplots draw four different donors,
#   so the cluster value is the mean of four real plots and behaves like an imputation with k = 4.
#   The variance ratio shows that loss.
#   Config 3 preserves k = 1 by matching one donor per cluster, but pays for it in sample size.
#
# The three effects, which pull against each other:
#   + the target becomes a multi-subplot mean, so reference noise falls         -> R2 up
#   - the sample falls from plots to clusters                                    -> R2 down
#   ? the target's own variance falls, and R2 is measured against that variance -> up if the
#     within-cluster variation was unpredictable noise, down if it was signal the map captured
#
# ! Targets are recomputed from the pooled stem list at both supports, with the same code.
#   Averaging plot-level summaries would be wrong for two of the four: conifer fraction is a
#   basal-area-weighted ratio, so the mean of ratios is not the pooled ratio, and Shannon over
#   pooled stems is strictly higher than the mean of subplot Shannons. The plot-level values are
#   checked against targets.csv (the GATE line below).
#
# Both routes are run at every configuration on the frozen 50 km folds. Subplots of one cluster
# share a block and a fold by construction (04a assigns them by cluster); any cluster that still
# straddles a fold is dropped rather than reassigned, because config 2 would otherwise average an
# out-of-fold prediction together with an in-fold one. The count is printed.
#
# FTM_ONLY <- TRUE runs the delivered predictor set and writes *_ftm files; without the switch the
# same set minus the forest-type map is run, for comparison, into the unsuffixed files.
#
# Inputs : data/interim/nfi8/{records.csv, rd_points.csv, plot_features_s2.csv}
#          outputs/nfi8/{cv_folds.csv, targets.csv}
# Output : outputs/nfi8/cluster_support<SUF>.csv   accuracy at both supports, both routes
#          outputs/nfi8/target_variance<SUF>.csv   the target's variance at each support
#          outputs/nfi8/cluster_support<SUF>.png
#   FTM_ONLY <- TRUE; source(here::here("pipeline", "07p_cluster_support.R"))
#   Runs sequentially: the stem table is large and the fits are few, so a worker pool buys little.
# ============================================================================================
suppressPackageStartupMessages({
  library(here); library(data.table); library(randomForest); library(yaImpute); library(foreach)
})
source(here::here("R","paths.R")); source(here::here("R","run_utils.R")); source(here::here("R","eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
SEED   <- 7214L
NTREE  <- 500L
MIN_SP <- 20L                                      # the rare-species rule 07a and 07b apply
FTM <- exists("FTM_ONLY") && isTRUE(FTM_ONLY)
SUF <- if (FTM) "_ftm" else ""
outd <- here("outputs", VERSION)
tm   <- run_timer("07p_cluster_support", version = VERSION, budget_min = 60)

# ---- 1. targets recomputed from the stem list, at plot and at cluster support --------------------
tm$stage("targets")
rec <- fread(require_file(here("data","interim",VERSION,"records.csv")), encoding = "UTF-8",
             colClasses = list(character = "신표본점번호"))
setnames(rec, c("신표본점번호","흉고직경","수종명","침활구분"), c("Plot_ID","DBH","SP","CB"))
rec[, Cluster := cluster_of(Plot_ID)]              # digits 1-12 of the plot id; see R/eval_utils.R
rec <- rec[is.finite(DBH) & DBH >= 6 & is.finite(EF) & EF > 0 & !is.na(SP) & SP != ""]
rec[, ba := pi * (DBH / 200)^2 * EF]               # basal area per hectare contributed by this stem
rec[SP %in% c("강원지방소나무", "중부지방소나무"), SP := "소나무"]   # the two regional pine races are one species, as in 04b

targets_by <- function(key, r = rec) {
  d <- r[, .(ba = sum(ba)), by = c(key, "SP")]
  d[, tot := sum(ba), by = key]
  sh  <- d[, .(shannon = -sum((ba / tot) * log(ba / tot))), by = key]
  dom <- d[d[, .I[which.max(ba)], by = key]$V1][, .(domsp = SP), by = key]
  cf  <- r[, .(coniferfrac = sum(ba[grepl("침", CB)]) / sum(ba)), by = key]
  Reduce(function(a, b) merge(a, b, by = key), list(sh, dom, cf))
}
TP <- targets_by("Plot_ID")
# Rare dominant species are pooled exactly as 07a and 07b pool them, so the accuracy stays
# comparable with every other stage.
rare <- names(which(table(TP$domsp) < MIN_SP))
TP[domsp %in% rare, domsp := "기타"]
cat(sprintf("dominant-species classes after pooling the rare ones (< %d plots): %d\n", MIN_SP, uniqueN(TP$domsp)))

# biomass is a per-hectare quantity, so the cluster value is the mean of its subplots
targets <- fread(file.path(outd, "targets.csv"), colClasses = list(character = "Plot_ID"))
TP <- merge(TP, targets[, .(Plot_ID, Bio_a_ha)], by = "Plot_ID")

# the recomputation must reproduce 04b's plot-level quantities, or the two supports are not comparable
chk <- merge(TP, targets[, .(Plot_ID, cf_targets = coniferfrac_full, sh_targets = shannon_full)], by = "Plot_ID")
cat(sprintf("GATE: recomputed vs targets.csv at plot level -- conifer r = %.4f, shannon r = %.4f (n = %d)\n",
            cor(chk$coniferfrac, chk$cf_targets, use = "complete.obs"),
            cor(chk$shannon, chk$sh_targets, use = "complete.obs"), nrow(chk)))

# ---- 2. predictors: the plot's own vector, and the cluster mean of its subplots' vectors ---------
tm$stage("predictors")
folds <- fread(file.path(outd, "cv_folds.csv"), colClasses = list(character = "Plot_ID"))
rd    <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding = "UTF-8",
               colClasses = list(character = "Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses = list(character = "Plot_ID"))
# rd_points carries its own copies of the targets (Bio_a_ha, dom_species, BAPH ...). TP supplies the
# recomputed ones, so rd's copies are dropped to keep the merge free of .x / .y duplicate columns.
rd <- rd[, setdiff(names(rd), c("Bio_a_ha","Carbon_a_ha","dom_species","dom_ba_share","BAPH",
                                names(TP)[names(TP) != "Plot_ID"])), with = FALSE]
D <- Reduce(function(a, b) merge(a, b, by = "Plot_ID"),
            list(folds[, .(Plot_ID, role, fold_spatial)], rd, feat, TP))
D[, Cluster := cluster_of(Plot_ID)]

preds_clim <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO",
                "sin_y","cos_x","C_PAR","Lon","Lat")
s2_med <- sort(intersect(names(D), S2_BAND_NAMES))   # the 42 declared bands; sort pins the order, see 05b
PRED   <- unique(c(preds_clim, s2_med))              # without the forest-type map: both supports use the same set
if (FTM) {
  # Add the forest-type map columns of the delivered set, typed as in 05a.
  D[, V_VG := factor(as.character(V_VG))]; D[, FROR := factor(as.character(FROR))]
  for (v in c("V_DC","V_DN")) D[[v]] <- as.numeric(D[[v]])
  PRED <- unique(c(PRED, V_PREDICTORS, "FROR"))
  cat("FTM_ONLY: predictors = climate/terrain + Sentinel-2 + V_VG/V_DC/V_DN/FROR\n")
}
MATCH <- c("fld_V_VG","fld_V_DC","fld_V_DN")         # the imputation's matching axis, as in 05b and 07b
D[, fld_V_VG := factor(as.character(fld_V_VG))]
for (v in c("fld_V_DC","fld_V_DN")) D[[v]] <- as.numeric(D[[v]])
D <- D[complete.cases(D[, c(PRED, MATCH), with = FALSE]) & !is.na(role) & role == "train"]
D <- D[is.finite(Bio_a_ha) & is.finite(coniferfrac) & is.finite(shannon) & !is.na(domsp)]
D[, domsp := factor(domsp)]

strad <- D[, .(nf = uniqueN(fold_spatial)), by = Cluster][nf > 1, Cluster]
cat(sprintf("dropping %d of %d clusters that straddle a spatial fold (%.2f%%)\n",
            length(strad), uniqueN(D$Cluster), 100 * length(strad) / uniqueN(D$Cluster)))
D <- D[!Cluster %in% strad]

# The cluster target pools the SAME subplots the prediction side keeps. Pooled from every subplot in
# records.csv it would include subplots with no predictors, outside the folds or in the test set, so
# configs 2 and 3 would score a different quantity than the one predicted.
TC <- targets_by("Cluster", rec[Plot_ID %in% D$Plot_ID])
TC[domsp %in% rare, domsp := "기타"]
TC <- merge(TC, D[, .(Bio_a_ha = mean(Bio_a_ha, na.rm = TRUE)), by = Cluster], by = "Cluster")
cat(sprintf("cluster targets from %d kept subplots in %d clusters\n", nrow(D), nrow(TC)))

# cluster level: mean predictor vector (the modal value for a categorical predictor), the cluster
# target, and the fold, which is identical within a cluster
PRED_NUM <- PRED[vapply(PRED, function(v) is.numeric(D[[v]]), logical(1))]
PRED_CAT <- setdiff(PRED, PRED_NUM)            # empty without the forest-type map
modal <- function(x) names(which.max(table(as.character(x))))
CL <- D[, c(lapply(.SD[, PRED_NUM, with = FALSE], mean), lapply(.SD[, PRED_CAT, with = FALSE], modal),
            .(n_sub = .N, fold_spatial = fold_spatial[1])),
        by = Cluster, .SDcols = PRED]
for (v in PRED_CAT) CL[[v]] <- factor(CL[[v]], levels = levels(D[[v]]))

CL <- merge(CL, TC, by = "Cluster")
CL <- CL[complete.cases(CL[, PRED, with = FALSE]) &
         is.finite(Bio_a_ha) & is.finite(coniferfrac) & is.finite(shannon) & !is.na(domsp)]
CL[, domsp := factor(domsp)]
# The imputation matches on the field structure, so CL takes those columns from the first subplot.
CL <- merge(CL, D[, .SD[1], by = Cluster, .SDcols = MATCH], by = "Cluster")
CL[, fld_V_VG := factor(as.character(fld_V_VG))]
stopifnot(all(MATCH %in% names(CL)), !anyNA(CL[[MATCH[1]]]))
cat(sprintf("plot support: %d plots | cluster support: %d clusters (mean %.2f subplots)\n",
            nrow(D), nrow(CL), mean(CL$n_sub)))

# ---- 3. out-of-fold predictions at PLOT support, fitted once and scored two ways ----------------
tm$stage("fit_plot")
REG <- c("coniferfrac","Bio_a_ha","shannon"); CLS <- "domsp"
FO  <- sort(unique(D$fold_spatial))
registerDoSEQ()
rm(list = intersect(ls(), c("rec", "feat", "rd", "targets", "folds", "chk", "TP", "strad")))
invisible(gc())

jobs <- expand.grid(tg = c(REG, CLS), f = FO, stringsAsFactors = FALSE)
cat(sprintf("plot support: %d jobs (%d targets x %d folds), sequential\n",
            nrow(jobs), length(c(REG, CLS)), length(FO)))

oof_list <- foreach(j = seq_len(nrow(jobs)), .packages = c("randomForest","data.table")) %dopar% {
  tg <- jobs$tg[j]; f <- jobs$f[j]
  tr <- which(D$fold_spatial != f); te <- which(D$fold_spatial == f)
  set.seed(SEED + f)
  y  <- if (tg %in% REG) D[[tg]][tr] else droplevels(D[[tg]][tr])
  rf <- randomForest(x = as.data.frame(D[tr, PRED, with = FALSE]), y = y, ntree = NTREE)
  p  <- predict(rf, as.data.frame(D[te, PRED, with = FALSE]))
  cat(sprintf("    [%s fold %d] %s\n", tg, f, format(Sys.time(), "%H:%M:%S"))); flush.console()
  list(tg = tg, te = te, pred = if (tg %in% REG) as.numeric(p) else as.character(p))
}
OOF <- list()
for (tg in c(REG, CLS)) OOF[[tg]] <- if (tg %in% REG) rep(NA_real_, nrow(D)) else rep(NA_character_, nrow(D))
for (o in oof_list) OOF[[o$tg]][o$te] <- o$pred
stopifnot(!any(sapply(OOF, function(v) any(is.na(v)))))

# the imputation donor, also once: one yai fit per fold
tm$stage("fit_impute")
don_list <- foreach(f = FO, .packages = c("yaImpute","data.table")) %dopar% {
  tr <- which(D$fold_spatial != f); te <- which(D$fold_spatial == f)
  set.seed(SEED + f)
  xr <- as.data.frame(D[tr, PRED, with = FALSE]); rownames(xr) <- as.character(tr)
  yr <- as.data.frame(D[tr, MATCH, with = FALSE]); rownames(yr) <- as.character(tr)
  yr$fld_V_VG <- droplevels(factor(yr$fld_V_VG))
  xp <- as.data.frame(D[te, PRED, with = FALSE]); rownames(xp) <- as.character(te)
  nt <- newtargets(yai(x = xr, y = yr, method = "randomForest", k = 1L, ntree = NTREE), newdata = xp)
  list(te = te, don = as.integer(nt$neiIdsTrgs[, 1]))
}
DON <- rep(NA_integer_, nrow(D))
for (o in don_list) DON[o$te] <- o$don
stopifnot(!anyNA(DON), !any(DON == seq_along(DON)))

# ---- 4. score the same predictions at plot support and at cluster support ------------------------
tm$stage("score")
row_reg <- function(sup, rt, tg, obs, p) data.table(
  support = sup, route = rt, target = tg, n = length(obs),
  point = round(r2(obs, p), 4), rmse = round(rmse(obs, p), 4),
  var_obs = round(var(obs), 4), varratio = round(var(p) / var(obs), 4))
row_cls <- function(sup, rt, tg, obs, p) data.table(
  support = sup, route = rt, target = tg, n = length(obs),
  point = round(mean(as.character(p) == as.character(obs)), 4), rmse = NA_real_,
  var_obs = NA_real_, varratio = NA_real_)

res <- list()
for (tg in c(REG, CLS)) {
  for (rt in c("direct","imputation")) {
    p <- if (rt == "direct") OOF[[tg]] else D[[tg]][DON]
    # CONFIG 1: plot support
    res[[length(res)+1L]] <- if (tg %in% REG) row_reg("1_plot", rt, tg, D[[tg]], p)
                             else            row_cls("1_plot", rt, tg, D[[tg]], p)
    # CONFIG 2: the same predictions, averaged inside each cluster
    agg <- if (tg %in% REG) D[, .(pred = mean(p[.I])), by = Cluster]
           else D[, .(pred = names(sort(table(p[.I]), decreasing = TRUE))[1]), by = Cluster]
    mm  <- merge(agg, TC[, c("Cluster", tg), with = FALSE], by = "Cluster")
    res[[length(res)+1L]] <- if (tg %in% REG) row_reg("2_plot_eval_cluster", rt, tg, mm[[tg]], mm$pred)
                             else            row_cls("2_plot_eval_cluster", rt, tg, mm[[tg]], mm$pred)
  }
}

# Configs 1 and 2 are written before config 3 runs, so a failure there does not discard them.
fwrite(rbindlist(res), file.path(outd, sprintf("cluster_support_partial%s.csv", SUF)))
cat(sprintf("configs 1-2 written to cluster_support_partial%s.csv\n", SUF)); flush.console()

# ---- 5. CONFIG 3: refit at cluster support (fewer, cleaner training units) ----------------------
tm$stage("fit_cluster")
FOC <- sort(unique(CL$fold_spatial))
jobs3 <- expand.grid(tg = c(REG, CLS), f = FOC, stringsAsFactors = FALSE)
oof3 <- foreach(j = seq_len(nrow(jobs3)), .packages = c("randomForest","data.table")) %dopar% {
  tg <- jobs3$tg[j]; f <- jobs3$f[j]
  tr <- which(CL$fold_spatial != f); te <- which(CL$fold_spatial == f)
  set.seed(SEED + f)
  y  <- if (tg %in% REG) CL[[tg]][tr] else droplevels(CL[[tg]][tr])
  rf <- randomForest(x = as.data.frame(CL[tr, PRED, with = FALSE]), y = y, ntree = NTREE)
  p  <- predict(rf, as.data.frame(CL[te, PRED, with = FALSE]))
  list(tg = tg, te = te, pred = if (tg %in% REG) as.numeric(p) else as.character(p))
}
O3 <- list()
for (tg in c(REG, CLS)) O3[[tg]] <- if (tg %in% REG) rep(NA_real_, nrow(CL)) else rep(NA_character_, nrow(CL))
for (o in oof3) O3[[o$tg]][o$te] <- o$pred

don3_list <- foreach(f = FOC, .packages = c("yaImpute","data.table")) %dopar% {
  tr <- which(CL$fold_spatial != f); te <- which(CL$fold_spatial == f)
  set.seed(SEED + f)
  xr <- as.data.frame(CL[tr, PRED, with = FALSE]); rownames(xr) <- as.character(tr)
  yr <- as.data.frame(CL[tr, MATCH, with = FALSE]); rownames(yr) <- as.character(tr)
  yr$fld_V_VG <- droplevels(factor(yr$fld_V_VG))
  xp <- as.data.frame(CL[te, PRED, with = FALSE]); rownames(xp) <- as.character(te)
  nt <- newtargets(yai(x = xr, y = yr, method = "randomForest", k = 1L, ntree = NTREE), newdata = xp)
  list(te = te, don = as.integer(nt$neiIdsTrgs[, 1]))
}
DON3 <- rep(NA_integer_, nrow(CL)); for (o in don3_list) DON3[o$te] <- o$don

for (tg in c(REG, CLS)) for (rt in c("direct","imputation")) {
  p <- if (rt == "direct") O3[[tg]] else CL[[tg]][DON3]
  res[[length(res)+1L]] <- if (tg %in% REG) row_reg("3_cluster", rt, tg, CL[[tg]], p)
                           else            row_cls("3_cluster", rt, tg, CL[[tg]], p)
}

RES <- rbindlist(res)
fwrite(RES, file.path(outd, sprintf("cluster_support%s.csv", SUF)))

cat(sprintf("\n===== accuracy by configuration (%s, 50 km spatial-block CV) =====\n",
            if (FTM) "the delivered predictor set" else "without the forest-type map"))
print(dcast(RES, target + route ~ support, value.var = "point"), row.names = FALSE)
cat("\n-- variance ratio: where config 2 costs the imputation its k=1 property --\n")
print(dcast(RES[target %in% REG], target + route ~ support, value.var = "varratio"), row.names = FALSE)
cat("\n-- the target's own variance at each support --\n")
v <- unique(RES[route == "direct" & target %in% REG, .(target, support, var_obs)])
print(dcast(v, target ~ support, value.var = "var_obs"), row.names = FALSE)
fwrite(v, file.path(outd, sprintf("target_variance%s.csv", SUF)))

# ---- 6. figure -----------------------------------------------------------------------------------
tm$stage("figure")
png(file.path(outd, sprintf("cluster_support%s.png", SUF)), width = 1900, height = 900, res = 170)
op <- par(mfrow = c(1, 2), mar = c(7.5, 4.4, 2.6, 0.8), mgp = c(2.5, 0.7, 0))
w <- dcast(RES, target + route ~ support, value.var = "point")
m <- t(as.matrix(w[, .(`1_plot`, `2_plot_eval_cluster`, `3_cluster`)]))
colnames(m) <- paste(w$target, ifelse(w$route == "direct", "dir", "imp"), sep = "\n")
barplot(m, beside = TRUE, col = c("grey35","#2166ac","#b2182b"), border = NA, las = 2,
        cex.names = 0.7, ylab = "spatial-CV R2 / accuracy", main = "accuracy by support configuration")
abline(h = 0); legend("topright", c("1 plot","2 plot-train, cluster-eval","3 cluster"),
                      fill = c("grey35","#2166ac","#b2182b"), border = NA, bty = "n", cex = 0.72)
wv <- dcast(RES[target %in% REG], target + route ~ support, value.var = "varratio")
mv <- t(as.matrix(wv[, .(`1_plot`, `2_plot_eval_cluster`, `3_cluster`)]))
colnames(mv) <- paste(wv$target, ifelse(wv$route == "direct", "dir", "imp"), sep = "\n")
barplot(mv, beside = TRUE, col = c("grey35","#2166ac","#b2182b"), border = NA, las = 2,
        cex.names = 0.7, ylab = "Var(pred)/Var(obs)", main = "variance retained")
abline(h = 1, lty = 2)
par(op); dev.off()

tm$done(n = nrow(D), params = sprintf("clusters=%d ntree=%d preds=%d", nrow(CL), NTREE, length(PRED)))
cat(sprintf("\nwrote cluster_support%s.csv / target_variance%s.csv / cluster_support%s.png\n", SUF, SUF, SUF))
