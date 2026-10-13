# 07b_compare_families.R — score the plot-index imputation on the same folds as direct prediction.
# ================================================================================================
# The dataset ships two product families, and they are only comparable if both are scored the same
# way. Stage 07a scores the direct random forests; this scores the imputation, reusing the
# frozen spatial folds and the same disjoint test set, and writing rows in the same schema.
#
# Engine: one yaImpute model per predictor set, scheme and fold. yai(method = "randomForest", k = 1)
# matches each held-out plot to its single most similar training plot in random-forest proximity space.
#
# The matching axis is the plot's field forest structure — forest type, diameter class and density
# class. The four published attributes are not matched on; they are read off whichever plot the
# structure selected, as in the delivered product.
#
# Rows are also emitted split by stand origin, natural against planted, mirroring the same split on
# the direct route.
#
# Predictor sets : full = the delivered set (Sentinel-2, climate and terrain, forest-type map
#                  attributes and stand origin) · clim_s2med = the same without the forest-type map.
# Scored : forest type, diameter class and density class (top-1, the axis the matching used) ·
#          dominant species, pooled to the direct map's legend and also unpooled (domsp_raw, what the
#          delivered index decodes to) · coniferous fraction, Shannon diversity and biomass (R2) ·
#          crown density at the one class boundary the field and the map share (fld_V_DN2) ·
#          forest-type-map baselines on the same axes.
#
# In : outputs/nfi8/{cv_folds.csv, m2_targets.csv} · data/interim/nfi8/rd_points.csv
#      data/interim/nfi8/plot_features_s2.csv
#      outputs/nfi8/m4a_ladder_ftmonly.csv and m4a_ladder_baseonly.csv, if 07a has run, for the
#      direct-against-imputation table (the six-set m4a_ladder.csv is read when neither exists)
# Out: outputs/nfi8/m4c_dualtrack.csv · outputs/nfi8/m4c_pred_points.csv
#
# Heavy. Runs the full pass. For a quick shakedown on fewer trees and fewer bootstrap resamples,
# set PILOT first; a pilot writes *_PILOT files.
#   source(here::here("pipeline", "07b_compare_families.R"))
#   PILOT <- TRUE; source(here::here("pipeline", "07b_compare_families.R"))
# ================================================================================================
suppressPackageStartupMessages({
  library(here); library(data.table); library(randomForest); library(yaImpute)
  library(foreach); library(doParallel)
})
source(here::here("R","paths.R")); source(here::here("R","run_utils.R")); source(here::here("R","eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
if (!exists("PILOT"))   PILOT   <- FALSE
SEED <- 7214L; MIN_SP <- 20L; KNN <- 1L                              # k=1 = TreeMap (single matched plot)
cfg  <- run_sizing(VERSION, "m4c_dualtrack", pilot = PILOT, budget_min = if (PILOT) 40 else 300)
K <- cfg$K; NTREE <- cfg$ntree; NBOOT <- cfg$nboot
outd <- here("outputs", VERSION); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
SUF <- if (PILOT) "_PILOT" else ""                                  # a pilot never overwrites the full run's files

# ---- 1. load + merge (identical to 07a) ----------------------------------------------------
tm <- run_timer("m4c_dualtrack", version = VERSION, budget_min = if (PILOT) 40 else 300)
tm$stage("load")
folds <- fread(file.path(outd, "cv_folds.csv"),  colClasses = list(character = "Plot_ID"))
m2    <- fread(file.path(outd, "m2_targets.csv"), colClasses = list(character = "Plot_ID"))
rd    <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding = "UTF-8",
               colClasses = list(character = "Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses = list(character = "Plot_ID"))
dup <- intersect(setdiff(names(m2), "Plot_ID"), names(rd))
if (length(dup)) rd <- rd[, setdiff(names(rd), dup), with = FALSE]
DT <- Reduce(function(a,b) merge(a, b, by = "Plot_ID", all.x = TRUE),
             list(folds[, .(Plot_ID, role, fold_spatial, fold_random, block_id)], rd, m2, feat))
cat(sprintf("merged %d plots | train=%d test=%d\n", nrow(DT), sum(DT$role=="train"), sum(DT$role=="test")))

# ---- 2. predictor sets + target typing (identical to 07a) ----------------------------------
preds_clim <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
s2_med <- sort(intersect(names(DT), S2_BAND_NAMES))   # sort() pins the column order: the yaImpute
# k=1 randomForest match samples mtry columns in a seed-driven order, so the same predictors in a
# different column order fit a different forest.
# s2_med's order otherwise inherits the input CSV's column order, which is regenerated per machine
# and not guaranteed identical between machines, so the scores would drift with it. Must match
# pipeline/05b_plot_index_imputation.R, or the validation and the deposited product would order predictors
# differently.
vstar  <- intersect(V_PREDICTORS, names(DT))            # map predictors; stand age is not among them
has_fror <- "FROR" %in% names(DT)
DT[, V_VG := factor(as.character(V_VG))]; for (v in c("V_DC","V_DN")) if (v %in% names(DT)) DT[[v]] <- as.numeric(DT[[v]])
# matching axis = FIELD structure {FT, DC, DN} (identical typing to pipeline/05b_plot_index_imputation.R)
DT[, fld_V_VG := factor(as.character(fld_V_VG))]; for (v in c("fld_V_DC","fld_V_DN")) if (v %in% names(DT)) DT[[v]] <- as.numeric(DT[[v]])
# yaImpute matches on the NUMERIC ordinals above; the SCORING below reads them as classes, so keep a
# character copy rather than re-coercing DT (which would change the matching axis).
if (has_fror) {
  DT[, FROR := factor(as.character(FROR))]
  # Stand-origin strata for the imputation route, mirroring 07f. Codes:
  # 1 = 인공림 planted, 2 = 천연림 natural, 0 = 무립목지/비산림 (non-forest). The script prints the distribution.
  DT[, FROR_lab := c("0"="other","1"="planted","2"="natural")[as.character(FROR)]]
  cat("\n-- FROR distribution (1=planted / 2=natural / 0=non-forest; mirrors 07f) --\n"); print(DT[, .N, by = FROR][order(FROR)])
  cat(sprintf("labelled: planted=%d natural=%d other/NA=%d\n",
              sum(DT$FROR_lab=="planted",na.rm=TRUE), sum(DT$FROR_lab=="natural",na.rm=TRUE),
              sum(is.na(DT$FROR_lab)|DT$FROR_lab=="other")))
}
TIERS <- list(clim_s2med = c(preds_clim, s2_med),                   # without the forest-type map
              full       = unique(c(preds_clim, s2_med, vstar, if (has_fror) "FROR")))   # the delivered set
cat(sprintf("predictor sets: clim_s2med=%d vars | full=%d vars (s2_med=%d vstar=%d fror=%s)\n",
            length(TIERS$clim_s2med), length(TIERS$full), length(s2_med), length(vstar), has_fror))

# MATCH = yaImpute matching axis (field structure); TARGETS = the SCORED set read off the matched plot.
MATCH <- c("fld_V_VG","fld_V_DC","fld_V_DN")                       # field FT + diameter/density classes (== 05b)
REG   <- c("coniferfrac_full","shannon_full","Bio_a_ha")          # regression passengers (read-off, not matched)
# All three matching variables are scored. DC and DN are ordinal codes, scored as exact classes; the
# value-tier band is left NA for them.
CLS   <- c("fld_V_VG","fld_V_DC","fld_V_DN","domsp_full","domsp_raw")   # classification scored: the full matching axis + dom species
TARGETS <- c(REG, CLS)                                            # scored targets
TOL <- c(coniferfrac_full = 0.1, shannon_full = 0.25)             # value-tier band (NA -> skip)
# domsp_raw is the matched plot's own dominant species, unpooled: the plot-index product decodes to
# these raw species, so that row scores the product as delivered. domsp_full pools the rare species
# into 기타 exactly as 07a does, and stays the like-for-like comparison with direct prediction.
DT[, domsp_raw := as.character(domsp_full)]
DT[, domsp_full := as.character(domsp_full)]; rare <- names(which(table(DT$domsp_full) < MIN_SP))
DT[domsp_full %in% rare, domsp_full := "기타"]; DT[, domsp_full := factor(domsp_full)]

run_banner("m4c_dualtrack", PILOT, cfg, n = nrow(DT),
           extra = sprintf("route=impute(yai k=%d) predsets=%d targets=%d", KNN, length(TIERS), length(TARGETS)))

# ---- 3. one imputation per (predset x scheme x fold); read all targets off the matched plot ----
# helper: fit yai on `ref` rows, find k=1 neighbour for `prd` rows, return matched ref Plot_ID per prd plot.
match_neighbour <- function(ref, prd, preds, seed) {
  set.seed(seed)
  xref <- as.data.frame(DT[ref, preds, with = FALSE]); rownames(xref) <- DT$Plot_ID[ref]
  yref <- as.data.frame(DT[ref, MATCH, with = FALSE]); rownames(yref) <- DT$Plot_ID[ref]   # match on FIELD structure
  yref$fld_V_VG <- droplevels(factor(yref$fld_V_VG))
  xprd <- as.data.frame(DT[prd, preds, with = FALSE]); rownames(xprd) <- DT$Plot_ID[prd]
  ya <- yai(x = xref, y = yref, method = "randomForest", k = KNN, ntree = NTREE)
  nt <- newtargets(ya, newdata = xprd)
  data.frame(plot_id = rownames(nt$neiIdsTrgs),
             ref_id  = as.character(nt$neiIdsTrgs[, 1]), stringsAsFactors = FALSE)
}

tm$stage("cluster"); cl <- register_cluster(16L)
jobs <- list()
for (tn in names(TIERS)) {
  for (sc in c("spatial","random")) for (f in 1:K) jobs[[length(jobs)+1L]] <- list(tier=tn, sc=sc, f=f)
  jobs[[length(jobs)+1L]] <- list(tier=tn, sc="test", f=NA_integer_)
}
tm$stage("impute")
nn <- foreach(j = jobs, .combine = rbind, .packages = c("yaImpute","randomForest","data.table"),
              .export = c("DT","TIERS","MATCH","TARGETS","CLS","KNN","NTREE","SEED","match_neighbour")) %dopar% {
  preds <- TIERS[[j$tier]]
  mask  <- complete.cases(DT[, c(preds, MATCH), with = FALSE])      # refs need complete predictors + matching axis (field FT/DC/DN)
  tr <- which(DT$role == "train" & mask); te <- which(DT$role == "test" & mask)
  if (j$sc == "test") { ref <- tr; prd <- te } else {
    fo <- if (j$sc == "spatial") DT$fold_spatial[tr] else DT$fold_random[tr]
    ref <- tr[fo != j$f]; prd <- tr[fo == j$f] }
  out <- tryCatch(match_neighbour(ref, prd, preds, SEED + (if (is.na(j$f)) 0L else j$f)),
                  error = function(e) { message("job fail ", j$tier, "/", j$sc, "/", j$f, ": ", conditionMessage(e)); NULL })
  if (is.null(out)) return(NULL)
  cbind(tier = j$tier, scheme = j$sc, out)
}
stopCluster(cl)
# A failed job returns NULL and rbind drops it, which would shrink n silently: every held-out plot must
# be matched exactly once per predictor set and scheme.
for (tn in names(TIERS)) {
  m_ok <- complete.cases(DT[, c(TIERS[[tn]], MATCH), with = FALSE])
  for (sc in c("spatial", "random", "test")) {
    want <- sum(m_ok & DT$role == (if (sc == "test") "test" else "train"))
    got  <- sum(nn$tier == tn & nn$scheme == sc)
    if (got != want) stop(sprintf("07b: %s/%s matched %d of %d held-out plots - a job failed", tn, sc, got, want), call. = FALSE)
  }
}

# ---- 4. score each (predset x scheme x target): obs vs the matched plot's value ----------------
tm$stage("score")
idx <- match(nn$plot_id, DT$Plot_ID); ridx <- match(nn$ref_id, DT$Plot_ID)
rows <- list()
# clu = the NFI cluster of each held-out plot: the bootstrap resamples clusters, not plots, because the
# subplots of one cluster sit 50 m apart and their residuals are correlated (R/eval_utils.R bci()).
metric_row <- function(target, predset, scheme, o, p, cell, is_class, tol, stratum = "full", clu = NULL) {
  if (is_class) { pt <- mean(p == o); vl <- NA_real_; wn <- NA_real_; nr <- nir(o)   # k=1 -> no top-3
                  ci <- bci(o, p, function(a,b) mean(a == b), NBOOT, cluster = clu); rm_ <- NA_real_ }
  else { o <- as.numeric(o); p <- as.numeric(p)
         pt <- r2(o, p); vl <- if (is.na(tol)) NA_real_ else mean(abs(o - p) <= tol)
         ao <- tapply(o, cell, mean); ap <- tapply(p, cell, mean)[names(ao)]
         wn <- if (length(ao) > 2) r2(as.numeric(ao), as.numeric(ap)) else NA_real_
         nr <- NA_real_; ci <- bci(o, p, r2, NBOOT, cluster = clu); rm_ <- rmse(o, p) }
  rows[[length(rows)+1L]] <<- data.frame(target = target, stratum = stratum, predset = predset, scheme = scheme,
    route = "impute", point = round(pt,4), value = round(vl,4), window = round(wn,4),
    rmse = if (is.na(rm_)) NA else round(rm_,4), nir = round(nr,4),
    ci_lo = round(ci[1],4), ci_hi = round(ci[2],4), n = length(o), testbed = VERSION)
}
# Crown density at the one class boundary both systems share: the field records 0-50/51-75/76-100 %,
# the map 0-50/51-70/71-100 %, so only <= 50 against > 50 compares like with like. Scored for the
# imputation and, in the baseline block, for the map.
bin50 <- function(x) ifelse(as.numeric(x) <= 50, "le50", "gt50")
for (tn in names(TIERS)) for (sc in c("spatial","random","test")) {
  sel <- which(nn$tier == tn & nn$scheme == sc); if (!length(sel)) next
  oi <- idx[sel]; ri <- ridx[sel]; cell <- DT$block_id[oi]; clu_all <- cluster_of(DT$Plot_ID[oi])
  for (tg in TARGETS) {
    is_cls <- tg %in% CLS
    tol <- if (tg %in% names(TOL)) TOL[[tg]] else NA_real_
    o <- if (is_cls) as.character(DT[[tg]][oi]) else DT[[tg]][oi]      # held-out truth
    p <- if (is_cls) as.character(DT[[tg]][ri]) else DT[[tg]][ri]      # matched plot's value
    keep <- !is.na(o) & !is.na(p)
    oo <- o[keep]; pp <- p[keep]; cc <- cell[keep]; ck <- clu_all[keep]
    metric_row(tg, tn, sc, oo, pp, cc, is_cls, tol, stratum = "full", clu = ck)
    # FROR-stratified imputation reliability (eval-stratum = held-out plot's origin; mirrors 07f pooled eval).
    # spatial + test only (random-CV stratum adds little); needs >2 plots per stratum for a stable metric.
    if (has_fror && sc %in% c("spatial","test")) {
      lab <- DT$FROR_lab[oi][keep]
      for (st in c("natural","planted")) { ii <- which(lab == st)
        if (length(ii) > 2) metric_row(tg, tn, sc, oo[ii], pp[ii], cc[ii], is_cls, tol, stratum = st, clu = ck[ii]) }
    }
  }
  o2 <- bin50(DT$fld_V_DN[oi]); p2 <- bin50(DT$fld_V_DN[ri]); k2 <- !is.na(o2) & !is.na(p2)
  if (sum(k2) > 2) metric_row("fld_V_DN2", tn, sc, o2[k2], p2[k2], cell[k2], TRUE, NA_real_, stratum = "full", clu = clu_all[k2])
  # same-axis forest-type-map baseline (imputation-INDEPENDENT): the field class (truth) against the map's
  # class at the held-out plot, for each of the three matching variables and for crown density at the 50 %
  # cut. Answers "does the plot-index imputation beat simply reading the class off the map?" Emitted once
  # (set-independent) under predset="FTM"; compare its point accuracy to the imputed field row of the same
  # target. The map is a second observation of the same place, and the full set receives the map's
  # classes as predictors.
  if (tn == "full" && "V_VG" %in% names(DT) && sc %in% c("spatial","test")) {
    for (v in c("VG", "DC", "DN")) {
      fv <- paste0("fld_V_", v); rv <- paste0("V_", v)
      if (!all(c(fv, rv) %in% names(DT))) next
      # V_DC/V_DN were coerced to numeric for the matching axis, so read them back as character here,
      # exactly as the fld_V_DC/fld_V_DN targets are scored.
      ov <- as.character(DT[[fv]][oi]); pv <- as.character(DT[[rv]][oi]); kk <- !is.na(ov) & !is.na(pv)
      if (sum(kk) > 2) metric_row(sprintf("fld_V_%s_ftmbaseline", v), "FTM", sc, ov[kk], pv[kk],
                                  cell[kk], TRUE, NA_real_, stratum = "full", clu = clu_all[kk])
    }
    if ("V_DN" %in% names(DT)) {
      o2 <- bin50(DT$fld_V_DN[oi]); p2 <- bin50(DT$V_DN[oi]); k2 <- !is.na(o2) & !is.na(p2)
      if (sum(k2) > 2) metric_row("fld_V_DN2_ftmbaseline", "FTM", sc, o2[k2], p2[k2], cell[k2], TRUE, NA_real_, stratum = "full", clu = clu_all[k2])
    }
  }
}

# ---- 4b. per-plot predictions (spatial-CV) for the accuracy figures (mirror m7_uncertainty_points) ----
# For the predicted-against-observed figure: each held-out plot's observed value and the
# value read off its k=1 matched TRAIN plot, per predictor set x target. obs/pred kept as CHARACTER so the
# categorical targets ride along with the 3 continuous ones; the figure coerces the continuous ones to numeric.
pp <- list()
for (tn in names(TIERS)) {
  sel <- which(nn$tier == tn & nn$scheme == "spatial"); if (!length(sel)) next
  oi <- idx[sel]; ri <- ridx[sel]
  for (tg in TARGETS)
    pp[[length(pp) + 1L]] <- data.table(Plot_ID = DT$Plot_ID[oi], predset = tn, target = tg,
                                        obs = as.character(DT[[tg]][oi]), pred = as.character(DT[[tg]][ri]))
}
pp <- rbindlist(pp)[!is.na(obs) & !is.na(pred)]
fwrite(pp, file.path(outd, sprintf("m4c_pred_points%s.csv", SUF)), bom = TRUE)
cat(sprintf("wrote m4c_pred_points%s.csv (%d plot-predictions; route=impute, spatial-CV, per predset x target)\n", SUF, nrow(pp)))

# ---- 5. write + headline (impute vs direct, spatial folds) --------------------------------
out <- rbindlist(rows); fwrite(out, file.path(outd, sprintf("m4c_dualtrack%s.csv", SUF)))
cat("\n===== plot-index imputation (route=impute), spatial folds =====\n")
print(out[scheme == "spatial" & stratum == "full", .(target, predset, point, ci_lo, ci_hi, n)][order(target, predset)], row.names = FALSE)
if (out[stratum %in% c("natural","planted"), .N] > 0) {
  cat("\n----- imputation route by stand origin (FROR; spatial-CV, eval-stratum) -----\n")
  print(out[scheme == "spatial" & stratum %in% c("natural","planted"),
            .(target, predset, stratum, point, n)][order(target, predset, stratum)], row.names = FALSE)
}
# 07a writes one table per run: *_ftmonly for the delivered set and *_baseonly for the sets without
# the forest-type map. The six-set ladder, m4a_ladder.csv, scores both on a common sample and is
# read only when neither of those exists.
lad <- list(c("m4a_ladder_ftmonly.csv", "full"), c("m4a_ladder_baseonly.csv", "clim_s2med"))
d <- rbindlist(lapply(lad, function(x) {
  f <- file.path(outd, x[1]); if (!file.exists(f)) return(NULL)
  fread(f)[scheme == "spatial" & stratum %in% c("full","-") & predset == x[2] & target %in% TARGETS,
           .(target, predset, direct = point)] }))
if (!nrow(d) && file.exists(file.path(outd, "m4a_ladder.csv")))
  d <- fread(file.path(outd, "m4a_ladder.csv"))[scheme == "spatial" & stratum %in% c("full","-") &
         predset %in% c("clim_s2med","full") & target %in% TARGETS, .(target, predset, direct = point)]
if (nrow(d)) {
  i <- out[scheme == "spatial" & stratum == "full", .(target, predset, impute = point)]
  cmp <- merge(d, i, by = c("target","predset"))
  cat("\n----- direct (07a) vs impute (07b), spatial-CV point -----\n"); print(cmp[order(target,predset)], row.names = FALSE)
} else cat("\n(no 07a table found — run 07a first for the side-by-side)\n")
tm$done(n = nrow(DT), params = sprintf("K=%d ntree=%d nboot=%d k=%d predsets=%d targets=%d",
        K, NTREE, NBOOT, KNN, length(TIERS), length(TARGETS)))
