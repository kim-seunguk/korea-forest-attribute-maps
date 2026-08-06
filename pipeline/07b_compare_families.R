# 07b_compare_families.R — score the plot-index imputation on the same folds as direct prediction.
# ================================================================================================
# The dataset ships two product families, and they are only comparable if both are scored the same
# way. Stage 07a scores the direct random forests; this scores the imputation, reusing the
# frozen spatial folds and the same disjoint test set, and writing rows in the same schema.
#
# Engine: one yaImpute model per tier, scheme and fold. yai(method = "randomForest", k = 1) matches
# each held-out plot to its single most similar training plot in random-forest proximity space.
#
# The matching axis is the plot's field forest structure — forest type, diameter class and density
# class. Those are strong, directly observed classes. The four published attributes are not matched
# on; they are read off whichever plot the structure selected, so the accuracy reported here is what
# a user gets rather than what the matching was tuned to produce, and it mirrors how the delivered
# product is consumed.
#
# Rows are also emitted split by stand origin, natural against planted, mirroring the same split on
# the direct route.
#
# Tiers  : base = Sentinel-2 with climate and terrain · FTM-added = base plus the forest-type map
#          attributes and stand origin.
# Scored : forest type (top-1, the same axis the matching used) · dominant species (top-1) ·
#          coniferous fraction, Shannon diversity and biomass (R2), plus a forest-type-map baseline.
#
# In : outputs/nfi8/{cv_folds.csv, m2_targets.csv} · data/interim/nfi8/rd_points.csv
#      outputs/nfi8/m4a_ladder.csv, if 07a has run, for the direct-against-imputation table
#      data/interim/nfi8/plot_features_s2.csv
# Out: outputs/nfi8/m4c_dualtrack.csv · outputs/nfi8/m4c_pred_points.csv
#
# Heavy. Runs the full pass, which is what the reported numbers come from. For a quick shakedown
# on fewer trees and fewer bootstrap resamples, set PILOT first:
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
s2_med <- sort(intersect(names(DT), S2_BAND_NAMES))   # sort() = REPRODUCIBILITY GUARD, not cosmetic: the yaImpute
# k=1 randomForest match samples mtry columns in a seed-driven order, so the SAME predictors in a
# different COLUMN ORDER fit a different forest and change a large share of the matched plots.
# s2_med's order otherwise inherits the input CSV's column order, which is regenerated per machine
# and not guaranteed identical between machines, so the scores would drift with it. Must match
# pipeline/05b_plot_index_imputation.R, or the validation and the deposited product would order predictors
# differently.
vstar  <- intersect(V_PREDICTORS, names(DT))            # map predictors; stand age is not among them
has_fror <- "FROR" %in% names(DT)
DT[, V_VG := factor(as.character(V_VG))]; for (v in c("V_DC","V_DN")) if (v %in% names(DT)) DT[[v]] <- as.numeric(DT[[v]])
# matching axis = FIELD structure {FT, DC, DN} (Riley-faithful; identical typing to pipeline/05b_plot_index_imputation.R)
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
TIERS <- list(clim_s2med = c(preds_clim, s2_med),                   # S2
              full       = unique(c(preds_clim, s2_med, vstar, if (has_fror) "FROR")))   # S2+FTM
cat(sprintf("tiers: clim_s2med=%d vars | full=%d vars (s2_med=%d vstar=%d fror=%s)\n",
            length(TIERS$clim_s2med), length(TIERS$full), length(s2_med), length(vstar), has_fror))

# MATCH = yaImpute matching axis (field structure); TARGETS = the SCORED set read off the matched plot.
MATCH <- c("fld_V_VG","fld_V_DC","fld_V_DN")                       # field FT + diameter/density classes (== 05b)
REG   <- c("coniferfrac_full","shannon_full","Bio_a_ha")          # regression passengers (read-off, not matched)
# All three matching variables are scored, not forest type alone. The matching runs on {FT, DC, DN},
# so scoring only FT would leave two thirds of the axis unvalidated, and diameter class is what
# carries stand size to a user. DC and DN are ordinal codes, scored as exact classes, which is the
# strict reading; the value-tier band is left NA for them.
CLS   <- c("fld_V_VG","fld_V_DC","fld_V_DN","domsp_full")         # classification scored: the full matching axis + dom species
TARGETS <- c(REG, CLS)                                            # scored targets
TOL <- c(coniferfrac_full = 0.1, shannon_full = 0.25)             # value-tier band (NA -> skip)
# rare dominant species -> 기타 (matches 07a so the legend is comparable)
DT[, domsp_full := as.character(domsp_full)]; rare <- names(which(table(DT$domsp_full) < MIN_SP))
DT[domsp_full %in% rare, domsp_full := "기타"]; DT[, domsp_full := factor(domsp_full)]

run_banner("m4c_dualtrack", PILOT, cfg, n = nrow(DT),
           extra = sprintf("route=impute(yai k=%d) tiers=%d targets=%d", KNN, length(TIERS), length(TARGETS)))

# ---- 3. one imputation per (tier x scheme x fold); read all targets off the matched plot ----
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

# ---- 4. score each (tier x scheme x target): obs vs the matched plot's value ----------------
tm$stage("score")
idx <- match(nn$plot_id, DT$Plot_ID); ridx <- match(nn$ref_id, DT$Plot_ID)
rows <- list()
metric_row <- function(target, predset, scheme, o, p, cell, is_class, tol, stratum = "full") {
  if (is_class) { pt <- mean(p == o); vl <- NA_real_; wn <- NA_real_; nr <- nir(o)   # k=1 -> no top-3
                  ci <- bci(o, p, function(a,b) mean(a == b), NBOOT); rm_ <- NA_real_ }
  else { o <- as.numeric(o); p <- as.numeric(p)
         pt <- r2(o, p); vl <- if (is.na(tol)) NA_real_ else mean(abs(o - p) <= tol)
         ao <- tapply(o, cell, mean); ap <- tapply(p, cell, mean)[names(ao)]
         wn <- if (length(ao) > 2) r2(as.numeric(ao), as.numeric(ap)) else NA_real_
         nr <- NA_real_; ci <- bci(o, p, r2, NBOOT); rm_ <- rmse(o, p) }
  rows[[length(rows)+1L]] <<- data.frame(target = target, stratum = stratum, predset = predset, scheme = scheme,
    route = "impute", point = round(pt,4), value = round(vl,4), window = round(wn,4),
    rmse = if (is.na(rm_)) NA else round(rm_,4), nir = round(nr,4),
    ci_lo = round(ci[1],4), ci_hi = round(ci[2],4), n = length(o), testbed = VERSION)
}
for (tn in names(TIERS)) for (sc in c("spatial","random","test")) {
  sel <- which(nn$tier == tn & nn$scheme == sc); if (!length(sel)) next
  oi <- idx[sel]; ri <- ridx[sel]; cell <- DT$block_id[oi]
  for (tg in TARGETS) {
    is_cls <- tg %in% CLS
    tol <- if (tg %in% names(TOL)) TOL[[tg]] else NA_real_
    o <- if (is_cls) as.character(DT[[tg]][oi]) else DT[[tg]][oi]      # held-out truth
    p <- if (is_cls) as.character(DT[[tg]][ri]) else DT[[tg]][ri]      # matched plot's value
    keep <- !is.na(o) & !is.na(p)
    oo <- o[keep]; pp <- p[keep]; cc <- cell[keep]
    metric_row(tg, tn, sc, oo, pp, cc, is_cls, tol, stratum = "full")
    # FROR-stratified imputation reliability (eval-stratum = held-out plot's origin; mirrors 07f pooled eval).
    # spatial + test only (random-CV stratum adds little); needs >2 plots per stratum for a stable metric.
    if (has_fror && sc %in% c("spatial","test")) {
      lab <- DT$FROR_lab[oi][keep]
      for (st in c("natural","planted")) { ii <- which(lab == st)
        if (length(ii) > 2) metric_row(tg, tn, sc, oo[ii], pp[ii], cc[ii], is_cls, tol, stratum = st) }
    }
  }
  # same-axis FTM baseline (imputation-INDEPENDENT): field FT (truth) vs FTM raster-X FT at the held-out plot.
  # Answers "does the plot-index imputation beat simply reading the forest type off the map?"
  # Emitted once (tier-independent) under predset="FTM"; compare its point acc to the imputed field-FT row.
  if (tn == "full" && "V_VG" %in% names(DT) && sc %in% c("spatial","test")) {
    of <- as.character(DT$fld_V_VG[oi]); pf <- as.character(DT$V_VG[oi]); kk <- !is.na(of) & !is.na(pf)
    if (sum(kk) > 2) metric_row("fld_V_VG_ftmbaseline", "FTM", sc, of[kk], pf[kk], cell[kk], TRUE, NA_real_, stratum = "full")
  }
}

# ---- 4b. per-plot predictions (spatial-CV) for the accuracy figures (mirror m7_uncertainty_points) ----
# For the predicted-against-observed figure: each held-out plot's observed value and the
# value read off its k=1 matched TRAIN plot, per tier x target. obs/pred kept as CHARACTER so the categorical
# domsp rides along with the 3 continuous targets; the figure coerces the continuous ones to numeric.
pp <- list()
for (tn in names(TIERS)) {
  sel <- which(nn$tier == tn & nn$scheme == "spatial"); if (!length(sel)) next
  oi <- idx[sel]; ri <- ridx[sel]
  for (tg in TARGETS)
    pp[[length(pp) + 1L]] <- data.table(Plot_ID = DT$Plot_ID[oi], predset = tn, target = tg,
                                        obs = as.character(DT[[tg]][oi]), pred = as.character(DT[[tg]][ri]))
}
pp <- rbindlist(pp)[!is.na(obs) & !is.na(pred)]
fwrite(pp, file.path(outd, "m4c_pred_points.csv"), bom = TRUE)
cat(sprintf("wrote m4c_pred_points.csv (%d plot-predictions; route=impute, spatial-CV, per tier x target)\n", nrow(pp)))

# ---- 5. write + headline (impute vs direct, spatial folds) --------------------------------
out <- rbindlist(rows); fwrite(out, file.path(outd, "m4c_dualtrack.csv"))
cat("\n===== plot-index imputation (route=impute), spatial folds =====\n")
print(out[scheme == "spatial" & stratum == "full", .(target, predset, point, ci_lo, ci_hi, n)][order(target, predset)], row.names = FALSE)
if (out[stratum %in% c("natural","planted"), .N] > 0) {
  cat("\n----- imputation route by stand origin (FROR; spatial-CV, eval-stratum) -----\n")
  print(out[scheme == "spatial" & stratum %in% c("natural","planted"),
            .(target, predset, stratum, point, n)][order(target, predset, stratum)], row.names = FALSE)
}
m4a_f <- file.path(outd, "m4a_ladder.csv")
if (file.exists(m4a_f)) {
  d <- fread(m4a_f)[scheme == "spatial" & stratum %in% c("full","-") & predset %in% c("clim_s2med","full") &
                      target %in% TARGETS, .(target, predset, direct = point)]
  i <- out[scheme == "spatial" & stratum == "full", .(target, predset, impute = point)]
  cmp <- merge(d, i, by = c("target","predset"))
  cat("\n----- direct (07a) vs impute (07b), spatial-CV point -----\n"); print(cmp[order(target,predset)], row.names = FALSE)
} else cat("\n(m4a_ladder.csv not found — run 07a first for the side-by-side)\n")
tm$done(n = nrow(DT), params = sprintf("K=%d ntree=%d nboot=%d k=%d tiers=%d targets=%d",
        K, NTREE, NBOOT, KNN, length(TIERS), length(TARGETS)))
