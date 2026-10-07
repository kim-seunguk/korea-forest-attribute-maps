# 07f_stratified_accuracy.R — accuracy split by stand origin, natural against planted.
# ============================================================================================
# On the frozen folds, splits the plots by the forest-type map's stand-origin field, planted against
# natural, and asks three things of the delivered layers:
#   (a) range restriction : is Var(target) wider in natural than planted? (the ceiling driver)
#   (b) POOLED reliability: does the ONE production model predict natural plots better than planted?
#                           (= where the shipped product is more reliable — a data-product value layer)
#   (c) WITHIN-stratum     : does a stratum-specific model do better than the pooled one?
# Hypothesis: natural forests (wider composition/diversity range) map better; plantations compress R².
#
# Engine reused from 07a: register_cluster, a job list per target through foreach, eval_utils metrics.
# Stand-origin codes: 1 = 인공림 planted, 2 = 천연림 natural, 0 = neither (non-forest).
# The script prints the distribution.
#
# Inputs : outputs/nfi8/{cv_folds.csv, m2_targets.csv} · data/interim/nfi8/rd_points.csv · data/interim/nfi8/plot_features_s2.csv
# Output : outputs/nfi8/m6_stratified<SUF>.csv (target,predset,model,fit_stratum,eval_stratum,scheme,point,n,var)
#          SUF names the run: "" for both sets on their common sample, "_baseonly" and "_ftm" for the
#          switches below, with "_PILOT" in front for a pilot.
#   FTM_ONLY <- TRUE; source(here::here("pipeline", "07f_stratified_accuracy.R"))   # the delivered set
#   source(...) alone runs both sets; PILOT <- TRUE first for a quick shakedown.
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table); library(randomForest); library(foreach); library(doParallel) })
source(here::here("R","paths.R")); source(here::here("R","run_utils.R")); source(here::here("R","eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
if (!exists("PILOT"))   PILOT   <- FALSE
SEED <- 7214L
cfg  <- run_sizing(VERSION, "m6_stratified", pilot = PILOT, budget_min = if (PILOT) 25 else 150)
K <- cfg$K; NTREE <- cfg$ntree
TARGETS <- c("coniferfrac_full", "shannon_full", "Bio_a_ha")        # headline composition + ceiling + structure
STRATA  <- c("planted", "natural")                                  # within-stratum + per-stratum eval (drop 'other')
outd <- here("outputs", VERSION); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
SUF  <- if (PILOT) "_PILOT" else ""                                 # a pilot writes its own file

tm <- run_timer("m6_stratified", version = VERSION, budget_min = if (PILOT) 25 else 150)
# ---- 1. load + merge (same as 07a) ---------------------------------------------------------
tm$stage("load")
folds <- fread(file.path(outd,"cv_folds.csv"), colClasses=list(character="Plot_ID"))
m2    <- fread(file.path(outd,"m2_targets.csv"), colClasses=list(character="Plot_ID"))
rd    <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding="UTF-8", colClasses=list(character="Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses=list(character="Plot_ID"))
dupc  <- intersect(setdiff(names(m2),"Plot_ID"), names(rd)); if (length(dupc)) rd <- rd[, setdiff(names(rd),dupc), with=FALSE]
DT <- Reduce(function(a,b) merge(a,b,by="Plot_ID",all.x=TRUE),
             list(folds[,.(Plot_ID,role,fold_spatial,block_id)], rd, m2, feat))

# ---- 2. predictor sets + FROR strata -------------------------------------------------------
preds_clim <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
s2_med <- intersect(names(DT), S2_BAND_NAMES)                      # the 42 declared bands (R/paths.R)
vstar  <- intersect(V_PREDICTORS, names(DT))   # stand age is not among them; matches 07a, 07b and 07d
DT[, V_VG := factor(as.character(V_VG))]; for(v in c("V_DC","V_DN")) DT[[v]] <- as.numeric(DT[[v]])
PSETS <- list(clim_s2med = c(preds_clim, s2_med),                   # no-FTM production set
              full       = unique(c(preds_clim, s2_med, vstar)))    # +FTM (FROR is the stratifier here, not a predictor)
cat("\n-- FROR distribution (1=planted / 2=natural / 0=non-forest) --\n"); print(DT[, .N, by=FROR][order(FROR)])
FROR_lab <- c("0"="other","1"="planted","2"="natural")[as.character(DT$FROR)]   # codebook-confirmed labels
DT[, FROR_lab := FROR_lab]
cat(sprintf("labelled: planted=%d natural=%d other/NA=%d\n",
            sum(FROR_lab=="planted",na.rm=TRUE), sum(FROR_lab=="natural",na.rm=TRUE), sum(is.na(FROR_lab)|FROR_lab=="other")))
fold_sp <- DT$fold_spatial
# The completeness mask below is built from the union of the sets in play. With both sets it is the
# common subset of plots complete under both; each switch keeps one set and so widens the sample to
# every plot complete under it, and writes its own file.
#   BASE_ONLY  the no-map set alone                                                  -> *_baseonly
#   FTM_ONLY   the delivered set alone. The pooled model must be the product's model, and the
#              product's set includes stand origin (05a), so FROR is added as a predictor as well as
#              being the stratifier. Within a stratum it is constant, so the stratum-specific models
#              lose nothing.                                                          -> *_ftm
if (exists("BASE_ONLY") && isTRUE(BASE_ONLY)) { PSETS <- PSETS["clim_s2med"]; SUF <- paste0(SUF, "_baseonly")
                                                 cat("BASE_ONLY: predsets = clim_s2med\n") }
if (exists("FTM_ONLY") && isTRUE(FTM_ONLY)) {
  DT[, FROR := factor(as.character(FROR))]
  PSETS <- list(full = unique(c(PSETS$full, "FROR"))); SUF <- paste0(SUF, "_ftm"); cat("FTM_ONLY: predsets = full (+FROR)\n")
}
all_pred <- unique(unlist(PSETS))
run_banner("m6_stratified", PILOT, cfg, n=nrow(DT), extra=sprintf("targets=%d strata=%s", length(TARGETS), paste(STRATA,collapse="/")))

# ---- 3. cluster + per-target job list (pooled + within) ------------------------------------
tm$stage("cluster"); cl <- register_cluster(16L)
rows <- list()
add <- function(...) rows[[length(rows)+1L]] <<- data.frame(..., stringsAsFactors=FALSE)
run_target <- function(tg) {
  mask <- complete.cases(DT[, c(all_pred, tg), with=FALSE]) & !is.na(DT$FROR_lab) & DT$FROR_lab %in% c(STRATA)
  tr <- which(DT$role=="train" & mask); te <- which(DT$role=="test" & mask)
  jobs <- list()
  for (ps in names(PSETS)) {
    for (f in 1:K) jobs[[length(jobs)+1L]] <- list(ps=ps, mode="pooled", strat=NA_character_, f=f)
    jobs[[length(jobs)+1L]] <- list(ps=ps, mode="pooled", strat=NA_character_, f=NA_integer_)        # pooled test
    for (s in STRATA) { for (f in 1:K) jobs[[length(jobs)+1L]] <- list(ps=ps, mode="within", strat=s, f=f)
                        jobs[[length(jobs)+1L]] <- list(ps=ps, mode="within", strat=s, f=NA_integer_) } # within test
  }
  res <- foreach(j=jobs, .combine=rbind, .packages=c("randomForest","data.table"),
                 .export=c("DT","NTREE","SEED","PSETS","FROR_lab","fold_sp")) %dopar% {
    preds <- PSETS[[j$ps]]
    base_tr <- if (j$mode=="pooled") tr else tr[FROR_lab[tr]==j$strat]
    if (is.na(j$f)) { trr <- base_tr; pidx <- if (j$mode=="pooled") te else te[FROR_lab[te]==j$strat]; sc <- "test" }
    else { trr <- base_tr[fold_sp[base_tr]!=j$f]; pidx <- base_tr[fold_sp[base_tr]==j$f]; sc <- "spatial" }
    if (!length(pidx) || !length(trr)) return(NULL)
    set.seed(SEED + (if (is.na(j$f)) 0L else j$f))
    rf <- randomForest(x=as.data.frame(DT[trr, preds, with=FALSE]), y=DT[[tg]][trr], ntree=NTREE)
    data.frame(ps=j$ps, mode=j$mode, strat=ifelse(is.na(j$strat),"-",j$strat), scheme=sc,
               pos=pidx, pred=as.numeric(predict(rf, as.data.frame(DT[pidx, preds, with=FALSE]))), stringsAsFactors=FALSE)
  }
  # ---- assemble ----
  for (ps in names(PSETS)) for (sc in c("spatial","test")) {
    # pooled: one model, evaluated overall + per eval-stratum
    dp <- res[res$ps==ps & res$mode=="pooled" & res$scheme==sc, ]
    if (nrow(dp)) { o <- DT[[tg]][dp$pos]; lab <- FROR_lab[dp$pos]
      add(target=tg, predset=ps, model="pooled", fit_stratum="all", eval_stratum="all", scheme=sc,
          point=round(r2(o, dp$pred),4), n=length(o), var=round(var(o),4))
      for (s in STRATA) { i <- lab==s
        if (sum(i)>2) add(target=tg, predset=ps, model="pooled", fit_stratum="all", eval_stratum=s, scheme=sc,
                          point=round(r2(o[i], dp$pred[i]),4), n=sum(i), var=round(var(o[i]),4)) } }
    # within: stratum-specific model, evaluated on its own stratum
    for (s in STRATA) { dw <- res[res$ps==ps & res$mode=="within" & res$strat==s & res$scheme==sc, ]
      if (nrow(dw)>2) { o <- DT[[tg]][dw$pos]
        add(target=tg, predset=ps, model="within", fit_stratum=s, eval_stratum=s, scheme=sc,
            point=round(r2(o, dw$pred),4), n=length(o), var=round(var(o),4)) } }
  }
}
for (i in seq_along(TARGETS)) { tm$stage(sprintf("t:%s",TARGETS[i])); progress(i,length(TARGETS),TARGETS[i]); run_target(TARGETS[i]) }
stopCluster(cl)

# ---- 4. write + headline -------------------------------------------------------------------
out <- rbindlist(rows); fwrite(out, file.path(outd, sprintf("m6_stratified%s.csv", SUF)))
cat(sprintf("\n===== stand origin: planted vs natural (spatial folds, predset %s) =====\n", names(PSETS)[1]))
print(out[scheme=="spatial" & predset==names(PSETS)[1], .(target,model,fit_stratum,eval_stratum,point,n,var)], row.names=FALSE)
cat("\nRange restriction: compare `var` planted vs natural (per target).\n")
cat("Reliability: pooled model eval_stratum=natural vs planted `point`.\n")
cat("Stratum model: model=within vs model=pooled(eval=same stratum) `point`.\n")
tm$done(n=nrow(DT), params=sprintf("K=%d ntree=%d targets=%d", K, NTREE, length(TARGETS)))
