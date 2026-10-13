# 07e_forest_type_map_tiers.R — how much the forest-type map adds, step by step.
# ============================================================================================
# Stage 07a compares the predictor sets on one common sample, the plots complete under all of them,
# which is restricted to plots the forest-type map covers. This re-evaluates each step on its own
# complete sample instead, so that a configuration needing no national map is credited with the
# plots it can actually use:
#   none    = climate and terrain + Sentinel-2 medians    (no forest-type map at all)
#   limited = + stand origin                              (the one map field that is reliable)
#   maximal = + the forest-type map stand attributes      (forest type, diameter class, density class)
# The three labels are local to this comparison, and none of them is the delivered predictor set,
# which adds stand origin AND the stand attributes; 07a scores that set with FTM_ONLY. What this
# measures is the portability-against-fit trade-off: R², skill and n per step, on spatial CV and on
# the test set, and how many plots each step costs.
#
# Inputs : outputs/nfi8/{cv_folds.csv,m2_targets.csv} · data/interim/nfi8/rd_points.csv · data/interim/nfi8/plot_features_s2.csv
# Output : outputs/nfi8/m5_tiers<SUF>.csv  (tier,target,scheme,point,value,nir,delta_vs_none,n); SUF = "_PILOT" for a pilot
#   source(...) runs the full pass; PILOT <- TRUE first for a quick shakedown, which writes its own file.
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table); library(randomForest); library(foreach); library(doParallel) })
source(here::here("R","paths.R")); source(here::here("R","run_utils.R")); source(here::here("R","eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
if (!exists("PILOT"))   PILOT   <- FALSE
SEED <- 7214L; MIN_SP <- 20L
cfg <- run_sizing(VERSION, "m5_tiers", pilot=PILOT, budget_min=if (PILOT) 25 else 150)
K <- cfg$K; NTREE <- cfg$ntree
outd <- here("outputs", VERSION)
SUF  <- if (PILOT) "_PILOT" else ""                                # a pilot writes its own file

tm <- run_timer("m5_tiers", version=VERSION, budget_min=if (PILOT) 25 else 150)
tm$stage("load")
folds <- fread(file.path(outd,"cv_folds.csv"), colClasses=list(character="Plot_ID"))
m2    <- fread(file.path(outd,"m2_targets.csv"), colClasses=list(character="Plot_ID"))
rd    <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding="UTF-8", colClasses=list(character="Plot_ID"))
feat  <- fread(require_file(here("data","interim",VERSION,"plot_features_s2.csv")), colClasses=list(character="Plot_ID"))
dupc  <- intersect(setdiff(names(m2),"Plot_ID"), names(rd)); if (length(dupc)) rd <- rd[, setdiff(names(rd),dupc), with=FALSE]
DT <- Reduce(function(a,b) merge(a,b,by="Plot_ID",all.x=TRUE),
             list(folds[,.(Plot_ID,role,fold_spatial,block_id)], rd, m2, feat))

preds_clim <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","sin_y","cos_x","C_PAR","Lon","Lat")
s2_med <- intersect(names(DT), S2_BAND_NAMES)                      # the 42 declared bands (R/paths.R)
if (sum(complete.cases(DT[, ..s2_med])) < 10000L)                  # guard: stale/mid-copy inputs -> abort loudly, not silent small-n
  stop("merged S2 completeness < 10k — inputs look stale/mid-update; Restart R + re-read files.", call. = FALSE)
vstar  <- intersect(V_PREDICTORS, names(DT))   # R/paths.R's three map stand attributes; no stand age
DT[, V_VG := factor(as.character(V_VG))]; for(v in c("V_DC","V_DN")) DT[[v]] <- as.numeric(DT[[v]])
if ("FROR" %in% names(DT)) DT[, FROR := factor(as.character(FROR))]
TIERS <- list(none    = c(preds_clim, s2_med),
              limited = c(preds_clim, s2_med, if ("FROR" %in% names(DT)) "FROR"),
              maximal = c(preds_clim, s2_med, vstar))
REG <- c("coniferfrac_full","shannon_full","Bio_a_ha"); CLS <- c("domsp_full")
for (clt in CLS) { DT[[clt]] <- as.character(DT[[clt]]); rare <- names(which(table(DT[[clt]]) < MIN_SP))
                   DT[DT[[clt]] %in% rare, (clt) := "기타"]; DT[[clt]] <- factor(DT[[clt]]) }
fold_sp <- DT$fold_spatial
run_banner("m5_tiers", PILOT, cfg, n=nrow(DT), extra=sprintf("tiers=%d targets=%d", length(TIERS), length(REG)+length(CLS)))

tm$stage("cluster"); cl <- register_cluster(16L)
rows <- list(); base_pt <- list()
add <- function(...) rows[[length(rows)+1L]] <<- data.frame(..., stringsAsFactors=FALSE)
run_one <- function(tg, is_class) {
  for (tname in names(TIERS)) {
    preds <- unique(TIERS[[tname]])
    mask <- complete.cases(DT[, c(preds, tg), with=FALSE])
    tr <- which(DT$role=="train" & mask); te <- which(DT$role=="test" & mask)
    # spatial OOF (parallel folds) + test
    jobs <- c(lapply(1:K, function(f) list(sc="spatial", f=f)), list(list(sc="test", f=NA_integer_)))
    res <- foreach(j=jobs, .combine=rbind, .packages=c("randomForest","data.table"),
                   .export=c("DT","NTREE","SEED","fold_sp")) %dopar% {       # globals explicit; tr/te/preds/tg/is_class auto-detected
      if (j$sc=="test") { trr <- tr; pidx <- te } else { trr <- tr[fold_sp[tr]!=j$f]; pidx <- tr[fold_sp[tr]==j$f] }
      set.seed(SEED + (if (is.na(j$f)) 0L else j$f)); y <- DT[[tg]][trr]; if (is_class) y <- droplevels(factor(y))
      rf <- randomForest(x=as.data.frame(DT[trr, preds, with=FALSE]), y=y, ntree=NTREE)
      if (is_class) { pb <- predict(rf, as.data.frame(DT[pidx, preds, with=FALSE]), type="prob")
        t3 <- apply(pb,1,function(p) names(sort(p,decreasing=TRUE))[1:3])
        data.frame(sc=j$sc, pos=pidx, pred=colnames(pb)[max.col(pb)],
                   in3=vapply(seq_along(pidx), function(i) as.character(DT[[tg]][pidx][i]) %in% t3[,i], logical(1)), stringsAsFactors=FALSE)
      } else data.frame(sc=j$sc, pos=pidx, pred=as.numeric(predict(rf, as.data.frame(DT[pidx, preds, with=FALSE]))), in3=NA) }
    for (sc in c("spatial","test")) {
      d <- res[res$sc==sc, ]; d <- d[order(d$pos), ]; o <- DT[[tg]][d$pos]
      pt <- if (is_class) mean(d$pred==o) else r2(o, d$pred)
      vl <- if (is_class) mean(d$in3) else NA_real_; nr <- if (is_class) nir(o) else NA_real_
      if (tname=="none") base_pt[[paste(tg,sc)]] <<- pt
      add(tier=tname, target=tg, scheme=sc, point=round(pt,4), value=round(vl,4), nir=round(nr,4),
          delta_vs_none=round(pt - base_pt[[paste(tg,sc)]], 4), n=length(o))
    }
  }
}
for (i in seq_along(c(REG,CLS))) { tg <- c(REG,CLS)[i]; tm$stage(sprintf("t:%s",tg)); progress(i,length(c(REG,CLS)),tg)
  run_one(tg, is_class = tg %in% CLS) }
stopCluster(cl)

out <- rbindlist(rows); fwrite(out, file.path(outd, sprintf("m5_tiers%s.csv", SUF)))
cat("\n===== forest-type map tiers (per-predset masked; spatial folds) =====\n")
print(out[scheme=="spatial", .(target,tier,point,delta_vs_none,n)], row.names=FALSE)
cat("\nPortability: tier=none n vs maximal n (how many more plots no-FTM uses).\n")
cat("Enhancement: delta_vs_none for maximal (what FTM adds), per target.\n")
tm$done(n=nrow(DT), params=sprintf("K=%d ntree=%d tiers=%d", K, NTREE, length(TIERS)))
