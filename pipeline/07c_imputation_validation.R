# 07c_imputation_validation.R — validate the delivered index raster itself, at the held-out plots.
# ================================================================================================
# Stage 07b scores the imputation as a model. This scores the raster that ships. It samples the
# delivered index raster at the spatially disjoint test plots, reads the matched plot's observed
# attributes — the values a user actually receives at that cell — and compares them against the test
# plots' own field measurements.
#
# ! Read the self-match fraction before reading the metrics. The delivered map uses every plot as a
#   donor, so a test plot is itself in the library, and the cell it sits in often matches it back to
#   itself. Where that happens the comparison is a plot against its own record and tells you
#   nothing. The metrics are therefore reported twice: over all held-out plots, and over only those
#   matched to a different plot. The second is the one that measures anything.
#
# This is supplementary. The reported accuracy of the plot-index family is the model-level figure
# from 07b, measured on held-out plots under spatial-block cross-validation. What this adds is the
# delivered raster's own behaviour, which is not the same measurement: the raster matches on
# predictors resampled to 30 m at the cell rather than on the plot's own extracted feature vector,
# and it draws on a library containing every plot rather than a training fold.
#
# For biomass and diversity the R2 is negative by construction: copying a single neighbour's value
# roughly doubles the variance of the estimate. That is the cost of returning a set of attributes
# that occurred together on one real plot.
#
#   source(here::here("pipeline", "07c_imputation_validation.R"))
# ================================================================================================
suppressPackageStartupMessages({ library(here); library(terra); library(data.table); library(sf) })
source(here::here("R","paths.R")); source(here::here("R","eval_utils.R"))   # r2(), nir()
# The delivered raster, imputation_ftm.tif, built by 05b from the predictor set 07b calls "full".
.tsuf     <- tier_suffix("full")                            # "_ftm"; see R/paths.R
.m4c_pred <- "full"
outd    <- here("outputs", VERSION)
idx_tif <- require_file(file.path(outd, "index",  sprintf("imputation%s.tif", .tsuf)))
ref_csv <- require_file(file.path(outd, "models", sprintf("train_reference%s.csv", .tsuf)))
r2m <- function(o, p) { ok <- is.finite(o) & is.finite(p); if (sum(ok) < 3) return(NA_real_); r2(o[ok], p[ok]) }

## --- observed TEST plots + targets + coords + plot_idx -> matched-plot lookup -----------------
folds <- fread(file.path(outd,"cv_folds.csv"),   colClasses=list(character="Plot_ID"))
m2    <- fread(file.path(outd,"m2_targets.csv"), encoding="UTF-8", colClasses=list(character="Plot_ID"))
rd    <- fread(require_file(here("data","interim",VERSION,"rd_points.csv")), encoding="UTF-8",
               select=c("Plot_ID","Lon","Lat"), colClasses=list(character="Plot_ID"))
ref   <- fread(ref_csv, encoding="UTF-8", colClasses=list(character="Plot_ID"))   # plot_idx (1..N) -> matched Plot_ID
setorder(ref, plot_idx); stopifnot(ref$plot_idx == seq_len(nrow(ref)))            # index raster value == row order

# domsp coding IDENTICAL to the direct-map validation: the classes the delivered map trained (05a's
# legend_used, the legend with the rare species already pooled), two regional pine races -> 소나무,
# anything outside the trained set -> 기타, so imputation & direct-map accuracies are comparable.
leg_f <- file.path(outd, "m8_maps_30m", "m8_domsp_legend_used.csv")
leg <- fread(if (file.exists(leg_f)) leg_f else here("data","species_codes.csv"), encoding="UTF-8")
SP_SET <- unique(leg$species)
code_domsp <- function(x) { x <- enc2utf8(as.character(x))
  x[x %in% c("강원지방소나무","중부지방소나무")] <- "소나무"; x[!(x %in% SP_SET)] <- "기타"; x }

obs <- Reduce(function(a,b) merge(a,b,by="Plot_ID"),
              list(folds[role=="test", .(Plot_ID)],
                   m2[, .(Plot_ID, coniferfrac=coniferfrac_full, shannon=shannon_full,
                          Bio_a_ha=Bio_a_ha, domsp=code_domsp(domsp_full))],
                   rd))
cat(sprintf("TEST plots: %d\n", nrow(obs)))

## --- extract the delivered index raster at the TEST plots -> matched TRAIN plot ---------------
im <- rast(idx_tif)
PLOT_CRS <- if (max(abs(obs$Lon), na.rm=TRUE) > 180) 5179L else 4326L             # rd_points Lon/Lat are 5179 metres despite the name
cat(sprintf("plot coords detected as EPSG:%d\n", PLOT_CRS))
pts <- vect(st_transform(st_as_sf(obs, coords=c("Lon","Lat"), crs=PLOT_CRS), crs(im)))
obs$plot_idx <- round(terra::extract(im, pts)[, 2])                               # matched-plot index at the test cell
hit <- is.finite(obs$plot_idx) & obs$plot_idx >= 1 & obs$plot_idx <= nrow(ref)
cat(sprintf("test plots on an imputed forest cell: %d / %d (off-forest/NA cells dropped)\n", sum(hit), nrow(obs)))

# the value delivered at that cell is the matched training plot's own measured target
mi <- match(ref$Plot_ID[obs$plot_idx], m2$Plot_ID)
obs$imp_coniferfrac <- m2$coniferfrac_full[mi]
obs$imp_shannon     <- m2$shannon_full[mi]
obs$imp_Bio_a_ha    <- m2$Bio_a_ha[mi]
obs$imp_domsp       <- code_domsp(m2$domsp_full[mi])

## --- self-match: is the cell returning the very plot that sits in it? -------------------------
obs[, matched_Plot_ID := NA_character_]
obs[hit, matched_Plot_ID := ref$Plot_ID[plot_idx]]
obs[, self_match := !is.na(matched_Plot_ID) & matched_Plot_ID == Plot_ID]
cat(sprintf("self-matched: %d / %d test plots on forest (%.1f%%) — excluded from the second block\n",
            sum(obs$self_match), sum(hit), 100 * sum(obs$self_match) / max(sum(hit), 1L)))

## --- metrics : R2(test) for regression, top-1 acc + NIR for domsp -----------------------------
# Scored twice, over the same code, so the two blocks are directly comparable.
score <- function(v, lab) {
  res <- rbindlist(list(
    data.table(target="coniferfrac", metric="R2_test",  value=round(r2m(v$coniferfrac, v$imp_coniferfrac),3)),
    data.table(target="shannon",     metric="R2_test",  value=round(r2m(v$shannon,     v$imp_shannon),3)),
    data.table(target="Bio_a_ha",    metric="R2_test",  value=round(r2m(v$Bio_a_ha,    v$imp_Bio_a_ha),3)),
    data.table(target="domsp",       metric="acc_test",
               value=round(mean(v$imp_domsp==v$domsp, na.rm=TRUE),3))))
  res[target=="domsp", nir := round(nir(v$domsp[!is.na(v$domsp)]),3)]
  res[, n := c(sum(is.finite(v$coniferfrac)&is.finite(v$imp_coniferfrac)),
               sum(is.finite(v$shannon)&is.finite(v$imp_shannon)),
               sum(is.finite(v$Bio_a_ha)&is.finite(v$imp_Bio_a_ha)),
               sum(!is.na(v$imp_domsp)&!is.na(v$domsp)))]
  res[, `:=`(tier="ftm", route="impute-raster", plots=lab)]
  res[]
}
res <- rbindlist(list(score(obs[hit], "all held-out"),
                      score(obs[hit & !self_match], "excluding self-matches")))

fwrite(obs, file.path(outd, sprintf("imputation_test_pred%s.csv", .tsuf)), bom=TRUE)   # per-plot
fwrite(res, file.path(outd, sprintf("imputation_map_validation%s.csv", .tsuf)))        # metrics
cat("\n===== delivered imputation, sampled at the disjoint test plots =====\n"); print(res)

## --- reconcile against the 07b plot-vector scores (test scheme, the delivered set) ------------
m4c_f <- file.path(outd, "m4c_dualtrack.csv")
if (file.exists(m4c_f)) {
  key <- c(coniferfrac_full="coniferfrac", shannon_full="shannon", Bio_a_ha="Bio_a_ha", domsp_full="domsp")
  m4c <- fread(m4c_f)[scheme=="test" & stratum=="full" & predset==.m4c_pred & target %in% names(key), .(target, m4c_plotvec=point)]
  m4c[, target := key[target]]
  # Compare against the self-match-excluded block. 07b holds its test plots out of the model, so
  # only the block that also excludes them is measuring a comparable thing.
  cmp <- merge(res[plots == "excluding self-matches", .(target, raster=value)], m4c,
               by="target", all.x=TRUE)
  cat(sprintf("\n----- delivered raster, self-matches excluded  vs  07b plot-vector (scheme=test, predset=%s) -----\n", .m4c_pred))
  print(cmp[order(match(target, c("coniferfrac","domsp","Bio_a_ha","shannon")))], row.names=FALSE)
  cat("(raster ~ 07b expected; small drop from 30 m predictor rasterisation is normal.)\n")
} else cat("\n(m4c_dualtrack.csv not found — run 07b for the side-by-side; raster metrics above still stand.)\n")

cat(sprintf("\nwrote:\n  %s/imputation_map_validation%s.csv   (metrics)\n  %s/imputation_test_pred%s.csv        (per-plot)\n",
            outd, .tsuf, outd, .tsuf))
