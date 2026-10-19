# 07s_environmental_representativeness.R — how well the inventory plots cover the country's environmental space.
# ----------------------------------------------------------------------------
# Two things follow from how well and how evenly the plots cover the climate and terrain of the
# forest they map: the extrapolation risk (prediction into environments with no nearby training
# plot), and whether a design-based estimate from the sample can differ from an area-weighted
# estimate from the map, which it will wherever the sample's environmental distribution differs from
# the area's.
#
# Compares the plots (rd_points) with the stocked-forest grid (the covariate rasters under the
# forest-type map's extent):
#   (a) per-variable distribution: sample against area quantiles, with the Kolmogorov-Smirnov distance;
#   (b) PCA environment space: the share of area cells inside the sample's convex hull (outside = extrapolation);
#   (c) representation ratio by PC1 decile: sample share / area share (> 1 over-, < 1 under-sampled);
#   (d) univariate extrapolation: area cells outside the sample's [1, 99] % range on any axis.
#
# Samples N_AREA grid cells and fits no model. Run with:
#   source(here::here("pipeline", "07s_environmental_representativeness.R"))
# Inputs : data/interim/nfi8/rd_points.csv · data/interim/covar_5179/<env>.tif · data/interim/fmap_30m/<edition>/V_VG.tif
# Outputs: outputs/nfi8/validation/env_repr_quantiles.csv · env_repr_pc1_ratio.csv   (no plot identifiers)
# ----------------------------------------------------------------------------
N_AREA <- 100000L; SEED <- 7214
ENV <- c("C_PREC","C_TMAX","C_TMIN","C_RH","C_GSL","C_FD0","C_DTR","T_ELE","T_SLO","C_PAR")

suppressPackageStartupMessages({ library(here); library(terra); library(data.table); library(sp) })
source(here::here("R", "paths.R"))

# ---- sample environment (the inventory plots, from rd_points) -------------------------------------
rd   <- fread(require_file(here::here("data", "interim", VERSION, "rd_points.csv")), encoding = "UTF-8")
samp <- rd[, ..ENV]; samp <- samp[complete.cases(samp)]
cat(sprintf("sample (inventory plots): %d plots with complete env\n", nrow(samp)))

# ---- area environment (the stocked-forest grid: the covariates under the forest-type map) -------
covdir <- here::here("data", "interim", "covar_5179")
cstk   <- terra::rast(file.path(covdir, paste0(ENV, ".tif"))); names(cstk) <- ENV
vmask  <- terra::rast(here::here("data", "interim", "fmap_30m", FTM_VINTAGE, "V_VG.tif")); names(vmask) <- "forest"
# The forest-type rasters cover a smaller extent than the covariate grid, with the same origin and
# cell size. Padding the mask with NA aligns the two without resampling.
if (!terra::compareGeom(cstk, vmask, stopOnError = FALSE)) vmask <- terra::extend(vmask, cstk)
stopifnot(terra::compareGeom(cstk, vmask, stopOnError = FALSE))
set.seed(SEED)
a <- as.data.table(terra::spatSample(c(cstk, vmask), N_AREA, method = "random", na.rm = TRUE, warn = FALSE))
area <- a[!is.na(forest), ..ENV]; area <- area[complete.cases(area)]
cat(sprintf("area (stocked-forest grid sample): %d cells with complete env\n", nrow(area)))

# ---- (a) per-variable sample vs area distribution + KS -------------------------------------------
qa <- rbindlist(lapply(ENV, function(v) data.table(
  var = v,
  samp_p05 = round(quantile(samp[[v]], .05), 2), samp_med = round(median(samp[[v]]), 2), samp_p95 = round(quantile(samp[[v]], .95), 2),
  area_p05 = round(quantile(area[[v]], .05), 2), area_med = round(median(area[[v]]), 2), area_p95 = round(quantile(area[[v]], .95), 2),
  KS_D = round(as.numeric(suppressWarnings(ks.test(samp[[v]], area[[v]])$statistic)), 3))))
cat("\n==== (a) per-variable: sample vs area (p05/median/p95) + KS distance ====\n"); print(qa)

# ---- PCA environment space (fit on area; project sample) -----------------------------------------
pca <- prcomp(area, scale. = TRUE)
ap  <- predict(pca, area)[, 1:2]; sp_ <- predict(pca, samp)[, 1:2]
varexp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
cat(sprintf("\nPCA: PC1 %.1f%% / PC2 %.1f%% of env variance\n", varexp[1], varexp[2]))

# ---- (b) sample convex-hull coverage of area points ---------------------------------------------
h    <- chull(sp_); poly <- sp_[c(h, h[1]), ]
inh  <- sp::point.in.polygon(ap[, 1], ap[, 2], poly[, 1], poly[, 2]) > 0
cat(sprintf("\n==== (b) extrapolation (PCA hull) ====\narea cells INSIDE sample env hull: %.1f%%  (outside = extrapolation: %.1f%%)\n",
            100 * mean(inh), 100 * mean(!inh)))

# ---- (c) representation ratio by PC1 decile ------------------------------------------------------
br  <- quantile(ap[, 1], seq(0, 1, .1)); br[1] <- br[1] - 1e-6; br[length(br)] <- br[length(br)] + 1e-6
abin <- cut(ap[, 1], br); sbin <- cut(sp_[, 1], br)
rr <- data.table(PC1_decile = levels(abin),
                 area_share = round(as.numeric(prop.table(table(abin))), 3),
                 samp_share = round(as.numeric(prop.table(table(factor(sbin, levels = levels(abin))))), 3))
rr[, repr_ratio := round(samp_share / area_share, 2)]
cat("\n==== (c) representation ratio by PC1 decile (samp_share / area_share) ====\n"); print(rr)

# ---- (d) univariate extrapolation: area outside sample [1,99]% on ANY axis -----------------------
qs  <- samp[, lapply(.SD, quantile, probs = c(.01, .99), na.rm = TRUE), .SDcols = ENV]
out <- rep(FALSE, nrow(area)); per <- numeric()
for (v in ENV) { lo <- qs[[v]][1]; hi <- qs[[v]][2]
  ov <- area[[v]] < lo | area[[v]] > hi; per[v] <- round(100 * mean(ov), 1); out <- out | ov }
cat(sprintf("\n==== (d) univariate extrapolation ====\narea cells outside sample [1,99]%% on ANY env axis: %.1f%%\n", 100 * mean(out)))
cat("  per-axis area% outside sample [1,99]%:\n"); print(per[order(-per)])

# ---- write + reading guide ----------------------------------------------------------------------
vdir <- here::here("outputs", VERSION, "validation"); dir.create(vdir, showWarnings = FALSE, recursive = TRUE)
fwrite(qa, file.path(vdir, "env_repr_quantiles.csv"), bom = TRUE)
fwrite(rr, file.path(vdir, "env_repr_pc1_ratio.csv"), bom = TRUE)
cat("\nReading the tables:\n",
    " - a large KS_D, or a repr_ratio far from 1, marks an environment the plots over- or under-sample\n",
    "   relative to its area, where a design-based estimate and an area-weighted map estimate can differ;\n",
    " - the share outside the hull, or outside the [1, 99] % range, is the forest area beyond the sampled\n",
    "   environmental range: the extrapolation risk of any prediction there.\n",
    sprintf("\nWrote %s/{env_repr_quantiles, env_repr_pc1_ratio}.csv\n", vdir))
