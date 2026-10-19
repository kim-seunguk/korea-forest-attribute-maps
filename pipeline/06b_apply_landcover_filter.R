# 06b_apply_landcover_filter.R — remove cells an independent land-cover map says cannot be forest.
# ================================================================================================
# The delivered maps inherit their extent from the forest-type map, which commits some non-forest to
# forest: airports, urban blocks, reservoirs, quarries, bare rock. The national land-cover map is an
# independent source, so its disagreement is informative.
#
# Removed  built, water and barren — surfaces that cannot carry trees.
# Kept     agriculture, grassland and wetland. Disagreement there often reflects real tree cover in
#          orchards and sparse or edge regeneration.
#
# Predicted biomass plays no part in the decision; the land-cover classes decide.
#
# The filter is applied in place to both product families and is idempotent: it only ever sets
# cells to NA and never restores them, so a second run changes nothing.
#
#   direct prediction      outputs/nfi8/maps_30m/*.tif and uncertainty_30m/*.tif
#   plot-index imputation  everything in outputs/nfi8/index/ — the index and the nearest-reference
#                          distance layer beside it — and attributes_ftm/spec30m_*.tif
#
# Each file is read, masked, written to a temporary file and then swapped in; the remove and the
# rename are both checked, so a file that cannot be replaced stops the run instead of leaving the
# old file in place. Stale .ovr sidecars are dropped so overviews rebuild.
#
# ! Re-run this after any step that rewrites a raster. The filter is lost on overwrite.
#   source(here::here("pipeline", "06b_apply_landcover_filter.R"))
# ================================================================================================
suppressPackageStartupMessages({ library(here); library(terra) })
source(here("R", "paths.R"))
terra::terraOptions(progress = 0)

REMOVE_CODES <- c(100L, 600L, 700L)                               # built, barren, water
CELL_HA <- 30 * 30 / 1e4
L1_TIF  <- here("data", "interim", "landcover_mcee_30m", "landcover_L1_30m.tif")
require_file(L1_TIF)

## ---- build the remove-mask (1 = remove, 0 = keep) ------------------------------------
L1    <- terra::rast(L1_TIF)
# TRUE = remove. A cell the land-cover map does not cover comes out FALSE, so nothing is removed
# without evidence. The resampled masks further down need an explicit NA -> 0, because resampling
# onto another grid does introduce NA outside this one's extent.
mask0 <- L1 %in% REMOVE_CODES

## ---- per-file apply (atomic, datatype/nodata-preserving) ---------------------------------------
apply_one <- function(f, m) {
  r <- terra::rast(f)
  # The mask and the target may sit on grids that are not cell-aligned, and cropping can only fix an
  # extent that already aligns. Nearest resampling regrids the 0/1 mask onto whatever grid the file
  # uses; a sub-cell shift does not matter for a built, water and barren removal mask.
  if (!terra::compareGeom(m, r, stopOnError = FALSE, messages = FALSE, crs = FALSE)) {
    m <- terra::subst(terra::resample(m, r, method = "near"), NA, 0L)
  }
  # Both are EPSG:5179, but imputation_ftm.tif is written by `raster` with different WKT text, which
  # makes terra::mask report a CRS mismatch on grids that already align.
  terra::crs(m) <- terra::crs(r)
  dt <- terra::datatype(r)[1]
  # nodata must match the dtype (uint8 class maps such as dominant species and V_* -> 255, the uint16
  # plot index -> 0), or the nodata tag is lost on write. A matching tag also keeps re-filtering
  # idempotent.
  na <- if (grepl("^FLT", dt)) NaN else if (dt == "INT1U") 255 else if (dt == "INT2U") 0 else -128
  nb <- terra::global(r, "notNA")[1, 1]
  r2 <- terra::mask(r, m, maskvalue = 1, updatevalue = NA)
  n_after <- terra::global(r2, "notNA")[1, 1]

  tmp <- paste0(f, ".tmp.tif")
  terra::writeRaster(r2, tmp, overwrite = TRUE, datatype = dt, NAflag = na,
                     gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
  rm(r, r2); gc(verbose = FALSE)
  if (file.exists(f)) stopifnot(file.remove(f))
  stopifnot(file.rename(tmp, f))
  for (ext in c(".ovr", ".aux.xml")) if (file.exists(paste0(f, ext))) file.remove(paste0(f, ext))
  cat(sprintf("  %-34s removed %s cells (%s ha)\n", basename(f),
              format(nb - n_after, big.mark = ","), format(round((nb - n_after) * CELL_HA), big.mark = ",")))
}

## ---- optional restriction: LC_ONLY = character vector of BASENAMES to filter -------------------
# A run that replaced a few layers need not rewrite the rest. Every rewrite of a GeoTIFF can damage
# its tile index, so re-filter only what changed.
#   LC_ONLY <- c("ftm_conifer_fraction_30m.tif"); source(...)
lc_keep <- function(v) if (exists("LC_ONLY")) v[basename(v) %in% LC_ONLY] else v

## ---- A) direct-RF maps -------------------------------------------------------------------------------
cat("[A] direct-RF maps:\n")
direct_maps <- c(list.files(here("outputs", "nfi8", "maps_30m"),        "\\.tif$", full.names = TRUE),
           list.files(here("outputs", "nfi8", "uncertainty_30m"), "\\.tif$", full.names = TRUE))
direct_maps <- lc_keep(direct_maps[!grepl("\\.tmp\\.tif$", direct_maps)])
for (f in direct_maps) apply_one(f, mask0)

## ---- B) plot-index imputation (the index, its distance layer, the decoded attributes) --------------
cat("[B] plot-index imputation (index/ and attributes_ftm/; mask nearest-resampled where needed):\n")
# The index directory is globbed so that every layer beside the index, such as the nearest-reference
# distance layer, loses the same cells.
imp_files <- c(list.files(here("outputs", "nfi8", "index"),          "\\.tif$", full.names = TRUE),
               list.files(here("outputs", "nfi8", "attributes_ftm"), "^spec30m_.*\\.tif$", full.names = TRUE))
imp_files <- lc_keep(imp_files[file.exists(imp_files) & !grepl("\\.tmp\\.tif$", imp_files)])
# one resampled mask per distinct grid; apply_one re-aligns per file.
mask_cache <- list()
mask_for <- function(r) {
  key <- paste(dim(r)[1:2], collapse = "x")
  if (is.null(mask_cache[[key]]))
    mask_cache[[key]] <<- terra::subst(terra::resample(mask0, r, method = "near"), NA, 0L)
  mask_cache[[key]]
}
for (f in imp_files) apply_one(f, mask_for(terra::rast(f)))

cat("\nDone. Land-cover filter (remove built/barren/water) applied in place to",
    length(direct_maps), "direct +", length(imp_files), "imputation rasters.\n",
    "Overviews dropped -> run pipeline/06d_finalize_rasters.py",
    "(rebuilds the overviews) before figures.\n")
