# 02d_forest_mask.R — rasterize the forest-type map's stocking field onto the delivered grid.
# ================================================================================================
# One national mask that the Sentinel-2 export in 03a and the predictions in 05 are both clipped to,
# so every product covers the same cells.
#
# The field equals 1 for stocked forest, which is the mapped extent. Other values are written
# through unchanged: their meaning is edition-specific, so check your codebook. No value at all
# means the map does not reach there, which is not the same as it saying there is no forest.
#
# ! This raster exists to be uploaded, and the upload is manual. Earth Engine cannot read a local
#   file, so 03a can only clip its export to this extent once the mask is an Image asset there.
#   Upload it, set EE_FTM_ASSET to the asset id. Nothing else in the pipeline reads this file.
#
# In : data/raw/forest_type_map/large_scale_5000/<edition>/*.shp   (override with FTM_SRC)
# Out: outputs/nfi8/ftm<edition>_storunst_30m.tif                  (INT1U, EPSG:5179, 30 m)
#
# Heavy. Rasterized province by province and merged, to bound memory.
#   source(here::here("pipeline", "02d_forest_mask.R"))
# ================================================================================================
suppressPackageStartupMessages({ library(here); library(terra); library(sf) })
source(here::here("R", "paths.R"))

FIELD <- "STORUNST"                                   # 1 = stocked forest; other values left as read
if (!exists("FTM_SRC"))
  FTM_SRC <- here("data", "raw", "forest_type_map", "large_scale_5000", FTM_VINTAGE)
if (!dir.exists(FTM_SRC))
  stop("Missing ", FTM_SRC, "\n  Stage the forest-type map there, or set FTM_SRC.", call. = FALSE)

outd  <- here("outputs", VERSION); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
outfn <- file.path(outd, sprintf("ftm%s_storunst_30m.tif", FTM_VINTAGE))

# The template is the stand-attribute grid from 02a, which the covariates in 02c were also built on.
# Taking it from there rather than declaring an extent keeps the mask cell-aligned with them.
V_REF <- here("data", "interim", "fmap_30m", FTM_VINTAGE, "V_VG.tif")
if (!file.exists(V_REF)) V_REF <- here("data", "interim", "covar_5179", "T_ELE.tif")   # same grid
tmpl <- rast(V_REF)[[1]]
tmpl[] <- NA_integer_

shps <- list.files(FTM_SRC, pattern = "\\.shp$", full.names = TRUE)
stopifnot(length(shps) > 0)
acc <- tmpl
for (s in shps) {
  v <- tryCatch(vect(s), error = function(e) NULL)
  if (is.null(v)) { message("skip (read fail): ", basename(s)); next }
  if (!(FIELD %in% names(v))) { message("skip (no ", FIELD, "): ", basename(s)); next }
  if (!same.crs(v, tmpl)) v <- project(v, crs(tmpl))          # already 5179; guard anyway
  v[[FIELD]] <- as.integer(v[[FIELD]][[1]])
  r <- rasterize(v, tmpl, field = FIELD, background = NA)
  acc <- merge(acc, r)                                        # provinces are disjoint, so merge fills NoData
  message("rasterized ", basename(s), "  (", sum(!is.na(values(r))), " cells)")
}
acc <- as.int(acc)
writeRaster(acc, outfn, overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
cat(sprintf("wrote %s\n", outfn))
cat("value counts:\n"); print(freq(acc))
cat("\nNext: upload this to Earth Engine as an Image asset and set EE_FTM_ASSET to it, so that\n",
    "03a_sentinel2_composites.py masks its export to the same extent.\n")
