# 03c_split_sentinel2_bands.R — split the 42-band S2 stack into 42 single-band tifs for stage 05b's fast-read path.
#
# Stage 05 reads the predictors in narrow blocks of raster rows, and pulling one band out of a
# 42-band file costs far more than reading a single-band file, because the bands are interleaved and
# the window drags the others along with it. Splitting once removes that cost from every block of
# the wall-to-wall run.
#
# ! Skipping this does not fail. 05b falls back to the merged file and simply becomes too slow to
#   finish. Run it once, after 03b and before 05b.
#
# INPUT : data/interim/covar_5179_s2/s2med42_5179.tif        (42-band, built by 03b_align_sentinel2.R)
# OUTPUT: data/interim/covar_5179_s2/s2_bands/s2band_01..42.tif  (single-band, physical band order)
#   source(here::here("pipeline","03c_split_sentinel2_bands.R"))
suppressPackageStartupMessages({ library(here); library(terra) })
source(here("R", "paths.R"))
SRC <- require_file(here("data", "interim", "covar_5179_s2", "s2med42_5179.tif"))
DST <- here("data", "interim", "covar_5179_s2", "s2_bands"); dir.create(DST, showWarnings = FALSE, recursive = TRUE)
r <- terra::rast(SRC)
stopifnot(terra::nlyr(r) == 42L)                                  # 42 bands = the 05b contract
for (i in seq_len(42L)) {
  f <- file.path(DST, sprintf("s2band_%02d.tif", i))
  terra::writeRaster(r[[i]], f, overwrite = TRUE,
                     gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
  cat(sprintf("  s2band_%02d.tif  (%s)\n", i, names(r)[i]))
}
cat(sprintf("wrote 42 single-band tifs -> %s\n", DST))
