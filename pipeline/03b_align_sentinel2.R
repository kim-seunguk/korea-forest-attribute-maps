# 03b_align_sentinel2.R — align the Sentinel-2 export onto the covariate grid that stage 05 predicts on.
#
# 03a exports the composite on Earth Engine's own 30 m grid, which is not the covariate grid the
# imputation runs on. This mosaics the downloaded tiles and resamples them onto the covariate
# template, bilinearly since the bands are continuous, so Sentinel-2 stacks cell for cell with the
# climate, terrain and stand-attribute layers in 05b.
#
# INPUT : data/raw/gee_s2med/m8_s2med42_30m*.tif   (1+ tiles from Drive; 42 bands, forest-masked)
# OUTPUT: data/interim/covar_5179_s2/s2med42_5179.tif  (42-band, covar_5179 grid, band names s2_<season>_<band>)
#   source(here::here("pipeline","03b_align_sentinel2.R"))

suppressPackageStartupMessages({ library(here); library(terra) })
source(here("R", "paths.R"))

SRC_DIR <- here("data", "raw", "gee_s2med")
TMPL    <- here("data", "interim", "covar_5179", "C_PREC.tif")        # covar_5179 grid template
OUT_DIR <- here("data", "interim", "covar_5179_s2"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_TIF <- file.path(OUT_DIR, "s2med42_5179.tif")

S2_NAMES <- S2_BAND_NAMES        # the 42 band names, season-major; defined in R/paths.R

tiles <- sort(list.files(SRC_DIR, "^m8_s2med42_30m.*\\.tif$", full.names = TRUE))
if (!length(tiles)) stop("No S2 tiles in ", SRC_DIR,
                         " — download m8_s2med42_30m*.tif from Drive first.", call. = FALSE)
# Tiles are mosaicked by georeference, so their numbering does not matter, but a leftover tile from an
# earlier export would be silently mixed in. Clear the folder before re-downloading; the spread of file
# times below is the warning if you did not.
mt <- file.info(tiles)$mtime
cat(sprintf("mosaicking %d S2 tile(s) from %s:\n", length(tiles), SRC_DIR))
for (i in seq_along(tiles)) cat(sprintf("  %-46s (%s)\n", basename(tiles[i]), format(mt[i], "%Y-%m-%d %H:%M")))
if (length(tiles) > 1 && as.numeric(diff(range(mt)), units = "hours") > 6)
  warning("S2 tiles span >6 h in mtime — an earlier export's tiles may be mixed in. ",
          "Clear ", SRC_DIR, " and re-download ONLY the latest export.", call. = FALSE)
# The export names its bands, so check them rather than assume them: the names below are reassigned
# after resampling, and if the export ever came back in a different order they would be attached to
# the wrong data with nothing to show for it. The check runs on the tiles because a terra VRT does
# not carry band names through, which also makes it a real test of the mtime warning above — a tile
# from another export fails here rather than being guessed at.
for (t in tiles) {
  nm <- names(terra::rast(t))
  if (!identical(nm, S2_NAMES)) {
    i <- which(nm != S2_NAMES)[1]
    stop(basename(t), " does not carry the 42 expected band names in order.",
         if (length(nm) != 42L) paste0("\n  it has ", length(nm), " bands, not 42")
         else paste0("\n  band ", i, ": expected ", S2_NAMES[i], ", found ", nm[i]),
         "\n  Re-export with 03a, and download only that export's tiles.", call. = FALSE)
  }
}

if (length(tiles) == 1) {
  s2 <- terra::rast(tiles[1])
} else {
  s2 <- terra::vrt(tiles, filename = file.path(OUT_DIR, "s2_src.vrt"), overwrite = TRUE)
}

tmpl <- terra::rast(require_file(TMPL))
s2a  <- terra::project(s2, tmpl, method = "bilinear")                 # same CRS -> resample/align to covar_5179
names(s2a) <- S2_NAMES
terra::writeRaster(s2a, OUT_TIF, overwrite = TRUE,
                   gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
cat("wrote", OUT_TIF, "\n  42 bands on the covar_5179 grid; names:", paste(head(S2_NAMES, 3), collapse = ", "),
    "...", tail(S2_NAMES, 1), "\n  Next: pipeline/05b_plot_index_imputation.R\n")
