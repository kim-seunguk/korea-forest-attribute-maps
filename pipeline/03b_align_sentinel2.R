# 03b_align_sentinel2.R — put the Sentinel-2 export on the covariate grid that stage 05 predicts on.
#
# 03a exports the composite on the covariate grid itself, so the tiles it writes already align cell
# for cell with the climate, terrain and stand-attribute layers. This mosaics the downloaded tiles,
# checks that they sit on that grid, labels the 42 bands, rebuilds EVI and writes one
# band-interleaved file. A tile set that is not on the grid — an export made with a scale instead
# of the grid transform — is resampled onto it bilinearly, since the bands are continuous.
#
# INPUT : data/raw/gee_s2med/s2_medians_grid30m*.tif    (1+ tiles from Drive; 42 bands, forest-masked)
# OUTPUT: data/interim/covar_5179_s2/s2med42_5179.tif  (42-band, covar_5179 grid, band names s2_<season>_<band>)
#   source(here::here("pipeline","03b_align_sentinel2.R"))                     # run in RStudio (HEAVY)
#   To read the tiles from another folder, set S2_RAW_DIR <- "<folder>" before source().

suppressPackageStartupMessages({ library(here); library(terra) })
source(here("R", "paths.R"))

SRC_DIR <- if (exists("S2_RAW_DIR")) S2_RAW_DIR else here("data", "raw", "gee_s2med")
TMPL    <- here("data", "interim", "covar_5179", "C_PREC.tif")        # covar_5179 grid template
OUT_DIR <- here("data", "interim", "covar_5179_s2"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_TIF <- file.path(OUT_DIR, "s2med42_5179.tif")

S2_NAMES <- S2_BAND_NAMES        # the 42 band names, season-major; defined in R/paths.R

tiles <- sort(list.files(SRC_DIR, "^s2_medians_grid30m.*\\.tif$", full.names = TRUE))
if (!length(tiles)) stop("No S2 tiles in ", SRC_DIR,
                         " — download s2_medians_grid30m*.tif from Drive first.", call. = FALSE)
# Tiles are mosaicked by georeference, so their numbering does not matter, but a leftover tile from an
# earlier export would be silently mixed in. Clear the folder before re-downloading; the spread of file
# times below is the warning if you did not.
mt <- file.info(tiles)$mtime
cat(sprintf("mosaicking %d S2 tile(s) from %s:\n", length(tiles), SRC_DIR))
for (i in seq_along(tiles)) cat(sprintf("  %-46s (%s)\n", basename(tiles[i]), format(mt[i], "%Y-%m-%d %H:%M")))
if (length(tiles) > 1 && as.numeric(diff(range(mt)), units = "hours") > 6)
  warning("S2 tiles span >6 h in mtime — an earlier export's tiles may be mixed in. ",
          "Clear ", SRC_DIR, " and re-download ONLY the latest export.", call. = FALSE)
# Check the band names each tile carries. The names are reassigned after the mosaic, so an export
# in a different band order would otherwise mislabel every band without an error. The check runs
# on each tile because a terra VRT drops band names, and so it also stops on a leftover tile from
# another export whose band names differ.
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
names(s2) <- S2_NAMES

# ---- EVI from the cell-mean reflectances --------------------------------------------------------
# The export's EVI band is the mean of the 10 m EVI values in each cell. Here each season's EVI is
# rebuilt from the cell's mean B2, B4 and B8 on the integer reflectance scale (L = 10000), so the
# EVI every later stage reads is the index of the reflectances beside it in the stack.
band   <- function(season, b) s2[[paste0("s2_", season, "_", b)]]
layers <- lapply(seq_len(terra::nlyr(s2)), function(i) s2[[i]])
for (season in c("spring", "summer", "fall")) {
  N <- band(season, "B8"); R <- band(season, "B4"); B <- band(season, "B2")
  evi <- 2.5 * (N - R) / (N + 6 * R - 7.5 * B + 10000)
  names(evi) <- paste0("s2_", season, "_EVI")
  layers[[match(names(evi), S2_NAMES)]] <- evi
}
s2 <- terra::rast(layers)
stopifnot(identical(names(s2), S2_NAMES))
cat("EVI rebuilt from the cell-mean B2/B4/B8 for spring, summer and fall\n")

tmpl <- terra::rast(require_file(TMPL))
# Written straight to disk and band-interleaved, so 03c reads one band without reading the other 41.
gdal_opts <- c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "INTERLEAVE=BAND")
if (terra::compareGeom(s2, tmpl, stopOnError = FALSE, messages = FALSE)) {
  s2a <- terra::writeRaster(s2, OUT_TIF, overwrite = TRUE, names = S2_NAMES, gdal = gdal_opts)
  cat("tiles already on the covariate grid: written without resampling\n")
} else {
  s2a <- terra::project(s2, tmpl, method = "bilinear", filename = OUT_TIF, overwrite = TRUE,
                        wopt = list(names = S2_NAMES, gdal = gdal_opts))
  cat("tiles on another lattice: resampled bilinearly onto the covariate grid\n")
}
cat("wrote", OUT_TIF, "\n  42 bands on the covar_5179 grid; names:", paste(head(S2_NAMES, 3), collapse = ", "),
    "...", tail(S2_NAMES, 1), "\n  Next: pipeline/03c_split_sentinel2_bands.R\n")
