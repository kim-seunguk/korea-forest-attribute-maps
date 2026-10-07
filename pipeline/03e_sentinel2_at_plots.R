# 03e_sentinel2_at_plots.R — read the Sentinel-2 training features at the plots from the aligned raster.
# ------------------------------------------------------------------------------------------------
# Stage 05 predicts onto the Sentinel-2 raster that 03a exported and 03b aligned. The models are
# trained on the same raster: this samples its 42 bands at the plot coordinates of rd_points.csv and
# writes the plot feature table that stages 05 and 07 read. Training and prediction therefore see
# one set of numbers by construction — the value at a plot is the value of the cell it falls in, not
# a second composite sampled some other way — and nothing here touches Earth Engine.
#
# The bands are read from 03c's single-band files, the fast path, or from 03b's merged file when 03c
# has not been run. Either way they are labelled in the raster's physical band order, which
# R/paths.R holds as S2_BAND_NAMES; any other labelling would scramble every predictor.
#
# The table carries Plot_ID and the 42 bands and nothing else. Stand origin (FROR) and the other
# raster predictors are already in rd_points.csv from 03d, and stages 05 and 07 join the two tables
# on Plot_ID, so a column present in both would collide there.
#
# In : data/interim/covar_5179_s2/s2_bands/s2band_01..42.tif   (03c)
#      or data/interim/covar_5179_s2/s2med42_5179.tif          (03b)
#      data/interim/nfi8/rd_points.csv                          (03d: Plot_ID, UTMK_X, UTMK_Y)
# Out: data/interim/nfi8/plot_features_s2.csv                   (Plot_ID + the 42 s2_* columns)
# An existing table is moved to data/interim/nfi8/_superseded_<date>/ first, never overwritten.
#
#   source(here::here("pipeline", "03e_sentinel2_at_plots.R"))   # run in RStudio, light
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(data.table); library(terra) })
source(here("R", "paths.R"))

s2_dir   <- here("data", "interim", "covar_5179_s2")
band_dir <- file.path(s2_dir, "s2_bands")
out_csv  <- here("data", "interim", VERSION, "plot_features_s2.csv")

s2_files <- if (dir.exists(band_dir)) {
  file.path(band_dir, sprintf("s2band_%02d.tif", seq_len(42L)))          # 03c's single-band files
} else {
  file.path(s2_dir, "s2med42_5179.tif")                                  # 03b's merged file
}
for (f in s2_files) require_file(f)
r <- terra::rast(s2_files)
stopifnot(terra::nlyr(r) == 42L)
names(r) <- S2_BAND_NAMES                                                 # the physical band order; see R/paths.R

rd <- fread(require_file(here("data", "interim", VERSION, "rd_points.csv")),
            encoding = "UTF-8", colClasses = list(character = "Plot_ID"))
stopifnot(all(c("Plot_ID", "UTMK_X", "UTMK_Y") %in% names(rd)), !anyDuplicated(rd$Plot_ID))
xy <- as.matrix(rd[, .(UTMK_X, UTMK_Y)])                                  # EPSG:5179 metres, the raster's CRS
cat(sprintf("extracting 42 Sentinel-2 bands at %s plots\n", format(nrow(xy), big.mark = ",")))

e    <- terra::extract(r, xy)                                             # rows follow rd's order
feat <- data.table(Plot_ID = rd$Plot_ID)
feat[, (S2_BAND_NAMES) := as.data.table(e[, S2_BAND_NAMES, drop = FALSE])]

n_complete <- sum(complete.cases(feat[, S2_BAND_NAMES, with = FALSE]))
cat(sprintf("plots with a complete 42-band record: %s / %s (%.1f%%)\n",
            format(n_complete, big.mark = ","), format(nrow(feat), big.mark = ","),
            100 * n_complete / nrow(feat)))
if (n_complete < 0.5 * nrow(feat))
  warning("fewer than half the plots got Sentinel-2 values; check the export's forest mask and extent.",
          call. = FALSE)

if (file.exists(out_csv)) {                                               # keep the previous table, never overwrite it
  a <- file.path(dirname(out_csv), sprintf("_superseded_%s", format(Sys.Date(), "%Y%m%d")))
  dir.create(a, showWarnings = FALSE)
  stopifnot(file.rename(out_csv, file.path(a, basename(out_csv))))
  cat("archived the previous table ->", basename(a), "\n")
}
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(feat, out_csv, bom = TRUE)
cat(sprintf("wrote %s\n  %d plots x %d features\n  Next: pipeline/04a_spatial_folds.R\n",
            out_csv, nrow(feat), length(S2_BAND_NAMES)))
