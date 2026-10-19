# 03d_extract_predictors_at_plots.R — sample every raster predictor at the inventory plot locations.
# ------------------------------------------------------------------------------------------------
# This is the join between the plot table (01b) and the covariate rasters (02). It produces the
# training table that stages 04-07 read.
#
# Predictors are taken from the RASTER at the plot, not from the plot's own field record — the
# forest-type attributes included. Training and prediction therefore see the same source of
# information.
#
# The plot's own field classification is kept alongside, prefixed fld_, because the imputation
# matches on it and 07b compares it against the raster value.
#
# In : data/interim/nfi8/plot_reference.csv
#      data/interim/fmap_30m/<edition>/{V_VG,V_DC,V_DN,V_FROR}.tif
#      data/interim/covar_5179/{C_*,T_*,sin_y,cos_x}.tif
# Out: data/interim/nfi8/rd_points.csv  and  rd_points.gpkg   (EPSG:5179)
#
#   source(here::here("pipeline", "03d_extract_predictors_at_plots.R"))   # run in RStudio, light
#   With REBUILD_RD <- TRUE set before source(), the script moves an existing table to
#   data/interim/nfi8/_superseded_<date>/ and rebuilds it. Rebuild after re-running 01b: this
#   script copies plot_reference's field classes and coordinates, so a stale plot_reference
#   propagates into every later stage.
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(terra); library(data.table) })
source(here("R", "paths.R"))

PT_CRS    <- "EPSG:5179"                                    # plot UTMK_X/Y are already on this CRS
covar_dir <- here("data", "interim", "covar_5179")
V_LAYERS  <- V_PREDICTORS                                   # the forest-type map attributes models read
COVARS    <- c("C_PREC", "C_TMAX", "C_TMIN", "C_RH", "C_GSL", "C_FD0", "C_DTR",
               "T_ELE", "T_SLO", "sin_y", "cos_x", "C_PAR")
log_t <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
if (!exists("REBUILD_RD")) REBUILD_RD <- FALSE

# OUT_COLUMNS lists the columns this script writes, in order. The checkpoint below stops on an
# existing file whose columns differ, so a stale table cannot pass silently to the later stages.
OUT_COLUMNS <- c("Plot_ID", "UTMK_X", "UTMK_Y", "Lon", "Lat",
                 V_LAYERS, "FROR", COVARS,
                 paste0("fld_", V_LAYERS),
                 "dom_species", "dom_ba_share", "BAPH", "Bio_a_ha", "Carbon_a_ha")

out_csv  <- here("data", "interim", VERSION, "rd_points.csv")
out_gpkg <- here("data", "interim", VERSION, "rd_points.gpkg")
if (file.exists(out_csv) && REBUILD_RD) {                   # keep the previous table, never overwrite it
  a <- here("data", "interim", VERSION, sprintf("_superseded_%s", format(Sys.Date(), "%Y%m%d")))
  dir.create(a, showWarnings = FALSE)
  for (f in c(out_csv, out_gpkg)) if (file.exists(f)) stopifnot(file.rename(f, file.path(a, basename(f))))
  log_t("archived the previous rd_points ->", basename(a))
}
if (file.exists(out_csv)) {
  have <- names(fread(out_csv, nrows = 0L))
  if (identical(have, OUT_COLUMNS)) {
    log_t("skip (exists):", out_csv)
  } else {
    stop(out_csv, " already exists, but its columns are not the ones this script writes now.",
         "\n  it has  : ", paste(setdiff(have, OUT_COLUMNS), collapse = ", "),
         "\n  it lacks: ", paste(setdiff(OUT_COLUMNS, have), collapse = ", "),
         "\n  Set REBUILD_RD <- TRUE and re-run, so that every later stage reads one table.", call. = FALSE)
  }
} else {
  # Plot_ID is read as character: read as a number, the 13-digit id loses its exact value and no
  # join matches it any more.
  pr <- fread(here("data", "interim", VERSION, "plot_reference.csv"), encoding = "UTF-8",
              colClasses = list(character = "Plot_ID"))
  log_t(sprintf("%s: %d plots (forest-type map %s)", VERSION, nrow(pr), FTM_VINTAGE))

  pts <- terra::vect(as.data.frame(pr), geom = c("UTMK_X", "UTMK_Y"), crs = PT_CRS, keepgeom = TRUE)

  # --- stand attributes from the period-matched forest-type map --------------------------------
  vg_dir <- here("data", "interim", "fmap_30m", FTM_VINTAGE)
  vstk   <- terra::rast(file.path(vg_dir, paste0(V_LAYERS, ".tif"))); names(vstk) <- V_LAYERS
  exV    <- terra::extract(vstk, pts, ID = FALSE)

  # --- stand origin (planted or natural), used as a predictor and as a stratum ------------------
  fror   <- terra::rast(file.path(vg_dir, "V_FROR.tif")); names(fror) <- "FROR"
  exF    <- terra::extract(fror, pts, ID = FALSE)

  # --- climate, terrain and position ------------------------------------------------------------
  cstk   <- terra::rast(file.path(covar_dir, paste0(COVARS, ".tif"))); names(cstk) <- COVARS
  exC    <- terra::extract(cstk, pts, ID = FALSE)

  # --- coverage check: plots must land on the grid ----------------------------------------------
  inext <- mean(!is.na(exV$V_VG))
  log_t(sprintf("  forest-type raster hit rate = %.1f%% (low means wrong CRS or off-grid points)", 100 * inext))
  if (inext < 0.5) warning("fewer than half the plots hit the forest-type raster; check coordinates and CRS")

  rd <- data.table(
    Plot_ID = pr$Plot_ID, UTMK_X = pr$UTMK_X, UTMK_Y = pr$UTMK_Y,
    Lon = pr$UTMK_X, Lat = pr$UTMK_Y                         # planar 5179 coordinates, used as predictors
  )
  rd <- cbind(rd, as.data.table(exV))                        # raster stand attributes
  rd <- cbind(rd, as.data.table(exF))                        # FROR
  rd <- cbind(rd, as.data.table(exC))                        # climate, terrain, position
  rd[, `:=`(fld_V_VG = pr$V_VG,                              # the plot's own field classification
            fld_V_DC = pr$V_DC, fld_V_DN = pr$V_DN)]
  rd[, `:=`(dom_species = pr$dom_species, dom_ba_share = pr$dom_ba_share,
            BAPH = pr$BAPH, Bio_a_ha = pr$Bio_a_ha, Carbon_a_ha = pr$Carbon_a_ha)]

  na_x <- sum(!complete.cases(rd[, c(V_LAYERS, COVARS), with = FALSE]))
  log_t(sprintf("  rows %d | plots with any missing predictor = %d", nrow(rd), na_x))

  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  fwrite(rd, out_csv, bom = TRUE)
  terra::writeVector(terra::vect(as.data.frame(rd), geom = c("UTMK_X", "UTMK_Y"), crs = PT_CRS),
                     out_gpkg, overwrite = TRUE)
  log_t(sprintf("  wrote %s (+ .gpkg)", out_csv))
}
