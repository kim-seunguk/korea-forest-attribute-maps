# 03d_extract_predictors_at_plots.R — sample every raster predictor at the inventory plot locations.
# ------------------------------------------------------------------------------------------------
# This is the join between the plot table (01b) and the covariate rasters (02). It produces the
# training table that stages 04-07 read.
#
# Predictors are taken from the RASTER at the plot, not from the plot's own field record — the
# forest-type attributes included. Training and prediction therefore see the same source of
# information. Reading the field record instead would train the model on something the map cannot
# supply at prediction time.
#
# The plot's own field classification is kept alongside, prefixed fld_, because the imputation
# matches on it and because comparing fld_ against the raster value is itself a reported result.
#
# The forest-type map carries five stand attributes; the three in R/paths.R's V_PREDICTORS are the
# ones read here. The map's own stand age and stand height are not: no model uses either, and
# neither is delivered.
#
# In : data/interim/nfi8/plot_reference.csv
#      data/interim/fmap_30m/<edition>/{V_VG,V_DC,V_DN,V_FROR}.tif
#      data/interim/covar_5179/{C_*,T_*,sin_y,cos_x}.tif
# Out: data/interim/nfi8/rd_points.csv  and  rd_points.gpkg   (EPSG:5179)
#
#   source(here::here("pipeline", "03d_extract_predictors_at_plots.R"))   # run in RStudio, light
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(terra); library(data.table) })
source(here("R", "paths.R"))

PT_CRS    <- "EPSG:5179"                                    # plot UTMK_X/Y are already on this CRS
covar_dir <- here("data", "interim", "covar_5179")
V_LAYERS  <- V_PREDICTORS                                   # the forest-type map attributes models read
COVARS    <- c("C_PREC", "C_TMAX", "C_TMIN", "C_RH", "C_GSL", "C_FD0", "C_DTR",
               "T_ELE", "T_SLO", "sin_y", "cos_x", "C_PAR")
log_t <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# The columns this script writes, in the order it writes them. The checkpoint below compares an
# existing file against this, so that a run after the column set has changed stops instead of
# quietly leaving the old table in place for every later stage to read.
OUT_COLUMNS <- c("Plot_ID", "UTMK_X", "UTMK_Y", "Lon", "Lat",
                 V_LAYERS, "FROR", COVARS,
                 paste0("fld_", V_LAYERS),
                 "dom_species", "dom_ba_share", "BAPH", "Bio_a_ha", "Carbon_a_ha")

out_csv <- here("data", "interim", VERSION, "rd_points.csv")
if (file.exists(out_csv)) {
  have <- names(fread(out_csv, nrows = 0L))
  if (identical(have, OUT_COLUMNS)) {
    log_t("skip (exists):", out_csv)
  } else {
    stop(out_csv, " already exists, but its columns are not the ones this script writes now.",
         "\n  it has  : ", paste(setdiff(have, OUT_COLUMNS), collapse = ", "),
         "\n  it lacks: ", paste(setdiff(OUT_COLUMNS, have), collapse = ", "),
         "\n  Delete it and re-run, so that every later stage reads one table.", call. = FALSE)
  }
} else {
  pr <- fread(here("data", "interim", VERSION, "plot_reference.csv"), encoding = "UTF-8")
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
                     here("data", "interim", VERSION, "rd_points.gpkg"), overwrite = TRUE)
  log_t(sprintf("  wrote %s (+ .gpkg)", out_csv))
}
