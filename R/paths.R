# R/paths.R — where the pipeline finds and writes things.
#
# Every path is resolved with here::here() from the repository root, so scripts run wherever the
# repository is checked out.
#
# Layout:
#   data/                      small lookup tables, tracked with the code
#   data/raw/<source>/         source data you supply (not in this repository)
#   data/interim/              regridded covariates and the plot table, written by stages 01-03
#   outputs/nfi8/              models, maps, validation tables — everything stages 04-08 produce
#
# Usage:
#   library(here)
#   source(here::here("R", "paths.R"))
#   rd <- data.table::fread(require_file(here::here("data", "interim", VERSION, "rd_points.csv")))

# Assert a file exists before reading it, with a message that says what to do about it.
require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path,
         "\n  Stage the source data under data/raw/ first.",
         call. = FALSE)
  }
  path
}

# ---- reference configuration ----------------------------------------------------------------
# The inventory cycle this dataset is built from. It names the output directory, outputs/nfi8/.
VERSION <- "nfi8"

# The forest-type map edition matched to the inventory period. One edition is used throughout.
FTM_VINTAGE <- "2025"

# The same three stand attributes in two roles:
#   V_PREDICTORS        read from the forest-type map by the models
#   V_STAND_ATTRIBUTES  carried from the matched plot into the delivered lookup table
V_PREDICTORS       <- c("V_VG", "V_DC", "V_DN")
V_STAND_ATTRIBUTES <- V_PREDICTORS

# ---- the Sentinel-2 band contract ---------------------------------------------------------------
# The 42 composite bands in the physical band order of the exported raster: spring's fourteen, then
# summer's, then fall's. These names label the raster, and 03e labels the plot feature table with
# them, so any other order mislabels every band. Code reading the table selects columns by name.
# The Earth Engine definition is pipeline/s2_composite.py.
#
# NDWI here is the NIR-SWIR1 moisture index (B8, B11) = NDMI (Wilson & Sader 2002), not Gao (1996)
# NDWI. The band is labeled NDWI in the pipeline and published as NDMI.
S2_BAND_NAMES <- as.vector(t(outer(
  c("spring", "summer", "fall"),
  c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12", "NDVI", "NDWI", "NDRE", "EVI"),
  function(season, band) paste0("s2_", season, "_", band))))

# ---- product file naming ----------------------------------------------------------------------
# Rasters are written as ftm_<attr>_30m.tif; 08d renames them to the published names.
#
# FILE_VAR maps the inventory column that keys a target in the analysis to its file name:
# `Bio_a_ha` is the column, `biomass` is the file.
#
# The products are built from one predictor set — the run parameter `full` — and published under
# the file token `ftm`, after the forest-type map it includes. FILE_PREDSET translates one to the other.
FILE_PREDSET <- c(full = "ftm")
FILE_VAR  <- c(coniferfrac = "conifer_fraction", Bio_a_ha = "biomass",
               shannon = "shannon", domsp = "dominant_species")

# predset: "full" (the run parameter)   var: a FILE_VAR key   kind: map | sd | uncertainty
product_file <- function(predset, var, kind = "map") {
  t <- unname(FILE_PREDSET[predset]); v <- unname(FILE_VAR[var])
  if (is.na(t) || is.na(v)) stop("unknown predset/var: ", predset, " / ", var, call. = FALSE)
  switch(kind,
         map         = sprintf("%s_%s_30m.tif", t, v),
         sd          = sprintf("%s_%s_sd_30m.tif", t, v),
         uncertainty = sprintf("%s_%s_uncertainty_30m.tif", t, v),
         stop("unknown kind: ", kind, call. = FALSE))
}

# The plot-index family takes the same token as an explicit suffix.
predset_suffix <- function(predset) {
  t <- unname(FILE_PREDSET[predset])
  if (is.na(t)) stop("unknown predset: ", predset, call. = FALSE)
  paste0("_", t)
}
