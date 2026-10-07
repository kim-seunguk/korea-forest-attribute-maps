# R/paths.R — where the pipeline finds and writes things.
#
# Every path is resolved with here::here() from the repository root, so scripts run wherever the
# repository is checked out. Nothing is absolute and nothing calls setwd().
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
#   pts <- sf::read_sf(require_file(here::here("data", "raw", "nfi_points", "plots.shp")))

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

# The stand attributes, and the two different roles they play.
#
# V_PREDICTORS are the three the models read from the forest-type map. Stand age is deliberately
# absent: it is largely implied by diameter class, and including it would add a variable the map
# supplies less reliably than the field does.
#
# V_STAND_ATTRIBUTES are the three the plot-index route carries out of the matched plot and into the
# delivered lookup table. They are the same three, because the matching axis is what that product
# can speak to. Stand age is not delivered: the field class is coarse and is not one of the
# attributes this dataset stands behind. A user who needs age can derive it from the plot record
# that nfi_plot_id resolves to.
V_PREDICTORS       <- c("V_VG", "V_DC", "V_DN")
V_STAND_ATTRIBUTES <- V_PREDICTORS

# ---- the Sentinel-2 band contract ---------------------------------------------------------------
# The 42 composite bands, in the order the export builds them: spring's fourteen, then summer's,
# then fall's. That is the physical band order of the raster, so it is what labels the raster, and
# it must not be inferred from anything else. The definition on the Earth Engine side is
# pipeline/s2_composite.py.
#
# 03e reads the plot feature table off that same raster and labels it with these names, so the
# table and the raster agree by construction; code reading the table still selects its columns by
# name. Assigning any other order to the raster mislabels every band.
#
# NDWI here is the NIR-SWIR1 moisture index (B8, B11) = NDMI (Wilson & Sader 2002), not Gao (1996)
# NDWI. The internal label stays NDWI to match the existing tables and models; it is reported and
# deposited as NDMI.
S2_BAND_NAMES <- as.vector(t(outer(
  c("spring", "summer", "fall"),
  c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12", "NDVI", "NDWI", "NDRE", "EVI"),
  function(season, band) paste0("s2_", season, "_", band))))

# ---- product file naming ----------------------------------------------------------------------
# Rasters are written under the names the dataset is published with.
#
# FILE_VAR translates because the two sides name things differently: the analysis keys a target by
# its inventory column, which cannot be renamed, while the file needs a name a reader can read.
# `Bio_a_ha` is the column, `biomass` is the file.
#
# The products are built from one predictor set — the run parameter `full` — and published under
# the file token `ftm`, after the forest-type map it includes. FILE_TIER translates one to the other.
FILE_TIER <- c(full = "ftm")
FILE_VAR  <- c(coniferfrac = "conifer_fraction", Bio_a_ha = "biomass",
               shannon = "shannon", domsp = "dominant_species")

# tier: "full" (the run parameter)   var: a FILE_VAR key   kind: map | sd | uncertainty
product_file <- function(tier, var, kind = "map") {
  t <- unname(FILE_TIER[tier]); v <- unname(FILE_VAR[var])
  if (is.na(t) || is.na(v)) stop("unknown tier/var: ", tier, " / ", var, call. = FALSE)
  switch(kind,
         map         = sprintf("%s_%s_30m.tif", t, v),
         sd          = sprintf("%s_%s_sd_30m.tif", t, v),
         uncertainty = sprintf("%s_%s_uncertainty_30m.tif", t, v),
         stop("unknown kind: ", kind, call. = FALSE))
}

# The plot-index family takes the same token as an explicit suffix, so that a filename never has to
# be read as "the one without a suffix".
tier_suffix <- function(tier) {
  t <- unname(FILE_TIER[tier])
  if (is.na(t)) stop("unknown tier: ", tier, call. = FALSE)
  paste0("_", t)
}
