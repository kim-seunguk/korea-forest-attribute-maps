# 05c_index_to_attributes.R — decode the plot-index raster into per-cell attribute rasters.
# ------------------------------------------------------------------------------------------------
# Stage 05b assigns every forest cell the index of one matched inventory plot. This script
# substitutes each index with that plot's measured values, so every attribute raster is a field
# observation carried to the cell rather than a model prediction.
#
# The imputation matched on the plot's field forest type, diameter class and density class. The
# attributes written here are read off the matched plot; they never took part in selecting it, so
# species, biomass and diversity all arrive as one internally consistent combination.
#
# In : outputs/nfi8/index/imputation{_base,_ftm}.tif
#      outputs/nfi8/models/train_reference{_base,_ftm}.csv   (plot_idx 1:1 with the raster)
# Out: outputs/nfi8/attributes{_base,_ftm}/spec30m_*.tif
#
# Run in RStudio, after 05b. Light — minutes, not hours.
#   IMPUTE_TIER <- "base"   # then "full"
#   source(here::here("pipeline", "05c_index_to_attributes.R"))
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(terra); library(data.table) })
source(here::here("R", "paths.R"))

# Which predictor set the index raster came from; mirrors 05b.
if (!exists("IMPUTE_TIER")) IMPUTE_TIER <- "full"
stopifnot(IMPUTE_TIER %in% c("base", "full"))
.tsuf <- tier_suffix(IMPUTE_TIER)                        # "_base" | "_ftm"

outd <- here("outputs", VERSION)
log_time <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
log_time(sprintf("index -> attributes | IMPUTE_TIER = '%s'", IMPUTE_TIER))

index_path <- require_file(file.path(outd, "index", sprintf("imputation%s.tif", .tsuf)))
out_dir    <- file.path(outd, sprintf("attributes%s", .tsuf))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
im <- terra::rast(index_path)
terra::crs(im) <- "EPSG:5179"          # carry the EPSG code onto the attribute rasters

stamp <- function(vals, nm, legend = NULL) {             # substitute plot_idx -> value, then write
  r <- terra::subst(im, from = seq_along(vals), to = vals, others = NA); names(r) <- nm
  terra::writeRaster(r, file.path(out_dir, paste0("spec30m_", nm, ".tif")), overwrite = TRUE)
  if (!is.null(legend)) fwrite(legend, file.path(out_dir, paste0("spec30m_", nm, "_legend.csv")), bom = TRUE)
  log_time(paste0("  wrote spec30m_", nm, ".tif"))
}

ref <- fread(file.path(outd, "models", sprintf("train_reference%s.csv", .tsuf)), encoding = "UTF-8")
setorder(ref, plot_idx)                                  # plot_idx must be 1..N in raster order
stopifnot(ref$plot_idx == seq_len(nrow(ref)))
cat("plots indexed:", nrow(ref), "\n")

# Field stand structure of the matched plot: forest type, age, diameter and density class.
for (v in V_STAND_ATTRIBUTES) stamp(as.numeric(ref[[v]]), v)

# Dominant species becomes an integer raster plus a legend, since raster values must be numeric.
sp_levels <- sort(unique(ref$dom_species))
stamp(match(ref$dom_species, sp_levels), "dom_species",
      legend = data.table(code = seq_along(sp_levels), species = sp_levels))

stamp(as.numeric(ref$dom_ba_share),     "dom_ba_share")   # basal-area share of the dominant species
stamp(as.numeric(ref$Bio_a_ha),         "Bio_a_ha")       # aboveground biomass, Mg/ha
stamp(as.numeric(ref$Carbon_a_ha),      "Carbon_a_ha")    # aboveground carbon, Mg/ha
stamp(as.numeric(ref$BAPH),             "BAPH")           # basal area per hectare
stamp(as.numeric(ref$shannon_full),     "shannon")        # tree-species diversity, Shannon H'
stamp(as.numeric(ref$coniferfrac_full), "coniferfrac")    # coniferous fraction of basal area

log_time("done")
cat("\nAttribute rasters in", out_dir, "\n")
