# 05c_index_to_attributes.R — decode the plot-index raster into per-cell attribute rasters.
# ------------------------------------------------------------------------------------------------
# Stage 05b assigns every forest cell the index of one matched inventory plot. This script
# substitutes each index with that plot's measured values, so every attribute raster is a field
# observation carried to the cell rather than a model prediction.
#
# The imputation matched on the plot's field forest type, diameter class and density class. The
# attributes written here never took part in selecting that plot.
#
# In : outputs/nfi8/index/imputation_ftm.tif
#      outputs/nfi8/models/train_reference_ftm.csv   (plot_idx 1:1 with the raster)
# Out: outputs/nfi8/attributes_ftm/spec30m_*.tif
#
# Run in RStudio, after 05b. It is light and runs in minutes.
#   source(here::here("pipeline", "05c_index_to_attributes.R"))
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(terra); library(data.table) })
source(here::here("R", "paths.R"))

.psuf <- predset_suffix("full")                          # "_ftm", the suffix 05b wrote; see R/paths.R

outd <- here("outputs", VERSION)
log_time <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
log_time("index -> attributes")

index_path <- require_file(file.path(outd, "index", sprintf("imputation%s.tif", .psuf)))
out_dir    <- file.path(outd, sprintf("attributes%s", .psuf))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
im <- terra::rast(index_path)
terra::crs(im) <- "EPSG:5179"          # carry the EPSG code onto the attribute rasters

stamp <- function(vals, nm, legend = NULL) {             # substitute plot_idx -> value, then write
  r <- terra::subst(im, from = seq_along(vals), to = vals, others = NA); names(r) <- nm
  terra::writeRaster(r, file.path(out_dir, paste0("spec30m_", nm, ".tif")), overwrite = TRUE)
  if (!is.null(legend)) fwrite(legend, file.path(out_dir, paste0("spec30m_", nm, "_legend.csv")), bom = TRUE)
  log_time(paste0("  wrote spec30m_", nm, ".tif"))
}

ref <- fread(file.path(outd, "models", sprintf("train_reference%s.csv", .psuf)), encoding = "UTF-8")
setorder(ref, plot_idx)                                  # plot_idx must be 1..N in raster order
stopifnot(ref$plot_idx == seq_len(nrow(ref)))
cat("plots indexed:", nrow(ref), "\n")

# Field stand structure of the matched plot: forest type, diameter and density class.
for (v in V_STAND_ATTRIBUTES) stamp(as.numeric(ref[[v]]), v)

# Dominant species becomes an integer raster plus a legend, since raster values must be numeric. This
# is the plot's field species name, unpooled, so the legend here is the inventory's own list and is
# longer than the direct map's: 07j puts the two on one basis before comparing them.
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
