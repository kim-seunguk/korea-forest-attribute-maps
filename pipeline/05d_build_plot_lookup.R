# 05d_build_plot_lookup.R — the table the plot-index rasters are read through.
# ------------------------------------------------------------------------------------------------
# The plot-index rasters carry an integer, and this is what the integer means. One row per matched
# reference plot, in raster order, holding that plot's attributes under the names the dataset is
# published with. Without it the index rasters cannot be interpreted at all.
#
# The table is a renaming of the reference library 05b writes. It is kept as its own stage because
# it is a deposited product and the published column names are decided here, not inside a modelling
# script.
#
# What is deliberately absent:
#   coordinates  the plot locations are confidential and are never distributed. A user who needs
#                more than this table joins nfi_plot_id to the open NFI survey database.
#   age class    the field class is coarse, and it is not one of the attributes this dataset stands
#                behind. It is available through nfi_plot_id for anyone who wants it.
#
# In : outputs/nfi8/models/train_reference_{base,ftm}.csv
# Out: outputs/nfi8/plot_lookup_{base,ftm}.csv
#
# The deposit copy is made by 08d, which harmonises the species columns to the published 39-class
# coding and drops the byte-order mark. Species here is the field name, as the imputation matched it.
#
#   source(here::here("pipeline", "05d_build_plot_lookup.R"))   # run in RStudio, light
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("R", "paths.R"))

# published name <- reference-library column. Order is the published column order.
COLUMNS <- c(plot_index                = "plot_idx",
             nfi_plot_id               = "Plot_ID",
             coniferous_fraction       = "coniferfrac_full",
             dominant_species          = "dom_species",
             aboveground_biomass_Mgha  = "Bio_a_ha",
             shannon_diversity         = "shannon_full",
             forest_type               = "V_VG",
             diameter_class            = "V_DC",
             crown_density             = "V_DN",
             dominant_species_ba_share = "dom_ba_share",
             carbon_Mgha               = "Carbon_a_ha",
             basal_area_m2ha           = "BAPH")

outd <- here("outputs", VERSION)

for (tier in c("base", "full")) {
  tsuf <- tier_suffix(tier)
  ref  <- fread(require_file(file.path(outd, "models", sprintf("train_reference%s.csv", tsuf))),
                encoding = "UTF-8")
  missing <- setdiff(unname(COLUMNS), names(ref))
  if (length(missing)) stop("train_reference", tsuf, " lacks: ", paste(missing, collapse = ", "),
                            call. = FALSE)

  setorder(ref, plot_idx)
  # The index is the raster value, so it has to be 1..N with nothing missing or repeated.
  stopifnot(identical(ref$plot_idx, seq_len(nrow(ref))))

  out <- ref[, unname(COLUMNS), with = FALSE]
  setnames(out, names(COLUMNS))
  fn <- file.path(outd, sprintf("plot_lookup%s.csv", tsuf))
  fwrite(out, fn, bom = TRUE)
  cat(sprintf("%-22s %d plots x %d columns\n", basename(fn), nrow(out), ncol(out)))
}

cat("\nNext: 08d stages these into the deposit with the species columns harmonised.\n")
