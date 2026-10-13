# 05d_build_plot_lookup.R — the table the plot-index raster is read through.
# ------------------------------------------------------------------------------------------------
# The plot-index raster carries an integer, and this is what the integer means. One row per matched
# reference plot, in raster order, holding that plot's attributes under the names the dataset is
# published with. Without it the index raster cannot be interpreted at all.
#
# The table is a renaming of the reference library 05b writes. The published column names are set
# here.
#
# Not in the table:
#   coordinates  the plot locations are confidential and are never distributed. A user who needs
#                more than this table joins nfi_plot_id to the inventory's own records.
#   age class    not delivered; it is available through nfi_plot_id.
#
# In : outputs/nfi8/models/train_reference_ftm.csv
# Out: outputs/nfi8/plot_lookup_ftm.csv   (12 columns)
#
# The deposit copy is made by 08d, which drops the byte-order mark and adds two species columns
# beside the field name — the code on the trained legend and the scientific name — for the 14
# published columns. Species here is the field name, as the imputation matched it.
#
# The values are copied as text, never parsed and re-formatted, and the byte-order mark, line ending
# and quoting are pinned, so the file is byte-identical whenever it is rebuilt from the same reference
# library on any machine. 08d relies on that: it rebuilds the table from the reference file and
# compares the two row by row, cell by cell (as text, or as numbers where the text differs).
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
tsuf <- tier_suffix("full")                                           # "_ftm"; see R/paths.R
ref  <- fread(require_file(file.path(outd, "models", sprintf("train_reference%s.csv", tsuf))),
              encoding = "UTF-8", colClasses = "character")           # text in, text out: no re-formatting
missing <- setdiff(unname(COLUMNS), names(ref))
if (length(missing)) stop("train_reference", tsuf, " lacks: ", paste(missing, collapse = ", "),
                          call. = FALSE)

# The index is the raster value, so the rows must already be 1..N in order: the row order carries the
# meaning and is never sorted or filtered here.
stopifnot(identical(as.integer(ref$plot_idx), seq_len(nrow(ref))))
# 05b leaves plots without biomass or carbon out of the library, so no row may be blank there.
stopifnot(all(nzchar(ref$Bio_a_ha)), all(nzchar(ref$Carbon_a_ha)))

out <- ref[, unname(COLUMNS), with = FALSE]
setnames(out, names(COLUMNS))
fn <- file.path(outd, sprintf("plot_lookup%s.csv", tsuf))
fwrite(out, fn, bom = TRUE, eol = "\r\n", quote = "auto", na = "")
cat(sprintf("%-22s %d plots x %d columns\n", basename(fn), nrow(out), ncol(out)))

cat("\nNext: 08d stages this into the deposit with the species columns added.\n")
