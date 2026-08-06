# 03f_collect_plot_features.R — assemble the Sentinel-2 feature table from the exported chunks.
# ------------------------------------------------------------------------------------------------
# 03e samples the composites at the plots and exports the result in chunks, because one request
# carrying every plot stalls on upload. This joins the downloaded chunks into the single table that
# stages 05 and 07 train on.
#
# It is not a concatenation. Earth Engine's table export adds its own bookkeeping columns and sorts
# the properties alphabetically, so the 42 features come back in an order that is not the order of
# the raster's bands. Downstream code selects them by name, and this writes them in the raster's
# order so that the two readings agree by eye as well as by name.
#
# In : data/raw/gee_plot_features/plot_features_s2_*.csv   (downloaded from Drive)
# Out: data/interim/nfi8/plot_features_s2.csv              (Plot_ID + the 42 s2_* columns)
#
#   source(here::here("pipeline", "03f_collect_plot_features.R"))   # run in RStudio, light
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("R", "paths.R"))

src_dir <- here("data", "raw", "gee_plot_features")
out_csv <- here("data", "interim", VERSION, "plot_features_s2.csv")

chunks <- sort(list.files(src_dir, "^plot_features_s2_.*\\.csv$", full.names = TRUE))
if (!length(chunks))
  stop("No chunk CSVs in ", src_dir,
       "\n  Run 03e, then download its exports from Drive into that folder.", call. = FALSE)
cat(sprintf("reading %d chunk(s) from %s\n", length(chunks), src_dir))

feat <- rbindlist(lapply(chunks, fread, colClasses = list(character = "Plot_ID")),
                  use.names = TRUE, fill = TRUE)

missing <- setdiff(c("Plot_ID", S2_BAND_NAMES), names(feat))
if (length(missing))
  stop("the chunks are missing ", length(missing), " expected column(s): ",
       paste(head(missing, 6), collapse = ", "), if (length(missing) > 6) " ...",
       "\n  Check that every export task finished before downloading.", call. = FALSE)

feat <- feat[, c("Plot_ID", S2_BAND_NAMES), with = FALSE]   # drops system:index, .geo and anything else
setorder(feat, Plot_ID)

dups <- feat[duplicated(Plot_ID), unique(Plot_ID)]
if (length(dups))
  stop(length(dups), " plot(s) appear in more than one chunk, starting with ", dups[1],
       "\n  Clear the folder and download only the latest run's exports.", call. = FALSE)

n_na <- sum(!complete.cases(feat[, ..S2_BAND_NAMES]))
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(feat, out_csv, bom = TRUE)
cat(sprintf("wrote %s\n  %d plots x %d features | plots with any missing feature: %d\n",
            out_csv, nrow(feat), length(S2_BAND_NAMES), n_na))
