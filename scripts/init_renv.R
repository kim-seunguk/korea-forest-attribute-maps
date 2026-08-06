# init_renv.R — pin the packages the pipeline depends on, into renv.lock.
# Run in RStudio from the project root. May prompt interactively.
#
# ! This is the script that PRODUCES the lockfile, not the one that consumes it. To reproduce the
#   pinned environment, call renv::restore() instead — running this re-pins against whatever happens
#   to be installed on the machine it runs on.

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# Core dependencies (install before snapshot if missing):
pkgs <- c(
  "here",          # project-root-relative paths
  "terra", "sf",   # raster and vector GIS, used throughout
  "raster", "sp",  # the older raster classes: a RasterStack serialises to the parallel workers in
                   # stage 05, where a terra SpatRaster cannot
  "yaImpute",      # random-forest nearest-neighbour imputation (the TreeMap engine)
  "gower",         # Gower distance — REQUIRED by yaImpute::newtargets (a Suggests, not auto-installed)
  "randomForest",  # the random-forest backend, for both product families
  "foreach", "doParallel", "doSNOW",   # the parallel harness in stages 05 and 07
  "data.table",    # the inventory tree records are large
  "dplyr", "purrr",        # 01b
  "ggplot2", "patchwork"   # figures/
)

# Create the renv structure without auto-discovery the first time.
if (!file.exists("renv/activate.R")) {
  renv::init(bare = TRUE, restart = FALSE)
}

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
}

renv::snapshot(packages = pkgs, prompt = FALSE)
message("renv.lock updated. Commit it alongside the code.")
