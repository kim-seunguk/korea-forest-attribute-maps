# 02c_climate_terrain_to_rasters.R — the climate (C_*) and terrain (T_*) covariate rasters, on the
# SAME EPSG:5179 / 30 m grid as the 02a V_* layers.
#
# Sibling of 02a, which rasterizes the forest-type map. Together they produce the full covariate
# stack; 03d then samples it at the 표본점 coordinates, and stage 05 predicts on the grid itself.
#
# SOURCES:
#   C_PREC/C_TMAX/C_TMIN/C_RH  <- KMA MK-PRISM v2.1 yearly txt, 2000-2019 mean   [data/raw/climate_kma]
#   C_GSL/C_FD0/C_DTR          <- MK-PRISM v2.1 bioclimatic indices               [data/raw/climate_kma]
#   T_ELE/T_SLO + sin_y/cos_x  <- SRTM elevation, exported onto this grid by 02f  [data/raw/dem_srtm_grid]
#   Lon/Lat                    <- planar EPSG:5179 cell coordinates (derived)
#   C_PAR                      <- GK2A satellite solar x 2.02 (W/m2 -> umol/m2/s) [data/raw/solar_gk2a]
#
# OUTPUT: data/interim/covar_5179/<LAYER>.tif   (EPSG:5179, 30 m, aligned to 02a V_*)
# CHECKPOINTED: a layer whose .tif already exists is SKIPPED (a re-run builds only what is missing),
# unless it is named in REBUILD; the old file is then moved to covar_5179/_superseded_<date>/ first,
# never deleted.
#   REBUILD <- c("T_ELE", "T_SLO", "sin_y", "cos_x")     # e.g. after a new DEM export, before source()
#
#   source(here::here("pipeline","02c_climate_terrain_to_rasters.R"))   # run in RStudio (HEAVY)
#
# ---- STAGING (do once; data/raw/ is not tracked, so nothing is fetched for you) -----------------
#   DEM (SRTM):   run 02f_srtm_to_grid.py and download its export into data/raw/dem_srtm_grid/
#   C_PAR (GK2A): hardlink the 60 monthly .nc into data/raw/solar_gk2a/  (supplied by you)
#   A hard link avoids a second copy when the source is already on the same volume.
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(terra) })
source(here::here("R", "paths.R"))

## ---- grid template = the 02a V_* grid, extended north ------------------------------------------
## The forest-type map ends at about 38.37N, which is as far as it was surveyed, while the DEM and
## the climate grids hold data further north. The template keeps the map's cell-for-cell alignment
## and raises its north edge to COVAR_NORTH_LAT, so the covariates reach as far as their sources do.
## This extended grid is the product grid: 02f and 03a export onto it, and every raster the pipeline
## writes shares it. The delivered products cover the forest-type map's stocked forest (03a masks the
## Sentinel-2 composite to it), so cells north of the map hold covariates but no prediction. Far
## offshore islands stay NA: the climate grids have no data there.
v_ref <- here("data", "interim", "fmap_30m", FTM_VINTAGE, "V_VG.tif")   # the edition set in R/paths.R
if (!file.exists(v_ref)) stop("Run 02a first — missing ", v_ref)
COVAR_NORTH_LAT <- 38.80                                           # north of both the map and the climate grids
.v0  <- terra::rast(v_ref)
.py  <- terra::geom(terra::project(terra::vect(cbind(128.4, COVAR_NORTH_LAT), crs = "EPSG:4326"),
                                   "EPSG:5179"))[1, "y"]           # 5179 y of the target north lat (east-coast lon)
.e   <- terra::ext(.v0); .e[4] <- .py                             # raise ONLY the north edge (ymax)
template <- terra::rast(terra::extend(.v0, .e, snap = "out"))      # geometry-only, cell-aligned, NA rows added on top
cat(sprintf("[covar template] north edge %.3f -> %.3fN (5179 y=%.0f); grid %d x %d (was %d rows)\n",
            terra::ymax(.v0), COVAR_NORTH_LAT, .py, terra::ncol(template), terra::nrow(template), terra::nrow(.v0)))
cat(sprintf("[covar template] west edge %.10f, north edge %.1f, cell %.15f x %.1f m — the GRID_TRANSFORM of 02f and 03a\n",
            terra::xmin(template), terra::ymax(template), terra::res(template)[1], terra::res(template)[2]))

out_dir <- here("data", "interim", "covar_5179"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
kma_dir <- here("data", "raw", "climate_kma")
log_t   <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
if (!exists("REBUILD")) REBUILD <- character(0)
FILL_PASSES <- 3L                                              # coastal fill of the climate grids, rings of cells; see the reader
have    <- function(nm) file.exists(file.path(out_dir, paste0(nm, ".tif"))) && !(nm %in% REBUILD)   # checkpoint
arch_dir <- file.path(out_dir, sprintf("_superseded_%s", format(Sys.Date(), "%Y%m%d")))
archive_old <- function(nm) {                                  # move, never delete: the old layer stays readable
  olds <- list.files(out_dir, pattern = sprintf("^%s[.]tif([.](ovr|aux[.]xml))?$", nm), full.names = TRUE)
  if (!length(olds)) return(invisible())
  dir.create(arch_dir, showWarnings = FALSE)
  for (f in olds) stopifnot(file.rename(f, file.path(arch_dir, basename(f))))
  log_t(sprintf("  archived old %s -> %s", nm, basename(arch_dir)))
}

write_layer <- function(r, nm) {
  names(r) <- nm
  if (nm %in% REBUILD) archive_old(nm)
  terra::writeRaster(r, file.path(out_dir, paste0(nm, ".tif")), overwrite = TRUE,
                     datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
  rng <- round(as.numeric(terra::global(r, "range", na.rm = TRUE)), 2)
  log_t(sprintf("  wrote %-7s range [%s, %s]", nm, rng[1], rng[2]))
}

## ---- MK-PRISM yearly text reader ---------------------------------------------------------------
# Each yearly grid = 20 lines (2000-2019), each line = the 751x601 grid flattened (byrow), NA=-9990.
read_mkprism_mean <- function(elem, nx = 751L, ny = 601L,
                              xmin = 124.5, ymin = 33.0, cell = 0.01, na_val = -9990, nyears = 20L) {
  path  <- file.path(kma_dir, sprintf("MKPRISM_MKPRISMv21_skorea_%s_gridraw_yearly_2000_2019.txt", elem))
  if (!file.exists(path)) stop("missing climate txt: ", path)
  vals  <- scan(path, what = numeric(), quiet = TRUE)
  ncell <- nx * ny; stopifnot(length(vals) == ncell * nyears)
  arr <- array(NA_real_, dim = c(ny, nx, nyears))
  for (i in seq_len(nyears)) arr[, , i] <- matrix(vals[((i - 1) * ncell + 1):(i * ncell)],
                                                  nrow = ny, ncol = nx, byrow = TRUE)
  arr <- arr[ny:1, , , drop = FALSE]                       # source rows run north to south
  arr[arr == na_val] <- NA
  xmax <- xmin + cell * (nx - 1); ymax <- ymin + cell * (ny - 1)   # the last grid point, not an edge
  # The grid points are cell centres, so the raster reaches half a cell beyond the first and last
  # point on every side. Taking the points as the outer edges would stretch the cells and shift the
  # data by up to half a cell at the far edge.
  e <- terra::ext(xmin - cell / 2, xmax + cell / 2, ymin - cell / 2, ymax + cell / 2)
  s <- terra::rast(arr, extent = e, crs = "EPSG:4326")
  stopifnot(abs(terra::res(s)[1] - cell) < 1e-9, abs(terra::res(s)[2] - cell) < 1e-9)
  m <- terra::app(s, fun = mean, na.rm = TRUE)              # 2000-2019 normal -> 1 layer (4326, ~1 km)
  # Coastal and island cells whose nearest climate cell is NA would otherwise get no climate at all,
  # and the plots on them would drop out of every model. NA cells next to data take the mean of their
  # valid neighbours, one ring of cells per pass.
  for (k in seq_len(FILL_PASSES)) m <- terra::focal(m, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
  m
}
to_5179 <- function(r) terra::project(r, template, method = "bilinear")   # ~1 km -> 5179/30 m

## ---- 1. climate C_* + bioclimatic indices (MK-PRISM) ----------------------------------------
log_t("=== climate (MK-PRISM v2.1 yearly, 2000-2019 mean) ===")
clim <- list(C_PREC = "RN", C_TMAX = "TAMAX", C_TMIN = "TAMIN", C_RH = "RHM",
             C_GSL = "GSL", C_FD0 = "FD0", C_DTR = "DTR")
for (nm in names(clim)) if (have(nm)) log_t("  skip (exists)", nm) else
  write_layer(to_5179(read_mkprism_mean(clim[[nm]])), nm)

## ---- 2. terrain T_* (SRTM, already on this grid) ------------------------------------------------
t_names <- c("T_ELE", "T_SLO", "sin_y", "cos_x")
if (all(vapply(t_names, have, logical(1)))) {
  log_t("=== terrain: skip (all exist) ===")
} else {
  log_t("=== terrain (SRTM on the product grid) ===")
  # 02f samples the native SRTM bilinearly, once, straight onto this grid, so nothing is reprojected
  # here: slope and aspect are derived from the elevation exactly as exported.
  dem_dir <- here("data", "raw", "dem_srtm_grid")
  tiles   <- list.files(dem_dir, pattern = "^srtm_elev_grid30m.*[.]tif$", full.names = TRUE)
  if (!length(tiles)) stop("Stage the grid-aligned SRTM first (02f_srtm_to_grid.py; see STAGING header): ", dem_dir)
  dem <- if (length(tiles) == 1L) terra::rast(tiles) else terra::vrt(tiles)
  if (terra::compareGeom(dem, template, stopOnError = FALSE, messages = FALSE)) {
    T_ELE <- dem                                              # already on the grid: nothing is resampled
  } else {                                                    # same lattice, other extent: crop or pad, still no resampling
    stopifnot(isTRUE(all.equal(terra::res(dem), terra::res(template))),
              abs(((terra::xmin(dem) - terra::xmin(template)) / terra::res(template)[1]) %% 1) < 1e-6)
    T_ELE <- terra::extend(terra::crop(dem, template, snap = "near"), template)
  }
  stopifnot(terra::compareGeom(T_ELE, template, stopOnError = FALSE, messages = FALSE))
  if (!have("T_ELE")) write_layer(T_ELE, "T_ELE")
  if (!have("T_SLO")) write_layer(terra::terrain(T_ELE, v = "slope", unit = "degrees"), "T_SLO")
  if (!have("sin_y") || !have("cos_x")) {
    asp  <- terra::terrain(T_ELE, v = "aspect", unit = "degrees")           # 0-360; terra gives flat cells 90
    flat <- terra::terrain(T_ELE, v = "slope", unit = "degrees") == 0       # no aspect on flat ground
    # Aspect enters as its sine and cosine so that north is continuous across 0 and 360 degrees. A
    # flat cell has no aspect, so both terms are 0 there rather than the sine and cosine of 90.
    # The layer names are historical and read like coordinate transforms; they are terrain.
    s_asp <- terra::mask(sin(pi / 180 * asp), flat, maskvalues = 1, updatevalue = 0)
    c_asp <- terra::mask(cos(pi / 180 * asp), flat, maskvalues = 1, updatevalue = 0)
    if (!have("sin_y")) write_layer(s_asp, "sin_y")
    if (!have("cos_x")) write_layer(c_asp, "cos_x")
  }
}

## ---- 3. planar coordinates Lon/Lat (EPSG:5179) -------------------------------------------------
if (have("Lon")) log_t("  skip (exists) Lon") else write_layer(terra::init(template, "x"), "Lon")
if (have("Lat")) log_t("  skip (exists) Lat") else write_layer(terra::init(template, "y"), "Lat")

## ---- 4. C_PAR (GK2A satellite solar -> PAR) ----------------------------------------------------
# The radiation covariate, from the GK2A satellite solar-radiation product in data/raw/solar_gk2a.
#
# Method: the mean of 60 monthly daytime-mean shortwave fields (W/m2, 2020-2024), converted to a
# photon-flux-density equivalent by a fixed factor of 2.02 umol/J. That factor is a shortwave-to-PAR
# energy fraction of 0.44 times a PAR photon conversion of 4.57 umol/J (McCree 1972).
#
# The files carry a degenerate LCC definition that has to be discarded and replaced: WGS84 Lambert
# conformal conic, standard parallels 30 and 60, origin 38N/126E, extent +/-900000 m at 2 km. The
# raw array is stored bottom-up, so it is flipped vertically on read.
if (have("C_PAR")) {
  log_t("  skip (exists) C_PAR")
} else {
  log_t("=== C_PAR (GK2A 2020-2024 daytime-mean SW -> PAR umol) ===")
  gk_crs <- "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=38 +lon_0=126 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  read_gk2a <- function(f) {
    r <- terra::flip(suppressWarnings(terra::rast(f)), "vertical")
    terra::ext(r) <- terra::ext(-900000, 900000, -900000, 900000); terra::crs(r) <- gk_crs; r
  }
  gk_files <- list.files(here("data", "raw", "solar_gk2a"), pattern = "\\.nc$", full.names = TRUE)
  if (length(gk_files) < 60) stop("Stage GK2A first (expect 60 .nc in data/raw/solar_gk2a): found ",
                                  length(gk_files))
  sw_mean <- terra::app(terra::rast(lapply(gk_files, read_gk2a)), "mean", na.rm = TRUE)  # W/m2 daytime mean
  write_layer(terra::project(sw_mean, template, method = "bilinear") * 2.02, "C_PAR")    # -> umol/m2/s
}

log_t("=== done. covar_5179/ layers (EPSG:5179, 30 m) aligned to 02a V_*. ===")
cat("Next: 03d_extract_predictors_at_plots.R — sample V_* and covar_5179 at the NFI 표본점 coordinates.\n")
