# 02b_forest_origin_to_raster.R — rasterize stand origin, planted against natural.
# ================================================================================================
# 02a reads the origin field but only reports it. This script writes it to a raster on the same
# grid, with the same forest filter and mosaic rules, so it co-registers with the other stand
# attributes. It serves both as a predictor and as the stratum the accuracy is broken down by in 07f.
#
# The 임종 (FROR_CD) codes are written through unchanged: 1 is planted, 2 is natural, and 0 is
# neither.
#
# Out: data/interim/fmap_30m/<edition>/V_FROR.tif   (INT2S, EPSG:5179, 30 m)
#
# Heavy — roughly 20 GB of vectors per edition. Resumable: finished tiles are cached and skipped.
#   source(here::here("pipeline", "02b_forest_origin_to_raster.R"))
# ================================================================================================
suppressPackageStartupMessages({ library(here); library(terra); library(sf); library(data.table) })
source(here("R", "paths.R"))

YEARS    <- FTM_VINTAGE                        # the forest-type map edition the dataset uses
RES      <- 30; OUT_CRS <- "EPSG:5179"
SRC_ROOT <- here("data", "raw", "forest_type_map", "large_scale_5000")
OUT_ROOT <- here("data", "interim", "fmap_30m")
NAT_EXT  <- terra::ext(740000, 1400010, 1440000, 2080020)   # identical to 02a → shared grid/origin
template <- terra::rast(NAT_EXT, resolution = RES, crs = OUT_CRS)
log_t <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

make_tile <- function(shp, tile_dir) {
  lyr <- tools::file_path_sans_ext(basename(shp))
  out <- file.path(tile_dir, sprintf("%s__V_FROR.tif", lyr))
  if (file.exists(out)) { log_t(" skip (cached):", lyr); return(invisible()) }
  x <- sf::st_read(shp, quiet = TRUE); flds <- names(x)
  if (!("FROR_CD" %in% flds)) { log_t("  WARN no FROR_CD field ->", lyr); return(invisible()) }
  if ("STORUNST" %in% flds) x <- x[as.character(x$STORUNST) == "1" & !is.na(x$STORUNST), ]  # 입목지 forest filter (= 02a)
  else log_t("  NOTE", lyr, "no STORUNST -> no forest filter")
  x <- x[, intersect("FROR_CD", names(x))]
  if (nrow(x) == 0) { log_t("  EMPTY after filter:", lyr); return(invisible()) }
  x <- sf::st_transform(x, OUT_CRS); v <- terra::vect(x)
  v$V_FROR <- suppressWarnings(as.integer(as.character(x$FROR_CD)))   # RAW code 0/1/2 (no recode)
  tt <- terra::crop(template, terra::ext(v), snap = "out")            # keep national origin/alignment
  r  <- terra::rasterize(v, tt, field = "V_FROR")
  terra::writeRaster(r, out, overwrite = TRUE, datatype = "INT2S", gdal = c("COMPRESS=LZW", "TILED=YES"))
  log_t(" tiled:", lyr, "(", nrow(x), "forest polys )"); rm(x, v); gc(verbose = FALSE)
}

build_year <- function(year) {
  src <- file.path(SRC_ROOT, year); require_file(src)
  out <- file.path(OUT_ROOT, year); td <- file.path(out, "_tiles"); dir.create(td, recursive = TRUE, showWarnings = FALSE)
  shps <- list.files(src, "\\.shp$", full.names = TRUE); shps <- shps[!grepl("\\.shp\\.xml$", shps)]
  log_t(sprintf("=== %s: %d shapefiles (FROR) ===", year, length(shps)))
  invisible(lapply(shps, make_tile, tile_dir = td))
  tiles <- list.files(td, "__V_FROR\\.tif$", full.names = TRUE)
  if (!length(tiles)) { log_t("  no FROR tiles for", year); return(invisible()) }
  m <- if (length(tiles) == 1) terra::rast(tiles[1]) else terra::mosaic(terra::sprc(lapply(tiles, terra::rast)), fun = "max")
  terra::writeRaster(m, file.path(out, "V_FROR.tif"), overwrite = TRUE, datatype = "INT2S",
                     gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
  log_t("  wrote national V_FROR ->", file.path(out, "V_FROR.tif"))
}

terra::terraOptions(progress = 0); for (y in YEARS) build_year(y)
cat("\nDone. V_FROR.tif per edition (EPSG:5179, 30 m, raw FROR_CD 0/1/2).\n",
    "Stage 03d extracts it at the plot locations as the FROR column.\n",
    "! The codes are written through unchanged, so a different edition's coding survives as-is.\n")
