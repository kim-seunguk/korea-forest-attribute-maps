# 02e_landcover_to_raster.R — rasterize the national land-cover map to its 7 top-level classes
# onto the delivered 30 m map grid (EPSG:5179), for the land-cover filter that 06b applies.
#
# The delivered maps take their forest extent from the forest-type map (STORUNST == 1), which
# commits some non-forest — airports, urban blocks — to forest, so the models predict forest values
# over it. The MCEE 세분류 토지피복지도 is an independent land-cover product and shows where that
# happened. The source is the Level-3 map at 1:5,000. Each polygon also holds its Level-1 group code,
# and that code is what is rasterized:
#   100 시가화·건조(built) / 200 농업(agri) / 300 산림(forest) / 400 초지(grass) /
#   500 습지(wetland)     / 600 나지(barren) / 700 수역(water)
# It is pixel-aligned to the delivered maps, so 06a_landcover_crosstab.R can tabulate the delivered
# forest extent against the MCEE class, and 06b then prunes the cells whose land-cover class alone
# marks them as errors. All seven classes are kept, so that the other direction — MCEE forest
# outside the delivered extent — stays visible too.
#
# SOURCE (staged under data/raw/landcover_mcee_national/, or reached by a junction or symlink):
#   one zip per map sheet, SG05_<도엽번호>_<yyyymmdd>.zip, each a Level-3 1:5,000 sheet,
#   EPSG:5186 (Korea 2000 / Central Belt), .dbf CP949, field L1_CODE = '100'..'700'.
#   Read straight out of the zip via GDAL /vsizip — the 30 GB is never fully extracted.
#
# OUTPUT:
#   data/interim/landcover_mcee_30m/landcover_L1_30m.tif   national 30 m INT2U, 0 = nodata, EPSG:5179
#   data/interim/landcover_mcee_30m/_tiles/<도엽번호>.tif       per-sheet resume cache (skipped if present)
#   data/interim/landcover_mcee_30m/tiles_index.csv             per-sheet poly count (coverage check)
#
# Heavy, and bound by disk rather than CPU: it reads one shapefile per map sheet, tens of thousands
# of them. Resumable — finished sheets are cached and skipped.
#
# ! Run after stage 05: the grid template is a delivered raster. Nothing before stage 06 reads the
#   output.
#
#   source(here::here("pipeline","02e_landcover_to_raster.R"))

suppressPackageStartupMessages({ library(here); library(terra); library(sf); library(data.table) })
source(here("R", "paths.R"))

## ---- config ------------------------------------------------------------------------------------
SRC_ROOT  <- here("data", "raw", "landcover_mcee_national")                 # the staged source
OUT_DIR   <- here("data", "interim", "landcover_mcee_30m")
TILE_DIR  <- file.path(OUT_DIR, "_tiles")
OUT_TIF   <- file.path(OUT_DIR, "landcover_L1_30m.tif")
# Any delivered raster works as the template: they all share one grid.
TMPL_TIF  <- here("outputs", "nfi8", "maps_30m", "ftm_conifer_fraction_30m.tif")

dir.create(TILE_DIR, recursive = TRUE, showWarnings = FALSE)
require_file(SRC_ROOT); require_file(TMPL_TIF)
template <- terra::rast(TMPL_TIF); terra::values(template) <- NA         # empty grid, delivered geometry
out_crs  <- terra::crs(template)

log_t <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

## ---- per-sheet rasterize -----------------------------------------------------------------------
# Each sheet's L1_CODE is rasterized onto the template cropped to the sheet extent. The crop keeps
# the delivered-grid origin, so the sheet tiles mosaic cell for cell.
make_tile <- function(zip) {
  id   <- sub("^SG05_([0-9]+)_.*\\.zip$", "\\1", basename(zip))          # inner shp basename = 도엽번호
  out  <- file.path(TILE_DIR, paste0(id, ".tif"))
  if (file.exists(out)) return(data.table(id = id, n = NA_integer_))     # resume

  vpath <- sprintf("/vsizip/%s/%s.shp", gsub("\\\\", "/", zip), id)
  x <- tryCatch(sf::st_read(vpath, quiet = TRUE, options = "ENCODING=CP949"),
                error = function(e) { log_t("  READ FAIL", id, conditionMessage(e)); NULL })
  if (is.null(x) || nrow(x) == 0) return(data.table(id = id, n = 0L))

  x$L1 <- suppressWarnings(as.integer(as.character(x$L1_CODE)))
  x <- x[!is.na(x$L1), "L1"]
  if (nrow(x) == 0) return(data.table(id = id, n = 0L))

  x <- sf::st_transform(x, out_crs)
  v <- terra::vect(x)
  tmpl <- terra::crop(template, terra::ext(v), snap = "out")             # keep delivered-grid origin/alignment
  r <- terra::rasterize(v, tmpl, field = "L1", background = NA)
  terra::writeRaster(r, out, overwrite = TRUE, datatype = "INT2U", NAflag = 0,
                     gdal = c("COMPRESS=LZW", "TILED=YES"))
  data.table(id = id, n = nrow(x))
}

## ---- run ---------------------------------------------------------------------------------------
terra::terraOptions(progress = 0)
zips <- list.files(SRC_ROOT, pattern = "\\.zip$", full.names = TRUE)
log_t(sprintf("=== %d land-cover sheets -> 30 m L1 raster (delivered grid) ===", length(zips)))

idx <- vector("list", length(zips))
for (i in seq_along(zips)) {
  idx[[i]] <- make_tile(zips[i])
  if (i %% 500 == 0) { gc(verbose = FALSE); log_t(sprintf("  %d / %d sheets", i, length(zips))) }
}
fwrite(rbindlist(idx, fill = TRUE), file.path(OUT_DIR, "tiles_index.csv"), bom = TRUE)

## ---- assemble national mosaic ------------------------------------------------------------------
tiles <- list.files(TILE_DIR, "\\.tif$", full.names = TRUE)
log_t(sprintf("=== assembling %d sheet tiles -> national L1 raster ===", length(tiles)))
vfile <- file.path(OUT_DIR, "landcover_L1.vrt")
vr <- terra::vrt(tiles, filename = vfile, overwrite = TRUE)             # scalable mosaic (no full load)
vr <- terra::extend(vr, template)                                      # pad to the full delivered extent (NA fill)
terra::writeRaster(vr, OUT_TIF, overwrite = TRUE, datatype = "INT2U", NAflag = 0,
                   gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
log_t("wrote", OUT_TIF)

cat("\nDone. National Level-1 land cover on the delivered 30 m grid:\n  ", OUT_TIF,
    "\n  values 100 built / 200 agri / 300 forest / 400 grass / 500 wetland / 600 barren / 700 water",
    "\n  (0 = nodata / sea). Next: source 06a_landcover_crosstab.R for the forest x land-cover table.\n")
