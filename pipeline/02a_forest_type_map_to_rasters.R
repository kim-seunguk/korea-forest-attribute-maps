# 02a_forest_type_map_to_rasters.R — regenerate the 30 m stand-attribute rasters (V_VG/V_DC/V_DN)
# from the forest-type map 1:5000 (large-scale forest-type map), per year, nationwide.
#
# The stand attributes the models use as predictors come from the forest-type map, which ships as
# 1:5,000 vector sheets. This rasterizes them onto the delivered 30 m grid, one archive per map
# edition, so that a run can be tied to a specific edition rather than to whatever was on disk.
#
# The map also carries stand age and stand height. Neither is rasterized: no model reads them, and
# no product delivers them.
#
# The source is already on EPSG:5179, so nothing is reprojected and no resampling loss is incurred.
#
# SOURCE (staged under data/raw/; see "STAGING" below):
#   data/raw/forest_type_map/large_scale_5000/<year>/<시도코드>_<시도>.shp   (17 시도 + 4 splits/yr)
#   source: 산림공간정보서비스 (map.forest.go.kr), Korea Forest Service.
#
# OUTPUT:
#   data/interim/fmap_30m/<year>/V_{VG,DC,DN}.tif         national 30 m, EPSG:5179
#   data/interim/fmap_30m/<year>/_tiles/                  per-시도 intermediate tiles (resume cache)
#   data/interim/fmap_30m/<year>/coverage.csv             per-시도 forest-polygon counts (2025 부분판 check)
#
# HEAVY GIS (≈20 GB of vectors per edition) → run in RStudio, one edition at a time. Resumable:
# finished tiles are cached and skipped, so an interrupted run continues where it stopped.
#
#   source(here::here("pipeline","02a_forest_type_map_to_rasters.R"))
#
# ---- STAGING (do once) -------------------------------------------------------------------------
# The source is large and need not be copied. A directory junction or symlink into data/raw/ works:
#   mklink /J "<repo>\data\raw\forest_type_map" "<where you keep the map>"     (Windows)
# (or copy/move the three year folders there). The script only reads it.
# ------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(here); library(terra); library(sf); library(data.table) })
source(here("R", "paths.R"))

## ---- config ------------------------------------------------------------------------------------
YEARS    <- FTM_VINTAGE                        # the edition the dataset uses; add more to bake them
RES      <- 30                                  # output resolution (m)
OUT_CRS  <- "EPSG:5179"                         # = source CRS (UTM-K, Korea 2000 Unified)
SRC_ROOT <- here("data", "raw", "forest_type_map", "large_scale_5000")
OUT_ROOT <- here("data", "interim", "fmap_30m")

# National 5179 extent (generous bbox covering South Korea; fixed so all 시도/years share one grid).
NAT_EXT  <- terra::ext(740000, 1400010, 1440000, 2080020)
template <- terra::rast(NAT_EXT, resolution = RES, crs = OUT_CRS)   # empty master grid

# Attribute source fields (FTM DBF) and their recodes -----------------------------------------
# ALL source fields are CHARACTER, across every edition checked: KOFTR_GROU
# "10".."77", DMCLS_CD "0/1/2/3", DNST_CD "A/B/C".
# CRS auto-reads as EPSG:5179 already, so st_transform below is a no-op safety.
# The ONLY unmapped KOFTR_GROU is 78 = 죽림 (bamboo, FRTP_CD=4) — a few per cent nationally,
# concentrated in the south, with its structural attributes blank. It is left NA: bamboo carries
# no structural attributes and has no counterpart among the plot forest types of 01b. The
# per-province "unmapped V_*" counts logged below are that bamboo, not a failure.
recl_vg <- fread(here("data","forest_type_recode.csv"), encoding = "UTF-8")
vg_map  <- setNames(recl_vg$recl, as.character(recl_vg$KOFTR_GROU))   # 원본코드 -> 10..60 수종군
dm_map  <- c(`0` = 3L, `1` = 12L, `2` = 24L, `3` = 36L)              # DMCLS_CD -> 대표경급(cm) 치수/소/중/대
dn_map  <- c(A = 25L, B = 60L, C = 85L)                              # DNST_CD  -> 수관밀도(%) 소/중/밀
SRC_FIELDS <- c("STORUNST", "FROR_CD", "KOFTR_GROU", "DMCLS_CD", "DNST_CD")
V_LAYERS   <- V_PREDICTORS                                           # V_VG, V_DC, V_DN — see R/paths.R

log_t <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

## ---- per-시도 -> per-attribute tiles -----------------------------------------------------------
# Reads only forest rows + needed cols (OGR-SQL at read time), maps attributes, rasterizes each V_*
# onto the national grid cropped to the 시도 extent (shared origin -> tiles mosaic cleanly).
make_tile <- function(shp, tile_dir) {
  lyr <- tools::file_path_sans_ext(basename(shp))
  done <- file.path(tile_dir, sprintf("%s__%s.tif", lyr, V_LAYERS))
  if (all(file.exists(done))) { log_t(" skip (cached):", lyr); return(data.table(sido = lyr, n = NA_integer_)) }

  # Read the full layer plainly (geometry guaranteed) then filter in R. OGR-SQL `query=` proved
  # unreliable here (drops the shapefile geometry column on this GDAL build). Memory is bounded
  # because we process one 시도 at a time and drop non-forest rows immediately.
  x <- sf::st_read(shp, quiet = TRUE)
  flds <- names(x)                                                    # schema differs 2019=14 vs 2024/25=20
  use  <- intersect(SRC_FIELDS, flds)
  miss <- setdiff(c("KOFTR_GROU","DMCLS_CD","DNST_CD"), use)
  if (length(miss)) log_t("  WARN", lyr, "missing fields ->", paste(miss, collapse = ", "))
  # forest = STORUNST == '1' (입목지, stocked; the polygons that carry species and structure).
  # Other values are left out rather than recoded — their meaning is edition-specific, so check the
  # codebook of the edition you stage. FROR_CD (임종) is not filtered on; STORUNST already isolates
  # stocked forest.
  if ("STORUNST" %in% flds) x <- x[as.character(x$STORUNST) == "1" & !is.na(x$STORUNST), ]
  else log_t("  NOTE", lyr, "no STORUNST field -> no forest filter applied")
  if (nrow(x) && "FROR_CD" %in% flds) {
    n0 <- sum(as.character(x$FROR_CD) == "0", na.rm = TRUE)
    if (n0) log_t(sprintf("  FROR_CD=='0' within forest: %d polys (kept)", n0))
  }
  x <- x[, intersect(use, names(x))]                                  # keep only V_-source columns (+ geometry)
  if (nrow(x) == 0) { log_t("  EMPTY after filter:", lyr); return(data.table(sido = lyr, n = 0L)) }

  x <- sf::st_transform(x, OUT_CRS)
  v <- terra::vect(x)
  v$V_VG <- unname(vg_map[as.character(x$KOFTR_GROU)])
  v$V_DC <- unname(dm_map[as.character(x$DMCLS_CD)])
  v$V_DN <- unname(dn_map[as.character(x$DNST_CD)])

  for (vl in V_LAYERS) {
    src_col <- c(V_VG="KOFTR_GROU", V_DC="DMCLS_CD", V_DN="DNST_CD")[[vl]]
    if (src_col %in% use && anyNA(v[[vl]]))
      log_t(sprintf("  unmapped %s (%s): %d polys", vl, src_col, sum(is.na(v[[vl]]))))
  }

  tile_tmpl <- terra::crop(template, terra::ext(v), snap = "out")     # keep national origin/alignment
  for (vl in V_LAYERS) {
    r <- terra::rasterize(v, tile_tmpl, field = vl)
    terra::writeRaster(r, file.path(tile_dir, sprintf("%s__%s.tif", lyr, vl)),
                       overwrite = TRUE, datatype = "INT2S",
                       gdal = c("COMPRESS=LZW", "TILED=YES"))
  }
  n <- nrow(x); rm(x, v); gc(verbose = FALSE)                         # bound memory across big 시도
  log_t(" tiled:", lyr, "(", n, "forest polys )")
  data.table(sido = lyr, n = n)
}

## ---- per-year mosaic ---------------------------------------------------------------------------
build_year <- function(year) {
  src_dir <- file.path(SRC_ROOT, year)
  require_file(src_dir)
  out_dir <- file.path(OUT_ROOT, year); tile_dir <- file.path(out_dir, "_tiles")
  dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)

  shps <- list.files(src_dir, "\\.shp$", full.names = TRUE)
  shps <- shps[!grepl("\\.shp\\.xml$", shps)]
  log_t(sprintf("=== %s: %d shapefiles ===", year, length(shps)))

  cov <- rbindlist(lapply(shps, make_tile, tile_dir = tile_dir), fill = TRUE)
  fwrite(cov, file.path(out_dir, "coverage.csv"), bom = TRUE)         # forest-poly count per 시도 (부분판 check)

  for (vl in V_LAYERS) {                                              # mosaic 시도 tiles -> national layer
    tiles <- list.files(tile_dir, sprintf("__%s\\.tif$", vl), full.names = TRUE)
    if (!length(tiles)) { log_t("  no tiles for", vl); next }
    m <- if (length(tiles) == 1) terra::rast(tiles[1])
         else terra::mosaic(terra::sprc(lapply(tiles, terra::rast)), fun = "max")
    terra::writeRaster(m, file.path(out_dir, paste0(vl, ".tif")),
                       overwrite = TRUE, datatype = "INT2S",
                       gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
    log_t("  wrote national", vl, "->", file.path(out_dir, paste0(vl, ".tif")))
  }
  log_t(sprintf("=== %s done. coverage.csv: %d 시도, %s forest polys ===",
                year, nrow(cov), format(sum(cov$n, na.rm = TRUE), big.mark = ",")))
}

## ---- run ---------------------------------------------------------------------------------------
terra::terraOptions(progress = 0)
for (y in YEARS) build_year(y)
cat("\nDone. Inspect each year's coverage.csv: a 2025 PARTIAL update shows far fewer / missing 시도\n",
    "than 2024. National V_* rasters: data/interim/fmap_30m/<year>/V_*.tif (EPSG:5179, 30 m).\n",
    "Next: pick the edition that matches your inventory cycle, then\n",
    "03d_extract_predictors_at_plots.R samples V_* at the NFI 표본점 coordinates.\n")
