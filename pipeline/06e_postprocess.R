# 06e_postprocess.R — summary statistics and area totals for the finished maps.
# ============================================================================================
# Post-processes the wall-to-wall maps (outputs/nfi8/maps_<RES>/): clip to the national boundary, per-layer
# sanity stats, and area/stock totals with an area-agreement check against the NFI design-based estimate.
# RES picks which map set to summarize. The delivered set is 30 m.
#   RES <- "30m"; source(here::here("pipeline","06e_postprocess.R"))
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(terra); library(data.table); library(sf) })
source(here::here("R", "paths.R"))
VERSION <- "nfi8"
if (!exists("RES"))       RES       <- "30m"
if (!exists("ENSURE_NODATA")) ENSURE_NODATA <- TRUE # set NoData=NaN on the delivered tifs (a float export can omit the tag)
mapdir <- here("outputs", VERSION, sprintf("maps_%s", RES)); stopifnot(dir.exists(mapdir))
outd   <- here("outputs", VERSION)
PREDSETS <- "full"; REG <- c("coniferfrac","shannon","Bio_a_ha")  # "full" = the delivered predictor set (R/paths.R::FILE_PREDSET)

# ---- NFI design-based reference (KFS 산림기본통계, 2020 basis; whole Korea) --------------------
# Design-based (sample-expansion) area estimates, the leakage-free reference for the agreement check.
# Source = Korea Forest Service 2021 (산림기본통계) / NIFoS 2020.
# Replace with the NFI-8 values when KFS publishes them. NB the mapped "conifer-dominant (cf>0.5)"
# area is only APPROXIMATELY the NFI 침엽수림 class (majority-conifer vs stand-type definition) ->
# report it as broad area agreement only.
NFI_FOREST_HA  <- 6298000     # official total forest area (ha), whole Korea
NFI_CONIFER_HA <- 2324000     # official conifer(침엽수림) area (ha)

# ---- boundary ---------------------------------------------------------------------------------
# The totals below are computed inside the national boundary, because the covariate grid reaches
# past the border and out to sea. This clips for the statistics only: the delivered rasters are not
# clipped. Stage the S. Korea 시도 boundary (SGIS ctp_rvn) under data/raw/boundaries/. To use
# another file set BOUNDARY_SHP <- "<path/to.shp>" before sourcing.
# honor an explicit BOUNDARY_SHP override only if it points to a real file; otherwise auto-resolve
# (guards against a stale BOUNDARY_SHP left in the session by an earlier run)
if (!exists("BOUNDARY_SHP") || !is.character(BOUNDARY_SHP) || !file.exists(BOUNDARY_SHP)) {
  BOUNDARY_SHP <- here("data","raw","boundaries","ctp_rvn.shp")
  if (!file.exists(BOUNDARY_SHP)) {                          # search nested subfolders too (e.g. boundaries/ctp_rvn/ctp_rvn.shp)
    hits <- list.files(here("data","raw","boundaries"), pattern="\\.shp$", full.names=TRUE, recursive=TRUE)
    if (length(hits)) BOUNDARY_SHP <- hits[[1]]
  }
}
if (!file.exists(BOUNDARY_SHP))
  stop("Boundary shapefile not found. Stage the S. Korea boundary under data/raw/boundaries/ ",
       "(ctp_rvn.shp + .shx/.dbf/.prj sidecars), or set BOUNDARY_SHP <- \"<path/to.shp>\" before sourcing.")
korea <- vect(st_transform(st_union(st_make_valid(st_read(BOUNDARY_SHP, quiet=TRUE))), 5179))
clip  <- function(tif) terra::mask(terra::crop(rast(file.path(mapdir, tif)), korea), korea)
# map_file() builds a delivered map's name from FILE_PREDSET and FILE_VAR in R/paths.R. product_file() names
# the 30 m set 05a writes, and for any other map set under maps_<RES>/ the resolution token is swapped for RES.
map_file <- function(predset, tgt) {
  f <- product_file(predset, tgt, "map")
  if (RES != "30m") f <- sub("_30m[.]tif$", sprintf("_%s.tif", RES), f)
  f
}

# ---- per-layer sanity stats + area/stock totals --------------------------------------------
stat <- list(); area <- list()
for (predset in PREDSETS) {
  cf    <- clip(map_file(predset, "coniferfrac"))
  pa_ha <- prod(res(cf)) / 1e4                                              # pixel area (ha)
  fmask <- !is.na(cf)
  forest_ha   <- as.numeric(global(fmask, "sum", na.rm=TRUE)) * pa_ha       # forest = valid (mapped) pixels
  conifer_ba  <- as.numeric(global(cf, "sum", na.rm=TRUE)) * pa_ha          # BA-fraction-weighted conifer area
  conifer_dom <- as.numeric(global(cf > 0.5, "sum", na.rm=TRUE)) * pa_ha    # conifer-dominant (>0.5) area
  for (tgt in REG) {
    r <- if (tgt=="coniferfrac") cf else clip(map_file(predset, tgt))
    stat[[length(stat)+1]] <- data.table(predset, target=tgt,
      mean=round(as.numeric(global(r,"mean",na.rm=TRUE)),3), sd=round(as.numeric(global(r,"sd",na.rm=TRUE)),3),
      min=round(as.numeric(global(r,"min",na.rm=TRUE)),3),  max=round(as.numeric(global(r,"max",na.rm=TRUE)),3),
      pct_valid=round(100*as.numeric(global(!is.na(r),"mean",na.rm=TRUE)),1))
  }
  total_AGB_Mt <- as.numeric(global(clip(map_file(predset, "Bio_a_ha")), "sum", na.rm=TRUE)) * pa_ha / 1e6  # Mg/ha * ha = Mg -> Mt
  area[[length(area)+1]] <- data.table(predset, forest_kha=round(forest_ha/1e3,1),
    conifer_dom_kha=round(conifer_dom/1e3,1), conifer_ba_kha=round(conifer_ba/1e3,1),
    total_AGB_Mt=round(total_AGB_Mt,1),
    forest_vs_NFI=ifelse(is.na(NFI_FOREST_HA), NA, round(100*forest_ha/NFI_FOREST_HA,1)),
    coniferdom_vs_NFI=ifelse(is.na(NFI_CONIFER_HA), NA, round(100*conifer_dom/NFI_CONIFER_HA,1)))
}
stat <- rbindlist(stat); area <- rbindlist(area)
fwrite(stat, file.path(outd, sprintf("postproc_stats_%s.csv", RES)))
fwrite(area, file.path(outd, sprintf("postproc_area_%s.csv",  RES)))
cat("== per-layer sanity ==\n"); print(stat)
cat("\n== area / stock totals (kha = 1000 ha; AGB in Mt) ==\n"); print(area)
if (is.na(NFI_FOREST_HA)) cat("\nNOTE: set NFI_FOREST_HA / NFI_CONIFER_HA (KFS/NFI design-based) for the *_vs_NFI agreement columns.\n")
cat(sprintf("wrote postproc_{stats,area}_%s.csv to %s\n", RES, outd))

# ---- ensure every delivered GeoTIFF declares NoData = NaN -----------------------------------
# A float32 export can omit the NoData tag, and a GIS then renders the off-forest NaN background as
# data. terra reads NaN as NA either way, so the statistics above are unaffected. This sets NoData to
# NaN on each tif in maps_<RES> and uncertainty_<RES> that lacks a NoData tag, which skips the uint8
# dominant-species map (nodata 255) and makes reruns idempotent. Run it after any step that rewrites a
# raster. The GDAL bundled with sf writes a Cloud-Optimized GeoTIFF with overviews, or a tiled GTiff
# if the COG driver is unavailable.
ensure_nodata_nan <- function(dirs) {
  for (d in dirs) {
    if (!dir.exists(d)) next
    for (f in list.files(d, pattern="\\.tif$", full.names=TRUE)) {
      info <- tryCatch(sf::gdal_utils("info", f, quiet=TRUE), error=function(e) "")
      if (is.character(info) && any(grepl("NoData Value", info, fixed=TRUE))) next   # already tagged
      tmp <- paste0(f, ".ndtmp")                                                     # temp name (won't match *.tif glob)
      done <- tryCatch({                                                             # prefer a proper COG
        sf::gdal_utils("translate", f, tmp, options=c("-a_nodata","nan","-of","COG","-co","COMPRESS=LZW")); TRUE
      }, error=function(e) tryCatch({                                                # fall back to tiled GTiff
        sf::gdal_utils("translate", f, tmp,
                       options=c("-a_nodata","nan","-co","COMPRESS=LZW","-co","TILED=YES","-co","BIGTIFF=IF_SAFER")); TRUE
      }, error=function(e2){ message("  NoData set FAILED: ", basename(f), " - ", conditionMessage(e2)); FALSE }))
      if (isTRUE(done)) { file.rename(tmp, f); cat(sprintf("  NoData=NaN: %s\n", basename(f))) }
      else if (file.exists(tmp)) unlink(tmp)
    }
  }
}
if (ENSURE_NODATA) {
  cat("\n== ensuring NoData=NaN on delivered tifs ==\n")
  ensure_nodata_nan(c(mapdir, here("outputs", VERSION, sprintf("uncertainty_%s", RES))))
}
