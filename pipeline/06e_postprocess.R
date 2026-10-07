# 06e_postprocess.R — summary statistics and area totals for the finished maps.
# ============================================================================================
# Post-processes the wall-to-wall maps (outputs/nfi8/m8_maps_<RES>/): clip to the national boundary, per-layer
# sanity stats, and area/stock totals with an NFI design-based AREA-AGREEMENT check (a 2nd, leakage-free QA beyond
# the test-plot R²: does mapped conifer/forest area match the official NFI estimate?).
# RES picks which map set to summarise. The delivered set is 30 m.
#   RES <- "30m"; source(here::here("pipeline","06e_postprocess.R"))
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(terra); library(data.table); library(sf) })
source(here::here("R", "paths.R"))
VERSION <- "nfi8"
if (!exists("RES"))       RES       <- "30m"
if (!exists("WRITE_CLIP")) WRITE_CLIP <- FALSE     # TRUE -> also write boundary-clipped tifs (clean product) to m8_maps_<RES>_clip/
if (!exists("ENSURE_NODATA")) ENSURE_NODATA <- TRUE # set NoData=NaN on the delivered tifs (a float export can omit the tag)
mapdir <- here("outputs", VERSION, sprintf("m8_maps_%s", RES)); stopifnot(dir.exists(mapdir))
outd   <- here("outputs", VERSION)
TIERS <- "full"; REG <- c("coniferfrac","shannon","Bio_a_ha")     # "full" = the delivered predictor set (R/paths.R::FILE_TIER)

# ---- NFI design-based reference (KFS 산림기본통계, 2020 basis; whole Korea) --------------------
# Design-based (sample-expansion) area estimates, the leakage-free reference for the agreement check.
# Source = Korea Forest Service 2021 (산림기본통계) / NIFoS 2020.
# SWAP for the NFI-8 yearbook values once published. NB the mapped "conifer-dominant (cf>0.5)" area is
# only APPROXIMATELY the NFI 침엽수림 class (majority-conifer vs stand-type definition) -> report as
# broad area agreement, not an exact match.
NFI_FOREST_HA  <- 6298000     # official total forest area (ha), whole Korea (2020 basis)
NFI_CONIFER_HA <- 2324000     # official conifer(침엽수림) area (ha) (2020 basis); cf. conifer-dominant mapped area

# ---- boundary ---------------------------------------------------------------------------------
# The totals below are computed inside the national boundary, because the covariate grid reaches
# past the border and out to sea. This clips for the statistics only: the delivered rasters are not
# clipped, and WRITE_CLIP is FALSE. Stage the S. Korea 시도 boundary (SGIS ctp_rvn) under
# data/raw/boundaries/. To use another file set BOUNDARY_SHP <- "<path/to.shp>" before sourcing.
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
# One place builds a delivered map's name, so this stage and R/paths.R cannot drift apart.
map_file <- function(tier, tgt) sprintf("%s_%s_%s.tif", unname(FILE_TIER[tier]), unname(FILE_VAR[tgt]), RES)
if (WRITE_CLIP) { clipdir <- here("outputs",VERSION,sprintf("m8_maps_%s_clip",RES)); dir.create(clipdir, showWarnings=FALSE) }

# ---- per-layer sanity stats + area/stock totals --------------------------------------------
stat <- list(); area <- list()
for (tier in TIERS) {
  cf    <- clip(map_file(tier, "coniferfrac"))
  pa_ha <- prod(res(cf)) / 1e4                                              # pixel area (ha)
  fmask <- !is.na(cf)
  forest_ha   <- as.numeric(global(fmask, "sum", na.rm=TRUE)) * pa_ha       # forest = valid (mapped) pixels
  conifer_ba  <- as.numeric(global(cf, "sum", na.rm=TRUE)) * pa_ha          # BA-fraction-weighted conifer area
  conifer_dom <- as.numeric(global(cf > 0.5, "sum", na.rm=TRUE)) * pa_ha    # conifer-dominant (>0.5) area
  for (tgt in REG) {
    r <- if (tgt=="coniferfrac") cf else clip(map_file(tier, tgt))
    if (WRITE_CLIP) writeRaster(r, file.path(clipdir, sub("\\.tif$", "_clip.tif", map_file(tier, tgt))), overwrite=TRUE, gdal="COMPRESS=LZW")
    stat[[length(stat)+1]] <- data.table(tier, target=tgt,
      mean=round(as.numeric(global(r,"mean",na.rm=TRUE)),3), sd=round(as.numeric(global(r,"sd",na.rm=TRUE)),3),
      min=round(as.numeric(global(r,"min",na.rm=TRUE)),3),  max=round(as.numeric(global(r,"max",na.rm=TRUE)),3),
      pct_valid=round(100*as.numeric(global(!is.na(r),"mean",na.rm=TRUE)),1))
  }
  total_AGB_Mt <- as.numeric(global(clip(map_file(tier, "Bio_a_ha")), "sum", na.rm=TRUE)) * pa_ha / 1e6  # Mg/ha * ha = Mg -> Mt
  area[[length(area)+1]] <- data.table(tier, forest_kha=round(forest_ha/1e3,1),
    conifer_dom_kha=round(conifer_dom/1e3,1), conifer_ba_kha=round(conifer_ba/1e3,1),
    total_AGB_Mt=round(total_AGB_Mt,1),
    forest_vs_NFI=ifelse(is.na(NFI_FOREST_HA), NA, round(100*forest_ha/NFI_FOREST_HA,1)),
    coniferdom_vs_NFI=ifelse(is.na(NFI_CONIFER_HA), NA, round(100*conifer_dom/NFI_CONIFER_HA,1)))
}
stat <- rbindlist(stat); area <- rbindlist(area)
fwrite(stat, file.path(outd, sprintf("m8_postproc_stats_%s.csv", RES)))
fwrite(area, file.path(outd, sprintf("m8_postproc_area_%s.csv",  RES)))
cat("== per-layer sanity ==\n"); print(stat)
cat("\n== area / stock totals (kha = 1000 ha; AGB in Mt) ==\n"); print(area)
if (is.na(NFI_FOREST_HA)) cat("\nNOTE: set NFI_FOREST_HA / NFI_CONIFER_HA (KFS/NFI design-based) for the *_vs_NFI agreement columns.\n")
cat(sprintf("wrote m8_postproc_{stats,area}_%s.csv to %s\n", RES, outd))

# ---- ensure every delivered GeoTIFF declares NoData = NaN -----------------------------------
# A float32 export can omit the NoData tag, in which case a GIS renders the off-forest NaN background as if
# it were data (terra reads NaN as NA regardless, so the stats above are unaffected — this is purely a
# delivered-file fix for downstream GIS). Bakes NoData=NaN into each float tif in m8_maps_<RES> and
# m8_uncertainty_<RES>. Idempotent: files that already declare NoData are skipped, which includes the
# uint8 dominant-species map (nodata 255) and the int8 Shannon display raster (nodata -128). Run it
# after any step that rewrites a raster. Prefers a Cloud-Optimized GeoTIFF (nodata + overviews, the
# delivery format); falls back to a plain tiled GTiff if the COG driver is unavailable. Uses the GDAL
# already bundled with sf.
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
  ensure_nodata_nan(c(mapdir, here("outputs", VERSION, sprintf("m8_uncertainty_%s", RES))))
}
