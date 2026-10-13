# 06a_landcover_crosstab.R — forest x land-cover confusion table for the mask refinement.
#
# Reads the Level-1 land-cover raster from 02e together with the delivered direct maps, both at 30 m,
# and tabulates the land-cover class distribution across the mapped forest cells — those with a
# finite delivered value, which is the forest-type map's stocked extent — with the mean coniferous
# fraction and mean biomass of each class. Those two local signals separate the two ways a
# disagreement can arise:
#   - a confident error : a class that cannot carry trees, with low mean biomass -> an airport or an
#                         urban block that the forest-type map committed to forest.
#   - a likely omission : a class that can carry trees, with forest-like biomass -> orchards, edge
#                         regeneration, sparse cover that is genuinely wooded.
# It also reports the reverse direction, land cover that says forest where the maps do not, which is
# informational only — nothing is added back.
#
# The output table is the evidence for the per-class remove/preserve decision.
# OUTPUT: outputs/nfi8/m8_landcover_crosstab.csv  (+ console summary; + archived to the DB _derived/).
#   source(here::here("pipeline","06a_landcover_crosstab.R"))

suppressPackageStartupMessages({ library(here); library(terra); library(data.table) })
source(here("R", "paths.R"))

L1_TIF  <- here("data", "interim", "landcover_mcee_30m", "landcover_L1_30m.tif")
CONIF   <- here("outputs", "nfi8", "m8_maps_30m", "ftm_conifer_fraction_30m.tif")   # forest extent + signal
# 06b filters on the land-cover class alone; the biomass map is tabulated here for the summary only.
BIO     <- here("outputs", "nfi8", "m8_maps_30m", "ftm_biomass_30m.tif")      # mean biomass per class
OUT_CSV <- here("outputs", "nfi8", "m8_landcover_crosstab.csv")
for (f in c(L1_TIF, CONIF, BIO)) require_file(f)

L1_NAME <- c("100"="built(시가화건조)","200"="agri(농업)","300"="forest(산림)","400"="grass(초지)",
             "500"="wetland(습지)","600"="barren(나지)","700"="water(수역)")
CELL_HA <- 30 * 30 / 1e4   # 0.09 ha per 30 m cell

L1  <- terra::rast(L1_TIF)
cf  <- terra::rast(CONIF)
bio <- terra::rast(BIO)

# The land-cover raster is rasterised from the source polygons on its own 30 m lattice, which is not
# cell-aligned with the delivered grid, so crop/extend cannot bring the two together. Regrid it by
# nearest neighbour onto the delivered grid, exactly as 06b does for the filter itself; a sub-cell
# shift does not matter for a class tabulation.
if (!terra::compareGeom(L1, cf, stopOnError = FALSE, messages = FALSE)) {
  L1 <- terra::resample(L1, cf, method = "near")
  stopifnot(terra::compareGeom(L1, cf, stopOnError = FALSE, messages = FALSE))
}

forest <- !is.na(cf)                                   # our mapped forest extent

## ---- distribution of OUR forest cells across MCEE L1 classes ------------------------------------
L1f  <- terra::mask(L1, forest, maskvalue = FALSE)     # L1 only where we mapped forest
fr   <- as.data.table(terra::freq(L1f))
setnames(fr, "value", "L1"); fr <- fr[!is.na(L1), .(L1, cells = count)]
zcf  <- as.data.table(terra::zonal(cf,  L1f, "mean", na.rm = TRUE)); setnames(zcf,  1:2, c("L1","mean_conif"))
zbio <- as.data.table(terra::zonal(bio, L1f, "mean", na.rm = TRUE)); setnames(zbio, 1:2, c("L1","mean_bio"))
tab  <- Reduce(function(a, b) merge(a, b, by = "L1", all.x = TRUE), list(fr, zcf, zbio))

totf <- sum(tab$cells)
tab[, `:=`(class      = L1_NAME[as.character(L1)],
           area_ha    = round(cells * CELL_HA, 0),
           pct_forest = round(100 * cells / totf, 2),
           mean_conif = round(mean_conif, 3),
           mean_bio   = round(mean_bio, 1))]
setorder(tab, -cells)
tab <- tab[, .(L1, class, cells, area_ha, pct_forest, mean_conif, mean_bio)]

## ---- reverse: MCEE forest (300) we did NOT map (omission; informational) -------------------------
L1_nf <- terra::mask(L1, forest, maskvalue = TRUE)     # L1 where our map is NA
frn   <- as.data.table(terra::freq(L1_nf)); setnames(frn, "value", "L1")
omit300 <- frn[L1 == 300, count]; if (length(omit300) == 0) omit300 <- 0L

fwrite(tab, OUT_CSV, bom = TRUE)
ARCHIVE <- here("data", "raw", "landcover_mcee_national", "_derived")
if (dir.exists(ARCHIVE)) file.copy(OUT_CSV, file.path(ARCHIVE, basename(OUT_CSV)), overwrite = TRUE)

## ---- console summary ---------------------------------------------------------------------------
cat(sprintf("\nOur forest cells: %s  (%s ha)\n",
            format(totf, big.mark = ","), format(round(totf * CELL_HA), big.mark = ",")))
cat(sprintf("  agreement  MCEE=forest(300): %.1f%%\n", 100 * tab[L1 == 300, cells] / totf))
cat(sprintf("  disagreement (non-300)    : %.1f%%\n", 100 * (totf - tab[L1 == 300, cells]) / totf))
print(tab)
cat(sprintf("\nReverse (informational): MCEE forest(300) OUTSIDE our map = %s cells (%s ha)\n",
            format(omit300, big.mark = ","), format(round(omit300 * CELL_HA), big.mark = ",")))
cat("\nThe decision this evidence led to, applied in 06b_apply_landcover_filter.R:\n",
    "  removed  built(100), barren(600), water(700) — surfaces that cannot carry trees\n",
    "  kept     agri(200), grass(400), wetland(500) — disagreement there is often real tree cover\n",
    "Wrote:", OUT_CSV, "\n")
