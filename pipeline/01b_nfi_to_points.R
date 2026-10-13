# 01b_nfi_to_points.R — build the per-plot reference table from the inventory tree records.
# ------------------------------------------------------------------------------------------------
# One row per plot, carrying everything later stages read off a plot: its coordinates, its forest
# type, diameter and density class, its dominant species, and its per-hectare biomass, carbon and
# basal area. Covariates are extracted at the plot coordinates later, in 03d.
#
# The unit is the 표본점 (subplot), not the 집락 (cluster). Each 표본점 is classified independently
# in the field, and the four 표본점 of a cluster sit 50 m apart, further than one 30 m cell.
#
# The plot coordinate comes from the plot register, except where the stand table, which records the
# position a second time, disagrees with it; the coordinate step below says how that is decided.
#
# In : data/raw/nfi_curated/NFI_Full.RData   (tree, stand, general, subplot_list)
#      data/species_forest_type_lookup.csv
# Out: data/interim/nfi8/plot_reference.csv
#      data/interim/nfi8/coordinate_corrections.csv   (the plots whose coordinate was replaced)
# An existing plot_reference.csv is moved to data/interim/nfi8/_superseded_<date>/ first, never
# overwritten. Both files hold exact plot coordinates and stay under data/interim/, which is not
# distributed.
#
#   source(here::here("pipeline", "01b_nfi_to_points.R"))   # run in RStudio
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(data.table); library(dplyr); library(purrr) })
source(here("R", "paths.R"))

CYCLE <- 8L                    # inventory cycle; 8 = 2021-2025, the cycle this dataset is built on

# Inventory codings for diameter and crown-density class, as recorded in the stand table.
#   diameter: 치수 = 3 cm, 소경목 = 12, 중경목 = 24, 대경목 = 36
#   density : 0-25 % and 26-50 % -> 25, 51-75 % -> 60, 76-100 % -> 85
recode_vdc <- function(dbh_class) {
  dplyr::recode(dbh_class,
                "치수" = 3L, "소경목" = 12L, "중경목" = 24L, "대경목" = 36L,
                .default = NA_integer_)
}
recode_vdn <- function(crown_density) {
  dplyr::case_when(
    crown_density %in% c("0-25%", "26-50%") ~ 25L,
    crown_density == "51-75%"               ~ 60L,
    crown_density == "76-100%"              ~ 85L,
    TRUE                                    ~ NA_integer_
  )
}

load(here("data", "raw", "nfi_curated", "NFI_Full.RData"))     # tree / stand / general / subplot_list
setDT(tree); setDT(stand); setDT(general); setDT(subplot_list)

## ---- species -> forest-type group -------------------------------------------------------------
## Maps the inventory species onto the forest-type groups of the forest-type map, so that plot
## composition and map composition are expressed in one vocabulary. Joined on Latin genus and
## species, falling back to the Korean name, and finally to "other broadleaf" (50).
species_ft  <- fread(here("data", "species_forest_type_lookup.csv"), encoding = "UTF-8")
ft_by_latin <- unique(species_ft[, .(Genus, Species, FT)], by = c("Genus", "Species"))
ft_by_ko    <- unique(species_ft[, .(Species_ko = Name_Korean, FT_ko = FT)])[!duplicated(Species_ko)]

## ---- forest-type dominance rule ---------------------------------------------------------------
## Evaluated in forest-type-group space rather than by naming species, so a species missing from
## the lookup cannot silently break a rule. Rules are tried in rank order; the first whose member
## groups exceed their threshold wins. A plot that clears none is mixed conifer-broadleaf (60).
# THR_SP  a single group must hold this share of the plot's importance value to name the plot.
# THR_GP  a set of groups (any conifer, any oak, any broadleaf) needs this higher share instead.
# DELTA   lowers both thresholds for the conifer rules and chestnut, which are tried first.
THR_SP <- 0.55; THR_GP <- 0.65; DELTA <- 0.05
classify_ft <- function(plot_ft_iv) {           # data.table(Plot_ID, FT, IV_sum)
  SP_LOW <- max(THR_SP - DELTA, 0); GP_LOW <- max(THR_GP - DELTA, 0)
  rules <- list(
    list(1L, 11L, 11L,   SP_LOW), list(2L, 12L, 12L, SP_LOW), list(3L, 13L, 13L, SP_LOW),
    list(4L, 14L, 14L,   SP_LOW), list(5L, 15L, 15L, SP_LOW), list(6L, 16L, 16L, SP_LOW),
    list(7L, 10L, 10:16, GP_LOW),                 # any conifer
    list(8L, 30L, 30L,   SP_LOW),                 # chestnut
    list(9L, 40L, 40L,   THR_GP),                 # evergreen broadleaf
    list(10L, 23L, 23L,  THR_SP), list(11L, 21L, 21L, THR_SP), list(12L, 22L, 22L, THR_SP),
    list(13L, 20L, 20:23, THR_GP),                # any oak
    list(14L, 50L, 20:50, THR_GP)                 # any broadleaf
  )
  hits <- data.table::rbindlist(lapply(rules, function(r) {
    h <- plot_ft_iv[FT %in% r[[3]], .(IV_sum = sum(IV_sum)), by = Plot_ID][IV_sum >= r[[4]]]
    if (nrow(h)) h[, `:=`(V_VG = r[[2]], pri = r[[1]])] else NULL
  }))
  res <- data.table(Plot_ID = character(0), V_VG = integer(0))
  if (nrow(hits)) res <- hits[order(Plot_ID, pri, -IV_sum), .(V_VG = V_VG[1]), by = Plot_ID]
  res <- res[unique(plot_ft_iv[, .(Plot_ID)]), on = "Plot_ID"]   # right-join so no plot is lost
  res[is.na(V_VG), V_VG := 60L]
  res[]
}

## ---- per-plot builder --------------------------------------------------------------------------
build_ref <- function(cyc) {
  yrs <- switch(as.character(cyc), "5" = 2006:2010, "6" = 2011:2015, "7" = 2016:2020, "8" = 2021:2025)

  # Importance value per forest-type group, from canopy trees only, weighted by the expansion
  # factor. The expansion factor is how many stems per hectare a tree represents, and weighting by
  # it corrects the over-sampling of large trees in the large-tree survey circle.
  tr <- tree[Nth == cyc & is.finite(DBH) & is.finite(EF) & EF > 0,
             .(Plot_ID = Subplot_new, Species_ko = Species_k, Genus, Species, Crown_class,
               EF = as.numeric(EF), BA = pi * (DBH / 2)^2)]
  tr <- merge(tr, ft_by_latin, by = c("Genus", "Species"), all.x = TRUE)
  if (anyNA(tr$FT)) {
    tr <- merge(tr, ft_by_ko, by = "Species_ko", all.x = TRUE)
    tr[is.na(FT), FT := FT_ko][, FT_ko := NULL]
  }
  unmatched <- sort(unique(tr[is.na(FT), Species_ko]))
  tr[is.na(FT), FT := 50L]
  trd <- tr[Crown_class %in% c("dominant", "codominant")]                # canopy trees
  iv  <- trd[, .(nef = sum(EF), baef = sum(BA * EF)), by = .(Plot_ID, FT)]
  iv[, `:=`(rel_n = nef / sum(nef), rel_BA = baef / sum(baef)), by = Plot_ID]
  iv[, IV := 0.5 * (rel_n + rel_BA)]                                     # mean of relative stems and basal area
  plot_ft_iv <- iv[, .(Plot_ID, FT, IV_sum = IV)]
  vg <- classify_ft(plot_ft_iv)

  # Diameter and density class as the field crew recorded them for this plot.
  st <- stand[Nth == cyc, .(Plot_ID = Subplot_new, DBH_class, Crown_density)]
  st[, `:=`(V_DC = recode_vdc(DBH_class), V_DN = recode_vdn(Crown_density))]
  st <- unique(st[, .(Plot_ID, V_DC, V_DN)])

  ge  <- unique(general[Nth == cyc, .(Plot_ID = Subplot_new, UTMK_X, UTMK_Y, F_or_NonF)])

  # ---- plot coordinates --------------------------------------------------------------------------
  # The register (general) and the stand table (stand) each record the plot position. From 2022 the
  # stand table's coordinates are in the register's coordinate system; the 2021 ones are on another
  # datum and are not compared. Where the two differ by more than 50 m, the register position is the
  # one found to be wrong (it disagrees with the plot's own survey records), so those plots take the
  # stand coordinate. A disagreement of 2 km or more stops the script, because then neither source
  # can be trusted for that plot.
  sc  <- unique(stand[Nth == cyc & Year >= 2022,
                      .(Plot_ID = Subplot_new, sX = as.numeric(Coord_E), sY = as.numeric(Coord_N))])
  sc  <- sc[is.finite(sX) & is.finite(sY)]
  stopifnot(!anyDuplicated(sc$Plot_ID))
  fix <- merge(ge[, .(Plot_ID, gX = UTMK_X, gY = UTMK_Y)], sc, by = "Plot_ID")
  fix[, d_gs := sqrt((sX - gX)^2 + (sY - gY)^2)]
  fix <- fix[d_gs > 50]
  if (any(fix$d_gs >= 2000))
    stop("register and stand coordinates are 2 km or more apart for: ",
         paste(fix[d_gs >= 2000, Plot_ID], collapse = ", "), call. = FALSE)
  ge[fix, on = "Plot_ID", `:=`(UTMK_X = i.sX, UTMK_Y = i.sY)]
  if (nrow(fix))
    cat(sprintf("coordinates: %d plots take the stand-table position (median shift %.0f m, max %.0f m)\n",
                nrow(fix), median(fix$d_gs), max(fix$d_gs)))

  sl  <- unique(subplot_list[Year %in% yrs, .(Plot_ID = Subplot_new, BAPH, Bio_a_ha, Carbon_a_ha)],
                by = "Plot_ID")                                          # this cycle's years only
  # Dominant species by expansion-factor-weighted basal area, describing the plot as surveyed. The
  # mapped attribute is a different quantity, rebuilt in 04b on the same per-hectare basis: live
  # trees only, and the two regional pine races merged into Pinus densiflora.
  dom <- tr[, .(baef = sum(BA * EF)), by = .(Plot_ID, Species_ko)][
            order(Plot_ID, -baef), .(dom_species = Species_ko[1], dom_ba_share = baef[1] / sum(baef)), by = Plot_ID]

  ref <- Reduce(function(a, b) merge(a, b, by = "Plot_ID", all.x = TRUE), list(ge, st, vg, dom, sl))
  # Forested plots that the tree records place in a forest-type group. V_VG is NA only where a plot
  # contributed no canopy tree. It selects the 17,099-plot analysis set, so changing it changes
  # every number that follows.
  ref <- ref[F_or_NonF == 1 & !is.na(V_VG)]
  list(ref = ref, unmatched = unmatched, fix = fix)
}

FT_NAME <- c("10"="기타침엽수","11"="곰솔","12"="낙엽송","13"="리기다소나무","14"="소나무","15"="잣나무",
             "16"="편백나무","20"="기타참나무","21"="굴참나무","22"="상수리나무","23"="신갈나무","30"="밤나무",
             "40"="상록활엽수","50"="기타활엽수","60"="침활혼효림")
distr <- function(dt, col, nm = NULL) {                      # count and proportion, in code order
  t <- dt[, .N, by = col][order(get(col))]
  t[, prop := round(N / sum(N), 4)]
  if (!is.null(nm)) t[, label := nm[as.character(get(col))]]
  t
}

r <- build_ref(CYCLE)
o <- here("data", "interim", VERSION); dir.create(o, showWarnings = FALSE, recursive = TRUE)
if (file.exists(file.path(o, "plot_reference.csv"))) {           # keep the previous table, never overwrite it
  a <- file.path(o, sprintf("_superseded_%s", format(Sys.Date(), "%Y%m%d"))); dir.create(a, showWarnings = FALSE)
  stopifnot(file.rename(file.path(o, "plot_reference.csv"), file.path(a, "plot_reference.csv")))
}
fwrite(r$ref, file.path(o, "plot_reference.csv"), bom = TRUE)
# The corrected plots, with both positions and the distance behind the decision (interim only: it
# carries plot coordinates and is never published).
fwrite(r$fix[, .(Plot_ID, UTMK_X_general = gX, UTMK_Y_general = gY, UTMK_X = sX, UTMK_Y = sY,
                 d_general_stand_m = round(d_gs, 1), in_reference = Plot_ID %in% r$ref$Plot_ID)],
       file.path(o, "coordinate_corrections.csv"), bom = TRUE)
cat(sprintf("\n%s -> plot_reference.csv  (%d forested plots); coordinate_corrections.csv (%d plots)\n",
            VERSION, nrow(r$ref), nrow(r$fix)))

# Print the class distributions. They are the first check that the codings above were read
# correctly, so look at them before running anything downstream.
cat("\nforest-type group (V_VG):\n");   print(distr(r$ref, "V_VG", FT_NAME))
cat("\ndiameter class, cm (V_DC):\n");  print(distr(r$ref, "V_DC"))
cat("\ncrown density, % (V_DN):\n");    print(distr(r$ref, "V_DN"))
if (length(r$unmatched))
  cat("\n! species absent from species_forest_type_lookup.csv, mapped to 50:\n  ",
      paste(head(r$unmatched, 40), collapse = ", "), "\n")
