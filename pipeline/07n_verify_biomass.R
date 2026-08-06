# 07n_verify_biomass.R — cross-check the biomass values against the inventory totals.
# ----------------------------------------------------------------------------
# The delivered per-ha biomass and carbon are 04b's aggregation of the curated per-tree values,
# Sum(Ba x EF) and Sum(Ca x EF), and the stages that build the maps deliberately keep those in place
# of the ones 01b carries over from the curated subplot_list. A units or expansion-factor error in
# that chain would make the direct route's weak biomass R2 an artefact of the data rather than a
# real result, and it would not show up anywhere else.
# This script re-derives both, independently of the pipeline, and cross-checks:
#
#   (A) per-tree coefficient chain (NIFoS):  Ba ?= Vol_eff*WD*BEF ;  Ca ?= Ba*CF ;
#                                            Bt ?= Ba*(1+R) ;        Ct ?= Bt*CF
#   (B) per-ha aggregation:  Bio_a_ha ?= Σ(Ba*EF) ;  Carbon_a_ha ?= Σ(Ca*EF)   (design-based estimator)
#       compared against subplot_list for a random sample of ~SAMPLE_N 표본점 (cycle 8 / NFI-8).
#       This is the same aggregation 04b performs, so it is also a check that the two agree.
#
# Small relative errors mean the units and expansion factors are right, so the weak biomass R² is a
# real property of the attribute — a state variable dominated by stand age and structure — and not a
# scaling artefact. Large errors mean the opposite and have to be resolved before anything is read
# from the biomass layer.
#
# Reads the curated RData only (light). Run in RStudio:
#   source(here::here("pipeline","07n_verify_biomass.R"))
# ----------------------------------------------------------------------------
CYC <- 8L; YRS <- 2021:2025          # NFI-8
SAMPLE_N <- 20L                      # 표본점 to cross-check on the per-ha estimator
TOL <- 0.01                          # 1% relative-error pass threshold

suppressPackageStartupMessages({ library(here); library(data.table) })
source(here::here("R", "paths.R"))

load(require_file(here::here("data", "raw", "nfi_curated", "NFI_Full.RData")))   # tree / subplot_list / ...
setDT(tree); setDT(subplot_list)
say  <- function(...) cat(sprintf(...), "\n")
relerr <- function(a, b) abs(a - b) / pmax(abs(b), .Machine$double.eps)
qsumm  <- function(x) sprintf("median %.3g | p95 %.3g | max %.3g", median(x, na.rm = TRUE),
                              quantile(x, .95, na.rm = TRUE), max(x, na.rm = TRUE))

tr <- tree[Nth == CYC & is.finite(Ba) & is.finite(EF)]
say("NFI-%d trees with finite Ba & EF: %s", CYC, format(nrow(tr), big.mark = ","))

# ---- (A) per-tree coefficient chain ----------------------------------------------------------------
chk <- tr[is.finite(Vol_eff) & is.finite(WD) & is.finite(BEF) & is.finite(CF)]
e_ba <- relerr(chk$Vol_eff * chk$WD * chk$BEF, chk$Ba)
e_ca <- relerr(chk$Ba * chk$CF,                chk$Ca)
e_bt <- if ("R"  %in% names(chk)) relerr(chk$Ba * (1 + chk$R), chk$Bt) else NA_real_
e_ct <- if ("Bt" %in% names(chk)) relerr(chk$Bt * chk$CF,      chk$Ct) else NA_real_
say("\n[A] per-tree coefficient chain (relative error, n=%s):", format(nrow(chk), big.mark = ","))
say("    Ba = Vol_eff*WD*BEF : %s", qsumm(e_ba))
say("    Ca = Ba*CF          : %s", qsumm(e_ca))
if (all(is.finite(e_bt))) say("    Bt = Ba*(1+R)       : %s", qsumm(e_bt))
if (all(is.finite(e_ct))) say("    Ct = Bt*CF          : %s", qsumm(e_ct))
chainA_ok <- median(e_ba, na.rm = TRUE) < TOL && median(e_ca, na.rm = TRUE) < TOL

# ---- (B) per-ha aggregation vs subplot_list --------------------------------------------------------
agg <- tr[, .(bio_recomp = sum(Ba * EF), car_recomp = sum(Ca * EF)), by = .(Plot_ID = Subplot_new)]
sl  <- unique(subplot_list[Year %in% YRS, .(Plot_ID = Subplot_new, BAPH, Bio_a_ha, Carbon_a_ha)],
              by = "Plot_ID")
cmp <- merge(agg, sl, by = "Plot_ID")
say("\n[B] per-ha cross-check: %s 표본점 present in BOTH tree(Nth==%d) and subplot_list(%d-%d)",
    format(nrow(cmp), big.mark = ","), CYC, min(YRS), max(YRS))

set.seed(7214)
samp <- cmp[sample(.N, min(SAMPLE_N, .N))]
samp[, `:=`(re_bio = relerr(bio_recomp, Bio_a_ha), re_car = relerr(car_recomp, Carbon_a_ha))]
print(samp[, .(Plot_ID, Bio_a_ha = round(Bio_a_ha, 1), bio_recomp = round(bio_recomp, 1),
               re_bio = round(re_bio, 4), Carbon_a_ha = round(Carbon_a_ha, 1),
               car_recomp = round(car_recomp, 1), re_car = round(re_car, 4))])
say("\n    Bio_a_ha    = Σ(Ba*EF) relerr (sample): %s", qsumm(samp$re_bio))
say("    Carbon_a_ha = Σ(Ca*EF) relerr (sample): %s", qsumm(samp$re_car))
# full-population summary too (cheap)
cmp[, `:=`(re_bio = relerr(bio_recomp, Bio_a_ha), re_car = relerr(car_recomp, Carbon_a_ha))]
say("    [all plots] Bio relerr: %s", qsumm(cmp$re_bio))
say("    [all plots] Carbon relerr: %s", qsumm(cmp$re_car))
aggB_ok <- median(cmp$re_bio, na.rm = TRUE) < TOL && median(cmp$re_car, na.rm = TRUE) < TOL

# ---- write result + verdict ------------------------------------------------------------------------
pass <- chainA_ok && aggB_ok
vdir <- here::here("outputs", "nfi8", "validation"); dir.create(vdir, showWarnings = FALSE, recursive = TRUE)
fwrite(samp, file.path(vdir, "verify_biomass.csv"), bom = TRUE)
out <- data.table(check = c("ba=vol*wd*bef_med", "ca=ba*cf_med", "bio_ha=sumBaEF_med",
                            "carbon_ha=sumCaEF_med", "PASS"),
                  value = c(round(median(e_ba, na.rm = TRUE), 5), round(median(e_ca, na.rm = TRUE), 5),
                            round(median(cmp$re_bio, na.rm = TRUE), 5), round(median(cmp$re_car, na.rm = TRUE), 5),
                            pass))
fwrite(out, file.path(vdir, "verify_biomass_summary.csv"), bom = TRUE)

cat("\n==== verify_biomass:", if (pass) "PASS ✅ (units/EF consistent → weak biomass R² is GENUINE)"
    else "FAIL ❌ (relative error > 1% → investigate units/EF; biomass weakness may be artifactual)", "====\n")
