# 04b_build_targets.R — build the four targets, whole-stand and canopy-only, on the frozen folds.
# ============================================================================================
# Per plot (신표본점번호 = Plot_ID): composition as a basal-area fraction, coniferous and evergreen;
# diversity as Shannon H', richness and evenness; the dominant species and its basal-area share; and
# stand structure — biomass, carbon, basal area per hectare, stems per hectare and quadratic mean
# diameter. Every basal-area quantity is per hectare — each tree's basal area times its expansion
# factor — and the unit is the 표본점.
#
# Composition and diversity are built for the whole stand and again for the canopy alone, under two
# definitions of canopy membership — the field crown class, and the largest trees making up 80 % of
# plot basal area — so that 07a can ask whether restricting to the canopy helps. Only the
# whole-stand versions are delivered.
#
# The script prints the distribution of every raw code it depends on (침활구분, 수관급,
# Leaf_phenology, and 수고 coverage) before using any of them. Read those tables: they are the check
# that the codings below match what the source means by them.
#
# Inputs : data/interim/nfi8/records.csv  + outputs/nfi8/cv_folds.csv (the Plot_ID universe)
# Output : outputs/nfi8/m2_targets.csv   (Plot_ID + <metric>_<stratum>; stratum in {full,crown,ba})
#   source(here::here("pipeline","04b_build_targets.R"))
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here::here("R","paths.R")); source(here::here("R","run_utils.R")); source(here::here("R","eval_utils.R"))

if (!exists("VERSION")) VERSION <- "nfi8"
BA_DOM_FRAC <- 0.80     # BA-dominance: trees making up the top 80% of plot basal area per hectare = canopy
outd <- here("outputs", VERSION); dir.create(outd, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("\n===== build targets [%s]  canopy rules={crown,ba}  ba_dom=%.2f =====\n",
            VERSION, BA_DOM_FRAC))
tm <- run_timer("m2_build_targets", version = VERSION, budget_min = 10)

# ---- 1. load + rename to stable internal names ---------------------------------------------
tm$stage("load")
sel <- c("신표본점번호","수종명","침활구분","흉고직경","수관급코드","수관급","수고","EF","Ba","Ca","Leaf_phenology")
rec <- fread(require_file(here("data","interim",VERSION,"records.csv")), encoding = "UTF-8",
             select = sel, colClasses = list(character = "신표본점번호"))
setnames(rec, c("신표본점번호","수종명","침활구분","흉고직경","수관급코드","수관급","수고","Leaf_phenology"),
              c("Plot_ID","sp","cb","dbh","crown_cd","crown_lbl","ht_meas","leaf"))
folds <- load_cv_folds(VERSION); universe <- unique(folds$Plot_ID)
cat(sprintf("records: %d trees / %d plots | cv_folds universe: %d plots\n",
            nrow(rec), uniqueN(rec$Plot_ID), length(universe)))

# ---- 2. DIAGNOSE raw codes (verify before trusting) ----------------------------------------
tm$stage("diagnose")
cat("\n-- 침활구분 (conifer/broadleaf) --\n");  print(rec[, .N, by = cb][order(-N)])
cat("\n-- Leaf_phenology --\n");               print(rec[, .N, by = leaf][order(-N)])
cat("\n-- 수관급 (crown class label) x 수관급코드 --\n"); print(rec[, .N, by = .(crown_cd, crown_lbl)][order(crown_cd)])
cat(sprintf("\n수고(measured) NA: %.1f%% | dbh NA/<=0: %.1f%%\n",
            100*mean(is.na(rec$ht_meas)), 100*mean(is.na(rec$dbh) | rec$dbh <= 0)))

# ---- 3. tree-level derived fields + canopy flags (robust label matching) -------------------
tm$stage("flags")
rec <- rec[!is.na(dbh) & dbh > 0 & !is.na(sp) & sp != ""]          # valid stems
rec <- rec[!grepl("고사|죽은|dead", crown_lbl, ignore.case = TRUE)] # drop dead (live-tree product)
# regional pine RACES are one species (Pinus densiflora) — merge so species/diversity targets are taxonomically correct
rec[sp %in% c("강원지방소나무","중부지방소나무"), sp := "소나무"]
rec[, ba := pi * (dbh/200)^2]                                      # basal area m^2 (dbh cm)
# Per-hectare weight. Trees of 6-30 cm are tallied on the 0.04 ha main circle and trees of 30 cm and
# over on the 0.08 ha large-tree circle, so a sum of unweighted basal areas gives each large tree
# twice the weight its per-hectare share warrants. Composition, diversity and the basal-area canopy
# rule use ba * EF, the same basis as the per-hectare structure below.
rec[, ba_w := ba * EF]
rec[, is_conifer   := grepl("침", cb)]                             # 침엽 -> conifer
rec[, is_evergreen := grepl("evergreen|상록", leaf, ignore.case = TRUE)]
# crown_cd labels are English, as the diagnostic printed above confirms:
# a = dominant, b = codominant, c = intermediate, d = suppressed, e = wolf.
# Canopy-forming = a + b, the trees that reach the main sunlit canopy; c and d sit below it.
CANOPY_CODES <- c("a", "b")
rec[, canopy_crown := crown_cd %in% CANOPY_CODES]
setorder(rec, Plot_ID, -ba)
rec[, canopy_ba := (cumsum(ba_w) - ba_w) / sum(ba_w) < BA_DOM_FRAC, by = Plot_ID]  # largest trees up to 80% of plot BA per ha
cat(sprintf("\ncanopy tree share - crown %.1f%% | ba %.1f%%  (of valid stems)\n",
            100*mean(rec$canopy_crown), 100*mean(rec$canopy_ba)))

# ---- 4. composition + diversity per stratum ------------------------------------------------
tm$stage("targets")
comp_div <- function(td, suf) {                                   # td = tree rows for one stratum
  if (!nrow(td)) return(data.table(Plot_ID = character(0)))
  sp_ba <- td[, .(ba = sum(ba_w)), by = .(Plot_ID, sp)][, tot := sum(ba), by = Plot_ID][, p := ba/tot]
  div   <- sp_ba[, .(shannon = -sum(p*log(p)), richness = .N), by = Plot_ID]
  div[, evenness := fifelse(richness > 1, shannon/log(richness), NA_real_)]
  dom   <- sp_ba[order(-ba), .(domsp = sp[1L], dom_ba_share = max(p)), by = Plot_ID]
  frac  <- td[, .(coniferfrac = sum(ba_w[is_conifer])/sum(ba_w),
                  evergreenfrac = sum(ba_w[is_evergreen])/sum(ba_w),
                  n_tree = .N, ba_tot = sum(ba_w)), by = Plot_ID]
  out <- Reduce(function(a,b) merge(a,b,by="Plot_ID",all=TRUE), list(frac, div, dom))
  setnames(out, setdiff(names(out),"Plot_ID"), paste0(setdiff(names(out),"Plot_ID"), "_", suf))
  out
}
strata <- list(full = rec, crown = rec[canopy_crown == TRUE], ba = rec[canopy_ba == TRUE])
comp <- Reduce(function(a,b) merge(a,b,by="Plot_ID",all=TRUE),
               Map(comp_div, strata, names(strata)))

# ---- 5. structure (FULL stand only) --------------------------------------------------------
stru <- rec[, .(Bio_a_ha = sum(Ba*EF), Carbon_a_ha = sum(Ca*EF), BAPH = sum(ba*EF),
                stems_ha = sum(EF), Dg = sqrt(sum(dbh^2*EF)/sum(EF))), by = Plot_ID]

# ---- 6. assemble onto the cv_folds universe + write ----------------------------------------
tm$stage("write")
tgt <- merge(data.table(Plot_ID = universe), merge(comp, stru, by = "Plot_ID", all = TRUE),
             by = "Plot_ID", all.x = TRUE)
fwrite(tgt, file.path(outd, "m2_targets.csv"))
# quick cross-check: the full-stand biomass vs the validated rd_points Bio_a_ha
rdp <- fread(here("data","interim",VERSION,"rd_points.csv"), select = c("Plot_ID","Bio_a_ha"),
             colClasses = list(character = "Plot_ID"))
chk <- merge(tgt[, .(Plot_ID, m2 = Bio_a_ha)], rdp[, .(Plot_ID, rd = Bio_a_ha)], by = "Plot_ID")
cat(sprintf("\nwrote %s/m2_targets.csv  (%d plots, %d target cols)\n", outd, nrow(tgt), ncol(tgt)-1L))
# Correlation alone cannot see a scale error: kg/ha against Mg/ha correlates at 1.0. The ratio of
# medians is what catches one, so both are printed and both have to pass.
cc    <- cor(chk$m2, chk$rd, use = "complete.obs")
ratio <- median(chk$m2, na.rm = TRUE) / median(chk$rd, na.rm = TRUE)
cat(sprintf("cross-check Bio_a_ha (targets against rd_points): cor=%.4f  median ratio=%.4f  over %d plots %s\n",
            cc, ratio, nrow(chk),
            if (cc > 0.99 && abs(ratio - 1) < 0.01) "[OK]" else "[CHECK recompute]"))
cat("\nheadline targets present (full + 3 canopy strata):\n")
print(grep("coniferfrac|shannon", names(tgt), value = TRUE))

tm$done(n = nrow(tgt), params = sprintf("strata=3 ba_dom=%.2f", BA_DOM_FRAC))
cat("\nVERIFY the printed 침활구분 / 수관급 / Leaf_phenology tables, then model these on the frozen folds.\n")
