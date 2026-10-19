# 01a_extract_nfi_records.R — extract the inventory tree records this dataset is built from.
# ------------------------------------------------------------------------------------------------
# The source is a curated multi-cycle integration of the Korean National Forest Inventory
# (제5–8차 국가산림자원조사, cycles 5 to 8, 2006-2025). Relative to the raw inventory tables it
# corrects species names, recomputes per-tree biomass and carbon with the official national
# coefficients, corrects the expansion factors, resolves the old-to-new plot-number crosswalk and
# removes duplicates. Everything downstream reads the result rather than the raw tables.
#
# This script selects one cycle and writes it as a flat table with the inventory's own Korean
# column names, so the columns can be checked directly against the field manual.
#
# In : data/raw/nfi_curated/NFI_Full.RData
# Out: data/interim/nfi8/records.csv
#
#   source(here::here("pipeline", "01a_extract_nfi_records.R"))   # run in RStudio
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("R", "paths.R"))

CYCLE <- 8L                    # 8 = 2021-2025

rdata <- require_file(here("data", "raw", "nfi_curated", "NFI_Full.RData"))
load(rdata); stopifnot(exists("tree")); setDT(tree)

build_records <- function(d) {
  data.table(
    # ---- measured fields, under the inventory's own column names -------------------------------
    집락번호                = d$Cluster,
    표본점번호              = d$Subplot,                 # original subplot id; the unified id is added below
    조사차기                = d$Nth,
    번호                    = d$Number,
    학명번호                = d$C_Species,
    수종명                  = d$Species_k,               # curated
    침활구분                = fifelse(d$Species_group == "conifer",  "침엽수",
                              fifelse(d$Species_group == "broadleaf", "활엽수", NA_character_)),
    흉고직경                = d$DBH,
    수관급코드              = d$C_Crown_class,           # letter code, language-neutral
    수관급                  = d$Crown_class,             # curated: English labels, not the Korean originals
    대경목조사원내존재여부  = d$Large_tree,              # raw flag; the corrected form is Large_tree_true below
    지하고                  = d$HCB,
    수고                    = d$HT,
    수령                    = d$Age,                     # sparse; measured on sample trees only
    표준목간재적            = d$Vol,                     # sparse; measured
    추정수고                = d$HT_est,
    추정간재적              = d$Vol_est,                 # estimated stem volume, populated for every tree
    # ---- values the curation adds ---------------------------------------------------------------
    신표본점번호            = d$Subplot_new,             # unified plot id across cycles 5 to 8
    Vol_eff = d$Vol_eff, EF = d$EF, Large_tree_true = d$Large_tree_true,
    WD = d$WD, BEF = d$BEF, CF = d$CF, R = d$R,
    Ba = d$Ba, Bt = d$Bt, Ca = d$Ca, Ct = d$Ct,
    Family = d$Family, Genus = d$Genus, Species = d$Species,
    Leaf_phenology = d$Leaf_phenology, Species_group = d$Species_group
  )
}

out <- here("data", "interim", VERSION)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
fwrite(build_records(tree[Nth == CYCLE]), file.path(out, "records.csv"), bom = TRUE)

cat(sprintf("cycle %d: %d tree records -> data/interim/%s/records.csv\n",
            CYCLE, tree[Nth == CYCLE, .N], VERSION))
print(tree[, .N, by = Nth][order(Nth)])
cat("\n  Check 수관급 and 침활구분 against the field manual before using them as categories.\n")
