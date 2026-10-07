# 07m_verify_index_alignment.R — confirm the decoded attribute rasters really carry the matched
# plot's values.
# ------------------------------------------------------------------------------------------------
# Stage 05c builds each attribute raster with terra::subst(index, from = plot_idx, to = attribute).
# If the plot table and the index raster ever disagreed on what plot_idx means, that substitution
# would still succeed and produce a plausible-looking but wrong map. This checks it did not.
#
#   (1) invariant: the plot table has plot_idx == 1..N, in order.
#   (2) round-trip: for a random sample of cells, the index value v must decode through the plot
#       table to exactly the value the attribute raster holds at that cell. Mismatches must be 0.
#   (3) coverage: the index raster must have a non-zero, plausible number of forest cells.
#
# Run it after 05c, on the delivered index raster.
#   source(here::here("pipeline", "07m_verify_index_alignment.R"))
# ------------------------------------------------------------------------------------------------
N_SAMPLE <- 50000L    # random non-NA cells to round-trip

suppressPackageStartupMessages({ library(here); library(terra); library(data.table) })
source(here::here("R", "paths.R"))
.tsuf  <- tier_suffix("full")                       # "_ftm", the suffix 05b wrote; see R/paths.R
outd   <- here("outputs", VERSION)
attr_d <- file.path(outd, sprintf("attributes%s", .tsuf))
ok_all <- TRUE
say <- function(...) cat(sprintf(...), "\n")
say("verifying imputation%s.tif against attributes%s/", .tsuf, .tsuf)

# ---- (1) plot-table invariant -----------------------------------------------------------------
ref <- fread(require_file(file.path(outd, "models", sprintf("train_reference%s.csv", .tsuf))), encoding = "UTF-8")
setorder(ref, plot_idx)
N <- nrow(ref)
inv_ok <- isTRUE(all(ref$plot_idx == seq_len(N)))
ok_all <- ok_all && inv_ok
say("[1] plot table: N=%d | plot_idx == 1..N : %s", N, inv_ok)

# ---- load the index and two attribute rasters (same grid, both derived from the index) -------------
im   <- terra::rast(require_file(file.path(outd, "index", sprintf("imputation%s.tif", .tsuf)))); names(im)   <- "idx"
vvg  <- terra::rast(require_file(file.path(attr_d, "spec30m_V_VG.tif")));        names(vvg)  <- "vvg"
doms <- terra::rast(require_file(file.path(attr_d, "spec30m_dom_species.tif"))); names(doms) <- "doms"
leg  <- fread(require_file(file.path(attr_d, "spec30m_dom_species_legend.csv")), encoding = "UTF-8")

# dom_species integer code as 05c built it: match(ref$dom_species, sort(unique(ref$dom_species)))
sp_levels <- sort(unique(ref$dom_species))
ref[, dom_code := match(dom_species, sp_levels)]
leg_ok <- isTRUE(all(leg$code == seq_along(sp_levels)) && all(leg$species == sp_levels))
ok_all <- ok_all && leg_ok
say("[2] dom_species legend matches sort(unique(ref$dom_species)) : %s  (%d species)", leg_ok, length(sp_levels))

# ---- (2) round-trip on a random non-NA cell sample -------------------------------------------------
set.seed(7214)
stk  <- c(im, vvg, doms)                                           # identical grid -> per-cell aligned
samp <- terra::spatSample(stk, size = N_SAMPLE, method = "random", na.rm = TRUE, warn = FALSE)
samp <- as.data.table(samp)[!is.na(idx)]
idx  <- as.integer(samp$idx)
say("    sampled %d non-NA index cells (idx range %d..%d, expect 1..%d)", nrow(samp), min(idx), max(idx), N)

in_range <- all(idx >= 1L & idx <= N)
exp_vvg  <- as.numeric(ref$V_VG[idx])
exp_code <- ref$dom_code[idx]
mis_vvg  <- sum(round(exp_vvg) != round(samp$vvg),  na.rm = TRUE)
mis_code <- sum(exp_code      != round(samp$doms), na.rm = TRUE)
rt_ok <- in_range && mis_vvg == 0L && mis_code == 0L
ok_all <- ok_all && rt_ok
say("[3] round-trip: idx in 1..N : %s | V_VG mismatches : %d | dom_species mismatches : %d",
    in_range, mis_vvg, mis_code)

# ---- (3) coverage (non-NA fraction of the index grid) ----------------------------------------------
ncell_tot <- terra::ncell(im)
nna <- terra::global(!is.na(im), "sum", na.rm = TRUE)[1, 1]        # chunked count, out-of-core safe
say("[4] coverage: non-NA index cells = %s of %s (%.1f%% of grid)",
    format(nna, big.mark = ","), format(ncell_tot, big.mark = ","), 100 * nna / ncell_tot)
cov_ok <- nna > 0
ok_all <- ok_all && cov_ok

# ---- write a tiny result file + verdict ------------------------------------------------------------
out <- data.table(check = c("train_ref_invariant", "legend_match", "roundtrip_vvg_mismatch",
                            "roundtrip_domsp_mismatch", "idx_in_range", "noncov_cells_nonNA", "PASS"),
                  value = c(inv_ok, leg_ok, mis_vvg, mis_code, in_range, nna, ok_all))
vdir <- file.path(outd, "validation"); dir.create(vdir, showWarnings = FALSE, recursive = TRUE)
fwrite(out, file.path(vdir, sprintf("verify_index_alignment%s.csv", .tsuf)), bom = TRUE)

cat(sprintf("\n==== verify_index_alignment%s: %s ====\n", .tsuf, if (ok_all) "PASS" else "FAIL"))
if (!ok_all) cat(sprintf("  Investigate before trusting outputs/%s/attributes%s.\n", VERSION, .tsuf))
