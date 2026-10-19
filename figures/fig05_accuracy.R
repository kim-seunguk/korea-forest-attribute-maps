# fig05_accuracy.R — predicted against observed, for both product families.
# ================================================================================================
# The figure has two rows, direct prediction and plot-index imputation, and three columns:
# coniferous fraction, aboveground biomass and tree-species diversity. Both rows use the delivered
# predictor set and spatial-block cross-validation.
#
# Each panel carries the 1:1 line, an OLS fit, and R2, RMSE and MBE, with %RMSE added for biomass.
#
# In : outputs/nfi8/uncertainty_points_ftm.csv      (direct, spatially out-of-fold; 07g with FTM_ONLY)
#      outputs/nfi8/family_comparison_points.csv    (imputation; 07b, predset full)
#
#   source(here::here("figures", "fig05_accuracy.R"))
# ================================================================================================
suppressPackageStartupMessages({
  library(here); library(data.table); library(ggplot2); library(patchwork); library(grid)
})
margin <- ggplot2::margin
source(here::here("figures","theme.R"))
VERSION <- "nfi8"; outd <- here("outputs", VERSION)

TARGETS <- list(
  list(key = "coniferfrac_full", title = "Coniferous fraction",                 pct = FALSE, lim = c(0, 1)),
  # NB: plotmath, not a literal ⁻ - Arial has no U+207B glyph and renders it as a tofu box.
  # bold() is required: a plotmath title ignores theme_fig's plot.title face = "bold".
  list(key = "Bio_a_ha",         title = expression(bold("Aboveground biomass (Mg ha"^"-1"*")")), pct = TRUE,  lim = NULL),
  list(key = "shannon_full",     title = "Tree-species diversity (Shannon H′)", pct = FALSE, lim = NULL))

## ---- both families' per-plot spatial-CV predictions (the delivered predictor set) ------------
dir_f <- file.path(outd, "uncertainty_points_ftm.csv")
if (!file.exists(dir_f))
  stop("this figure needs ", dir_f, " — run pipeline/07g_uncertainty.R with FTM_ONLY <- TRUE first.")
d_dir <- fread(dir_f)[target %in% c("coniferfrac_full", "Bio_a_ha", "shannon_full")]  # 07g cols: target/Plot_ID/pred/obs (no predset)
d_dir <- d_dir[, .(route = "Direct prediction", target,
                   obs = suppressWarnings(as.numeric(obs)), pred = suppressWarnings(as.numeric(pred)))]
imp_f <- file.path(outd, "family_comparison_points.csv")
if (!file.exists(imp_f))
  stop("this figure needs ", imp_f, " — run pipeline/07b_compare_families.R, which writes it.")
d_imp <- fread(imp_f)[predset == "full"]                                # full = the delivered set, 07b's token
d_imp <- d_imp[, .(route = "Plot-index imputation", target,
                   obs = suppressWarnings(as.numeric(obs)), pred = suppressWarnings(as.numeric(pred)))]
D <- rbind(d_dir, d_imp)[is.finite(obs) & is.finite(pred)]

r2  <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
fmt <- function(x, d = 2) formatC(ifelse(abs(x) < 0.5 * 10^-d, 0, x), format = "f", digits = d)   # no "-0.00"
# Each panel has a bold letter (a to f) at its top left and its attribute as a column title. Each
# row is headed by its method.
panel <- function(rt, spec, show_title, tag) {
  d <- D[route == rt & target == spec$key]
  if (!nrow(d)) { warning("no rows for ", rt, " / ", spec$key); return(patchwork::plot_spacer()) }
  R2 <- r2(d$obs, d$pred); rmse <- sqrt(mean((d$pred - d$obs)^2)); mbe <- mean(d$pred - d$obs)
  lab <- sprintf("R² = %s\nRMSE = %s\nMBE = %s%s", fmt(R2), fmt(rmse), fmt(mbe),
                 if (spec$pct) sprintf("\n%%RMSE = %s%%", fmt(100 * rmse / mean(d$obs), 0)) else "")
  lim <- if (is.null(spec$lim)) range(c(d$obs, d$pred), na.rm = TRUE) else spec$lim
  g <- ggplot(d, aes(obs, pred)) +
    geom_point(alpha = 0.10, size = 0.32, colour = "#2166AC", shape = 16) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey35", linewidth = 0.5) +
    geom_smooth(method = "lm", se = FALSE, colour = "#B2182B", linewidth = 0.8, formula = y ~ x) +
    annotate("text", x = lim[2], y = lim[1], hjust = 1, vjust = 0, size = 3, colour = "grey20", label = lab) +   # stats: BOTTOM-RIGHT
    annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 5, fontface = "bold", label = tag) +   # letter: TOP-LEFT (opposite corner -> never overlaps the stats)
    coord_equal(xlim = lim, ylim = lim) + labs(x = "Observed", y = "Predicted") +
    theme_fig(12) + theme(plot.margin = margin(6, 6, 4, 4),           # top/bottom breathing room between the two rows
                         plot.title = element_text(size = 11, hjust = 0.5))
  if (show_title) g + labs(title = spec$title) else g
}

## ---- assemble: a–f tagged panels; column titles on both rows; descriptive row headers --------
row_hdr <- function(txt) wrap_elements(full = textGrob(txt, x = 0.004, hjust = 0,
             gp = gpar(fontface = "bold", fontsize = 12)))
dP <- Map(function(s, t) panel("Direct prediction",     s, TRUE, t), TARGETS, c("a", "b", "c"))
iP <- Map(function(s, t) panel("Plot-index imputation", s, TRUE, t), TARGETS, c("d", "e", "f"))   # column titles on BOTH rows
fig <- row_hdr("Direct prediction")     / (dP[[1]] | dP[[2]] | dP[[3]]) /
       row_hdr("Plot-index imputation") / (iP[[1]] | iP[[2]] | iP[[3]]) +
       plot_layout(heights = c(0.05, 1, 0.05, 1))

outdir <- here("figures","rendered"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fn <- file.path(outdir, "fig05_accuracy.png")
save_png_pdf(fn, fig, width = 10.5, height = 7.8, dpi = 600)
cat(sprintf("wrote %s  (2x3 a–f: direct vs plot-index imputation x 3 continuous targets; the delivered set, spatial-CV)\n", fn))
