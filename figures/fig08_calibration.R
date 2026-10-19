# fig08_calibration.R — reliability curves for the four uncertainty layers.
# ============================================================================================
# One 2x2 reliability figure from outputs/nfi8/uncertainty_calibration_ftm.csv (07g with FTM_ONLY):
#   (a) coniferous fraction · (c) aboveground biomass · (d) Shannon diversity
#       -> predicted ensemble SD vs mean absolute error (monotonic increase = reliable uncertainty)
#   (b) dominant species -> vote margin vs observed accuracy
#   source(here::here("figures","fig08_calibration.R"))
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2); library(patchwork) })
source(here::here("figures","theme.R"))
VERSION <- "nfi8"
outd <- here("figures","rendered"); dir.create(outd, showWarnings=FALSE, recursive=TRUE)
cal <- fread(here("outputs", VERSION, "uncertainty_calibration_ftm.csv"))

PRETTY <- c(coniferfrac_full="Coniferous fraction", Bio_a_ha="Aboveground biomass",
            shannon_full="Tree-species diversity")
small_theme <- theme_fig(12) +
  theme(plot.title=element_text(face="bold", size=12, hjust=0, margin=margin(b=6)),
        plot.margin=margin(8,10,4,8))

# ---- (a, c, d) continuous targets: predicted SD vs mean absolute error ----------------------------
cont_panel <- function(RT, tag) {
  d <- cal[metric=="tree_sd" & target==RT]
  ggplot(d, aes(mean_uncert, mean_abs_err)) +
    geom_line(colour="#2166AC", linewidth=0.8) +
    geom_point(colour="#2166AC", size=2.2) +
    scale_x_continuous(expand=expansion(mult=0.06)) +
    scale_y_continuous(expand=expansion(mult=0.08)) +
    labs(title=sprintf("(%s) %s", tag, PRETTY[[RT]]),
         x="Predicted ensemble SD (decile mean)", y="Mean absolute error") +
    small_theme
}
# ---- (b) dominant species: vote margin vs observed accuracy ------------------------------------
cm <- cal[metric=="vote_margin"]
p_d <- ggplot(cm, aes(mean_uncert, accuracy)) +
  geom_line(colour="#B2182B", linewidth=0.8) +
  geom_point(colour="#B2182B", size=2.2) +
  scale_y_continuous(limits=c(0,1)) +
  labs(title="(b) Dominant species",
       x="Vote margin (decile mean)", y="Observed accuracy") +
  small_theme

fig <- (cont_panel("coniferfrac_full","a") | p_d) /
       (cont_panel("Bio_a_ha","c")          | cont_panel("shannon_full","d"))

fn <- file.path(outd, "fig08_calibration.png")
save_png_pdf(fn, fig, width=8.6, height=8.0, dpi=600)
cat(sprintf("wrote %s  (2x2 reliability curves; continuous SD-vs-error + domsp margin-vs-accuracy)\n", fn))
