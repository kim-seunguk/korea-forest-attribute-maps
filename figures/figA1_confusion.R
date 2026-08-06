# figA1_confusion.R — dominant-species confusion matrix, row-normalised.
# ============================================================================================
# Reads the 07d confusion matrix (observed x predicted counts) and draws a row-normalized heatmap
# (each row = an observed species, cells = fraction predicted as each class), species labelled with
# scientific names and ordered by observed abundance.
#   source(here::here("figures","figA1_confusion.R"))          # as published
#   TIER <- "full"; source(here::here("figures","figA1_confusion.R"))   # the FTM-added tier
# ============================================================================================
suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2); library(patchwork) })
margin <- ggplot2::margin
source(here::here("figures","theme.R"))
VERSION <- "nfi8"
# The defaults reproduce the figure as published. Override any of the three before sourcing.
if (!exists("TIER"))   TIER   <- "base"      # base | full  (predictor set)
if (!exists("SCHEME")) SCHEME <- "spatial"   # spatial | test
# SUF picks the evaluation SAMPLE, which is a different axis from TIER (TIER = predictor set):
#   ""            = common subset, n = 12,032 / 2,348 (base and FTM-added comparable like for like)
#   "_wholebase"  = every plot with a complete base predictor record, n = 12,926 / 2,521 (base predictor
#                   set only; the manuscript standard). "whole", not "full", so it cannot be read as the
#                   FTM-added predictor set.
if (!exists("SUF"))    SUF    <- "_wholebase"
outd <- here("outputs", VERSION)

cmfile <- file.path(outd, sprintf("m4f_domsp_confusion_%s_%s%s.csv", TIER, SCHEME, SUF))
cm <- fread(cmfile, encoding = "UTF-8")                            # first col = obs (Korean), rest = predicted classes
setnames(cm, 1, "obs")

# ---- Korean -> scientific name (legend); "기타" -> "Other" ---------------------------------
leg <- fread(here("data","species_codes.csv"), encoding = "UTF-8")
kr2sci <- setNames(leg$sci_name, leg$species); kr2sci["기타"] <- "Other"
sci <- function(k) ifelse(is.na(kr2sci[k]) | k == "Other", ifelse(k == "기타", "Other", k), kr2sci[k])

# ---- long + row-normalize ------------------------------------------------------------------
long <- melt(cm, id.vars = "obs", variable.name = "pred", value.name = "n")
long[, pred := as.character(pred)]
long[, rowsum := sum(n), by = obs]
long[, frac := ifelse(rowsum > 0, n / rowsum, 0)]

# order both axes fully taxonomically: CLASS (문/강: angiosperms then gymnosperms) -> ORDER (목) -> FAMILY
# (과) -> GENUS (속) -> SPECIES (종), all alphabetical; the residual "Other" is forced LAST. Because
# "Magnoliopsida" < "Pinopsida", plain alphabetical on Class puts all angiosperm orders first, then the
# gymnosperm orders (Cupressales, Pinales) as a trailing conifer block. Two auditable lookups drive this:
#   data/species_nfi_lookup.csv  (Shin NFI: Korean -> Family/Genus/Species; authoritative)
#   data/family_order_lookup.csv (Family -> Order + Class; APG IV orders, Kew-consistent)
# genus+species come from the domsp legend scientific name (family name only drives the join, not display).
tax  <- fread(here("data","species_nfi_lookup.csv"),  encoding = "UTF-8")
ford <- fread(here("data","family_order_lookup.csv"), encoding = "UTF-8")
g2f  <- setNames(tax$Family, tax$Genus)                         # genus  -> family (Shin NFI)
f2o  <- setNames(ford$Order, ford$Family)                       # family -> order
f2c  <- setNames(ford$Class, ford$Family)                       # family -> class (angio/gymno)
scinm  <- sci(cm$obs); genus <- sub(" .*$", "", scinm)
family <- unname(g2f[genus]); order_ <- unname(f2o[family]); class_ <- unname(f2c[family])
is_oth <- (scinm == "Other") | is.na(family)                   # 기타 or unknown-family -> forced last
if (any(is_oth & scinm != "Other"))
  message("figS1: no family for: ",  paste(unique(scinm[is_oth & scinm != "Other"]), collapse = ", "))
if (any(!is.na(family) & is.na(order_)))
  message("figS1: no order for family: ", paste(unique(family[!is.na(family) & is.na(order_)]), collapse = ", "))
ord_kr <- cm$obs[order(is_oth, paste(class_, order_, family, scinm))]  # Class -> Order -> Family -> Genus -> species; Other last
lab <- setNames(sci(ord_kr), ord_kr)
long[, obs  := factor(obs,  levels = ord_kr, labels = lab[ord_kr])]
long[, pred := factor(pred, levels = ord_kr, labels = lab[ord_kr])]

p <- ggplot(long, aes(pred, obs, fill = frac)) +
  geom_tile() +
  geom_tile(data = long[as.character(obs) == as.character(pred)],           # emphasise the 1:1 diagonal
            fill = NA, colour = "#B2182B", linewidth = 0.45) +              # (correct predictions); red contrasts the dark mako
  scale_fill_viridis_c(name = "Row fraction (observed → predicted)", option = "G", direction = -1,
    limits = c(0, 1), breaks = c(0, 0.5, 1),
    guide = guide_colorbar(title.position = "top", ticks.colour = "grey30")) +
  scale_x_discrete(position = "top") + scale_y_discrete(limits = rev) +
  labs(x = "Predicted", y = "Observed") + coord_equal() + theme_fig(11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6, face = "italic"),
        axis.text.y = element_text(size = 6, face = "italic"), panel.grid = element_blank(),
        legend.position = "bottom", legend.key.width = unit(18, "mm"), legend.key.height = unit(4, "mm"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7))

# ---- row marginal: observed plots per species (n, sqrt scale), aligned to the heatmap rows ----------
# The row-normalized matrix hides sample size; this bar shows how many plots each observed row rests on
# (P. densiflora ~thousands vs rare species ~tens), so the reader can tell reliable rows from noisy ones.
rown <- long[, .(n = sum(n)), by = obs]                           # obs already carries the taxonomic factor + sci labels
pbar <- ggplot(rown, aes(n, obs)) +
  geom_col(fill = "grey55", width = 0.82) +
  scale_x_sqrt(expand = expansion(mult = c(0, 0.05)), breaks = c(0, 500, 2000, 4000)) +
  scale_y_discrete(limits = rev) +                                # SAME rev as the heatmap -> rows align
  labs(x = "Plots (n)", y = NULL) + theme_fig(11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 6))
fig <- p + pbar + plot_layout(widths = c(1, 0.2), guides = "collect") & theme(legend.position = "bottom")

outdir <- here("figures","rendered"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fn <- file.path(outdir, sprintf("figA1_confusion_%s_%s%s.png", TIER, SCHEME, SUF))
# NB with coord_equal the heatmap panel is square; the canvas WIDTH sets how tight the bar sits against it.
# If a gap remains, DROP `width`; if the bar overlaps the matrix, RAISE it. (height fixes the square size.)
save_png_pdf(fn, fig, width = 9, height = 8.5, dpi = 600)
cat(sprintf("wrote %s  (tier=%s scheme=%s, %d classes; row marginal = observed plots, sqrt)\n", fn, TIER, SCHEME, length(ord_kr)))
