# theme.R — shared figure theme, palettes and save helpers.
# --------------------------------------------------------------------------
# One source of style for every figure generated from code, so the set is
# visually consistent: theme_bw as the base, no minor grid, bold left-aligned
# panel tags, in-panel statistics, patchwork for assembly, and fixed physical
# dimensions at export. Figure text is English and species are labelled by
# scientific name.
#
# Usage:
#   source(here::here("figures", "theme.R"))
#   ... build a ggplot `p` ...
#   save_png_pdf(file.path(here::here("figures", "rendered"), "fig05_accuracy.png"),
#                p, width = 8.5, height = 3.6)
#
# This file defines functions and constants only. It reads no data and draws
# nothing.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

# ---- 1. Base theme ------------------------------------------------------- #
# theme_bw as the spine, with the invariant overrides on top: no minor grid,
# bold left-aligned panel title, tight margins, clean legend.
theme_fig <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0,
                                       size = rel(1.0)),
      plot.subtitle    = element_text(hjust = 0, size = rel(0.85),
                                       colour = "grey30"),
      plot.margin      = margin(3, 4, 3, 3),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      legend.position  = "bottom",
      legend.key.size  = unit(0.4, "cm"),
      legend.title     = element_text(size = rel(0.85)),
      legend.text      = element_text(size = rel(0.8))
    )
}

# Journal-style panel tag, e.g. panel_tag("a", "Species top-1") -> "(a) Species top-1"
panel_tag <- function(letter, title = NULL) {
  if (is.null(title)) sprintf("(%s)", letter) else sprintf("(%s) %s", letter, title)
}

# ---- 2. Colours ---------------------------------------------------------- #
# Species colours live in species_palette.R, which reads the same table the maps
# are symbolised from, so a species has one colour everywhere.

# ---- 3. Save helper ------------------------------------------------------ #
# Rendered figures land in figures/rendered/, at 600 dpi.
#
# Write the same figure as PNG *and* PDF. Journals prefer vector art, and these panels are text/line
# work (points, tiles, labels) that stays sharp and small as PDF; only the raster MAP figures
# (the maps) are composed in ArcGIS and exported from there. cairo_pdf, not the default pdf(),
# so UTF-8 glyphs (×, ′, ⁻¹) and fonts embed correctly on Windows. `png_path` may be given with or
# without the .png suffix; the PDF is written beside it.
save_png_pdf <- function(png_path, plot, width, height, dpi = 600) {
  png_path <- sub("\\.png$", "", png_path)
  for (ext in c(".png", ".pdf")) {
    ggsave(paste0(png_path, ext), plot, width = width, height = height, bg = "white",
           dpi = dpi, units = "in", device = if (ext == ".pdf") grDevices::cairo_pdf else NULL)
  }
  cat(sprintf("wrote %s.{png,pdf}\n", png_path))
  invisible(paste0(png_path, c(".png", ".pdf")))
}
