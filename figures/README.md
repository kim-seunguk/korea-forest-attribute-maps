# Figures

Three figures are generated from code. Each is named for the number it carries in the paper.

| script | figure | what it shows |
|---|---|---|
| [`fig05_accuracy.R`](fig05_accuracy.R) | 5 | predicted against observed, both product families |
| [`fig08_calibration.R`](fig08_calibration.R) | 8 | how well the uncertainty layers are calibrated |
| [`figA1_confusion.R`](figA1_confusion.R) | A1 | dominant-species confusion matrix |

The remaining figures in the paper are maps and schematics composed in GIS and drawing software over
the released rasters, so they have no script here.

`theme.R` and `species_palette.R` are shared. Source the figure script; it pulls in the theme itself.

```r
source(here::here("figures", "fig05_accuracy.R"))
```

Output goes to `figures/rendered/`, which is not included here. The figures regenerate from the
released data, so only the code is kept.

## Conventions

- Figures read from the CSV and raster outputs on disk. They never recompute a statistic, so a
  number in a figure and the same number in a table cannot disagree.
- `ggplot2` with `patchwork`, on the shared `theme_fig(base_size = 11)`.
- Colour-blind-safe palettes throughout. `species_palette.R` holds the species symbology, read from
  the same table the published maps are coloured by, so the legends cannot drift apart. It is used
  by the map figures, which are composed in GIS rather than here.
- Panel labels are lower case in brackets, `(a)`, `(b)`. Legends sit inside the panel.
- Captions live in the manuscript, not in the image.
- Species are labelled by scientific name.
