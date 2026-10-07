# Figures

Two of the paper's figures are generated from code, and one more script is kept as an optional view.
The scripts keep the names they were released under, so the number in a script's name is not the
number the figure carries in the paper; the table gives both.

| script | number in the paper | what it shows |
|---|---|---|
| [`fig05_accuracy.R`](fig05_accuracy.R) | 4 | predicted against observed, both product families |
| [`fig08_calibration.R`](fig08_calibration.R) | 5 | how well the uncertainty layers are calibrated |
| [`figA1_confusion.R`](figA1_confusion.R) | not in the paper | the dominant-species confusion matrix 07d writes, an optional view |

The paper's other figures have no script here: the workflow diagram is drawn in presentation
software, and the study area, with the plots drawn as 1 km cluster points, and the national maps are
composed in GIS software over the released rasters.

`theme.R` and `species_palette.R` are shared. Source the figure script; it pulls in the theme itself.

```r
source(here::here("figures", "fig05_accuracy.R"))
```

Both paper figures read the delivered predictor set's tables: `fig05_accuracy.R` needs
`m7_uncertainty_points_ftm.csv` from `07g_uncertainty.R` run with `FTM_ONLY <- TRUE` and
`m4c_pred_points.csv` from `07b_compare_families.R`; `fig08_calibration.R` needs
`m7_calibration_ftm.csv` from the same 07g run.

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
