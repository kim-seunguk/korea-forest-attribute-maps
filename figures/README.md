# Figures

| script | what it shows |
|---|---|
| [`fig05_accuracy.R`](fig05_accuracy.R) | predicted against observed, both product families |
| [`fig08_calibration.R`](fig08_calibration.R) | how well the uncertainty layers are calibrated |
| [`figA1_confusion.R`](figA1_confusion.R) | the dominant-species confusion matrix 07d writes, an optional view |

The workflow diagram and the maps have no script here: the diagram is drawn in presentation
software, and the study-area and national maps are composed in GIS software over the released
rasters.

`theme.R` is shared. Source the figure script, which sources `theme.R` itself.

```r
source(here::here("figures", "fig05_accuracy.R"))
```

`fig05_accuracy.R` and `fig08_calibration.R` read validation tables for the delivered predictor
set. `fig05_accuracy.R` reads `uncertainty_points_ftm.csv`, which `07g_uncertainty.R` writes when it
is run with `FTM_ONLY <- TRUE`, and `family_comparison_points.csv`, which `07b_compare_families.R`
writes. `fig08_calibration.R` reads `uncertainty_calibration_ftm.csv`, which the same `07g` run
writes.

Output goes to `figures/rendered/`, which is not included here. The scripts read their tables from
`outputs/nfi8/`. To draw the `fig08_calibration.R` figure without the pipeline, copy
`validation/uncertainty_calibration_ftm.csv` to `outputs/nfi8/` first. `fig05_accuracy.R` reads the
two per-plot tables, which are not distributed.

## Conventions

- `fig05_accuracy.R` computes the R², RMSE and MBE it prints from the two per-plot tables.
- `ggplot2` with `patchwork`, on the shared `theme_fig()`.
- Color-blind-safe palettes throughout. The species colors come from `data/species_colors.csv`,
  the table the delivered species map and the viewer are colored by.
- Captions live in the manuscript, not in the image.
- Species are labeled by scientific name.
