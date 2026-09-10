# Validation outputs

The numbers reported in the paper's Technical Validation section, as the files the pipeline wrote.
They are here so a reader can check any figure without re-running the pipeline, which needs the
inventory plot data and several days of compute.

Every file is a plain UTF-8 CSV. **None contains a plot identifier or a coordinate**; the inventory
plot locations are confidential and are not distributed in any form.

| file | what it holds | produced by |
|---|---|---|
| `m4a_ladder.csv` | accuracy of the direct-prediction maps for every predictor set and every validation scheme (spatial-block, random, disjoint test), with bootstrap intervals | `pipeline/07a_predictor_ladder.R` |
| `m4c_dualtrack.csv` | the same targets scored for the plot-index family, plus the three field structure variables the matching is performed on and the forest-type-map baselines for them | `pipeline/07b_compare_families.R` |
| `m9_consistency.csv` | variance ratios, Bhattacharyya distribution overlap and attribute-space plausibility distance for both product families | `pipeline/07j_attribute_consistency.py` |
| `m11_support_scale.csv` | accuracy as a function of aggregation support, from 100 m to 100 km, with bootstrap intervals and a permutation control | `pipeline/07o_support_scale.R` |
| `m7_calibration.csv` | mean absolute error by decile of the per-cell uncertainty layers | `pipeline/07g_uncertainty.R` |
| `m7_coverage.csv` | interval coverage against nominal 50, 80 and 95 %, overall and by biomass third | `pipeline/07h_uncertainty_coverage.py` |
| `m8_map_summary_plotuse.csv` | how many cells each reference plot is assigned in the plot-index map | `pipeline/07l_map_summaries.py` |

## Reading them

Rows are long-format. The columns that appear across files:

- `target` — `coniferfrac_full`, `domsp_full`, `Bio_a_ha`, `shannon_full`, and the field structure
  variables `fld_V_VG`, `fld_V_DC`, `fld_V_DN`
- `predset` — `clim_s2med` is the base predictor set and `full` is the FTM-added set; `clim` and
  `clim_nolonlat` are diagnostic subsets that are not released as products
- `scheme` — `spatial` is the 50 km spatial-block cross-validation, `random` is random folds over the
  same plots, `test` is the spatially disjoint hold-out
- `point` / `r2` — the coefficient of determination for continuous targets, overall accuracy for
  categorical ones; `nir` is the majority-class baseline
- `ci_lo` / `ci_hi` — 95 % bootstrap interval
- `n` — plots contributing to that row. Sample sizes differ between rows because each predictor set
  masks to its own complete cases; the paper's Methods gives the lineage.

## One thing to know before comparing rows

`scheme` matters more than any other column. Random folds place training and validation plots in the
same neighbourhood on a systematic 4 km inventory grid, so they measure spatial interpolation as well
as prediction. The same biomass model reads 0.21 under spatial blocks and 0.33 under random folds, and
a model built from climate, terrain and coordinates alone, with no satellite data, reaches 0.31 under
random folds and 0.01 under blocks. Figures from different schemes are not comparable.
