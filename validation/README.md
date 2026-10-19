# Validation outputs

This directory holds the tables behind the accuracy figures in the accompanying paper, as the
pipeline wrote them. They are here so that a user can check any number without re-running the
pipeline, which needs the inventory plot data and a long compute run.

Every file is a plain UTF-8 CSV. **None contains a plot identifier or a coordinate**; the inventory
plot locations are confidential and are not distributed in any form. The per-plot tables the same
stages write beside these (predictions, uncertainty and local spacing at each plot) carry plot
identifiers and are not here.

## Files

| file | what it holds | produced by |
|---|---|---|
| `predictor_ladder_ftmonly.csv` | the direct random forests under the delivered predictor set: the four attributes, three validation schemes, cluster-bootstrap intervals. The headline table | `07a_predictor_ladder.R` with `FTM_ONLY` |
| `predictor_ladder_baseonly.csv` | the same for the two sets without the forest-type map, `clim` and `clim_s2med`, on their own complete-case sample | `07a_predictor_ladder.R` with `BASE_ONLY` |
| `predictor_ladder_rich.csv` | `clim_s2med` against `clim_s2rich` on their common sample | `07a_predictor_ladder.R` with `RICH_ONLY` |
| `family_comparison.csv` | the plot-index imputation scored on the same folds: the four attributes, the field structure variables it matches on, crown density at the 50 % cut, the forest-type-map baselines, and the split by stand origin | `07b_compare_families.R` |
| `continuous_metrics_wholeftm.csv` | RMSE, relative RMSE and mean bias of the three continuous direct maps | `07d_extra_metrics.R` with `FTM_ONLY` |
| `stratified_accuracy_ftm.csv` | direct-map accuracy by stand origin: the pooled model evaluated within each stratum, and models fitted within one stratum | `07f_stratified_accuracy.R` |
| `uncertainty_calibration_ftm.csv` | mean absolute error, or species accuracy, by decile of the per-cell uncertainty layers | `07g_uncertainty.R` |
| `uncertainty_coverage_ftm.csv` | coverage of the 50, 80 and 95 % intervals built from the ensemble SD, overall and by biomass third | `07h_uncertainty_coverage.py` |
| `temporal_sensitivity_ftm.csv` | direct-map accuracy by the plot's survey year | `07i_temporal_sensitivity.py` |
| `attribute_consistency_ftm.csv` | variance ratios, attribute-space plausibility distances and Bhattacharyya overlap of both families against the plot population | `07j_attribute_consistency.py` |
| `dominant_species_table_ftm.csv` | the 37 trained classes: code, names, reference-plot count and prevalence, and producer's and user's accuracy and F1 under spatial-block cross-validation. The table deposited beside the maps | `07k_species_reliability.py` |
| `map_summary_plotuse_ftm.csv` | how the plot-index map uses the reference library: mapped cells, plots assigned, cells and hectares per plot | `07l_map_summaries.py` |
| `support_scale_ftm.csv` | accuracy against aggregation support, from the cell to 100 km, with a permutation control and a minimum number of plots per unit | `07o_support_scale.R` |
| `cluster_support_ftm.csv` | both families at the support of the plot and at the support of the inventory cluster | `07p_cluster_support.R` |
| `density_weighted.csv` | accuracy with each plot weighted by the area its cluster represents, beside the unweighted figure | `07q_density_weighting.py` |
| `block_size_test.csv` | a climate-only biomass model under spatial blocks of 15 to 100 km, beside random folds | `07r_block_size_test.R` |
| `variogram_range.csv` | semivariogram ranges of plot biomass, basal area and the dominant species' share of basal area, with the variogram's own block size (`block_km_variogram`) and the block size used (`block_km_used`) | `04a_spatial_folds.R` |
| `env_repr_quantiles.csv`, `env_repr_pc1_ratio.csv` | how the plots cover the climate and terrain of the mapped forest: sample and area quantiles with a Kolmogorov–Smirnov distance per variable, and the sample-to-area ratio by decile of the first principal component | `07s_environmental_representativeness.R` |

The texture and temporal Sentinel-2 features behind `predictor_ladder_rich.csv` came from an
exploratory export outside the released pipeline.

## Reading them

Rows are long-format. The columns that recur across files:

- `target` — `coniferfrac_full` (coniferous fraction of basal area), `domsp_full` (dominant species,
  pooled to the 37-class legend of the map), `Bio_a_ha` (aboveground biomass, Mg ha⁻¹) and
  `shannon_full` (Shannon *H*′). `family_comparison.csv` adds `domsp_raw`, the plot's field-recorded
  dominant species without pooling, which is what the delivered index decodes to; `fld_V_VG`,
  `fld_V_DC` and `fld_V_DN`, the field forest type, diameter class and crown-density class the
  imputation matches on; `fld_V_DN2`, crown density at the one class boundary the field record and
  the map share, the 50 % cut; and `*_ftmbaseline`, the forest-type map's own class at the plot
  scored as if it were a prediction (`predset` `FTM`).
- `predset` — `full` is the delivered set: Sentinel-2 seasonal composites, climate, terrain and
  coordinates, and the forest-type map's forest type, diameter class, density class and stand
  origin. `clim_s2med` is the same without the forest-type map; `clim` is climate, terrain and
  coordinates only; `clim_s2rich` is `clim_s2med` plus the texture and temporal Sentinel-2 features.
  `clim_nolonlat`, in `block_size_test.csv`, is `clim` without the planar coordinates and the aspect
  terms.
- `scheme` — `spatial` is the 50 km spatial-block cross-validation (K = 10), `random` is ten random
  folds over the same plots, `test` is the spatially disjoint hold-out. Blocks, folds and the
  hold-out are assigned by inventory cluster, so the subplots of one cluster are always on the same
  side of a split.
- `route` — the product family, `direct` or `impute` / `imputation`.
- `point` or `r2` — the coefficient of determination for a continuous target, top-1 accuracy for a
  categorical one. `value` and `window` are diagnostics written beside it, accuracy within a
  tolerance (the top three classes for species) and R² after averaging to the 50 km block; they are
  not comparable across targets. `nir` is the majority-class rate.
- `ci_lo`, `ci_hi` — a 95 % bootstrap interval that resamples inventory clusters, not plots, with
  seed 7214.
- `n` — the plots in that row. The delivered set is scored on one complete-case sample: 13,023 plots
  under cross-validation and 2,989 in the hold-out, 13,013 and 2,986 for biomass, which a few plots
  lack. The sets without the forest-type map and the richer Sentinel-2 comparison each have a
  complete-case sample of their own, given in their `n`.
- `stratum` — `full` is every plot; `natural` and `planted` split by stand origin where a stage does
  so; `-` marks a row the stage does not stratify.

Columns particular to one file:

- `stratified_accuracy_ftm.csv` — `model` `pooled` is one forest fitted on every plot and evaluated
  within `eval_stratum`; `within` is fitted and evaluated inside `fit_stratum`. `var` is the
  variance of the observed target in that evaluation set.
- `uncertainty_calibration_ftm.csv` — `metric` `tree_sd` is the ensemble standard deviation of the
  continuous layers, `vote_margin` the species layer; `bin` is the decile of the uncertainty value,
  1 the most certain; `accuracy` is filled for species, `mean_abs_err` for the rest.
- `uncertainty_coverage_ftm.csv` — `stratum` is `all`, a third of observed biomass (`low_third`,
  `mid_third`, `high_third`) or a third of predicted biomass (`pred_low_third` and so on).
- `support_scale_ftm.csv` — `support_km` 0 is the cell itself; `n_units` the aggregation units
  holding at least `min_plots` plots; `ctl_mean`, `ctl_lo` and `ctl_hi` the permutation control,
  the same statistic after the plots are shuffled among units, over `nperm` permutations;
  `admissible` is 1 where the bootstrap interval (`ci_width`) is narrow enough for the row to be
  reported.
- `cluster_support_ftm.csv` — `support` `1_plot` trains and evaluates at the plot, the delivered
  product; `2_plot_eval_cluster` trains at the plot and averages the subplots' out-of-fold
  predictions inside each cluster; `3_cluster` trains and evaluates at the cluster, one prediction
  per cluster. `varratio` is the variance of the predictions relative to `var_obs`, the variance of
  the observations.
- `density_weighted.csv` — `unweighted` is the published figure; `weighted` gives each plot the
  square of its cluster's local grid spacing (1, 2 or 4 km) as a weight, the area it represents.
- `block_size_test.csv` — biomass only; `gap` is `random_r2` minus `spatial_r2`, `n_blocks` the
  number of blocks at that size.
- `env_repr_quantiles.csv` — `samp_*` are quantiles over the plots, `area_*` over a random sample of
  mapped forest cells, `KS_D` the Kolmogorov–Smirnov distance between the two distributions.
  `env_repr_pc1_ratio.csv` — `repr_ratio` is `samp_share` over `area_share`; above 1 the decile is
  over-sampled, below 1 under-sampled.
- `dominant_species_table_ftm.csv` — `users_acc` is empty for the 14 classes the classifier never
  assigned, where it is undefined rather than zero.
- `attribute_consistency_ftm.csv` — the `note` column defines each metric.

## The validation design

Random folds place training and validation plots in the same neighborhood, so they measure spatial
interpolation as well as prediction, and figures from different schemes are not comparable. With the
clusters kept whole, the gap between random and spatial-block folds under the delivered set is 0.011
to 0.023 across the four attributes (`predictor_ladder_ftmonly.csv`, `random` minus `spatial`). A
model built from climate, terrain and coordinates alone scores an R² of 0.07 for biomass under
random folds and 0.02 under blocks (`predictor_ladder_baseonly.csv`, `predset` `clim`);
`block_size_test.csv` gives the gap across block sizes.

The inventory is not a uniform 4 km grid: 24.7 % of the plots sit on the 1 and 2 km supplementary
spacing of the metropolitan areas and Jeju, so an unweighted figure counts those regions more
heavily per unit of forest. `density_weighted.csv` reports every figure with each plot weighted
by the area it represents; the change is at most 0.016.

## Where the headline numbers come from

| number | value | file | row |
|---|---|---|---|
| direct maps, spatial-block cross-validation: coniferous fraction *R*², species accuracy, biomass *R*², Shannon *R*² | 0.72, 0.52, 0.31, 0.32 | `predictor_ladder_ftmonly.csv` | `scheme` `spatial`, column `point` |
| the same, disjoint test | 0.70, 0.57, 0.32, 0.35 | `predictor_ladder_ftmonly.csv` | `scheme` `test` |
| the same, random folds | 0.73, 0.55, 0.34, 0.34 | `predictor_ladder_ftmonly.csv` | `scheme` `random` |
| majority-class baseline for species | 0.29 (cross-validation), 0.33 (test) | `predictor_ladder_ftmonly.csv` | `domsp_full`, column `nir` |
| RMSE of the continuous maps, cross-validation | 0.18, 63.9 Mg ha⁻¹, 0.42 | `continuous_metrics_wholeftm.csv` | `scheme` `spatial`, column `rmse` |
| without the forest-type map, cross-validation | 0.70, 0.49, 0.27, 0.30 | `predictor_ladder_baseonly.csv` | `predset` `clim_s2med`, `scheme` `spatial` |
| climate, terrain and coordinates only, cross-validation | 0.13, 0.36, 0.02, 0.03 | `predictor_ladder_baseonly.csv` | `predset` `clim`, `scheme` `spatial` |
| largest change from the richer Sentinel-2 features | 0.004 (cross-validation), 0.009 (test) | `predictor_ladder_rich.csv` | `clim_s2rich` minus `clim_s2med`, same target and scheme |
| plot-index maps, cross-validation | 0.48, 0.42 (0.42 unpooled), −0.29, −0.27 | `family_comparison.csv` | `predset` `full`, `stratum` `full`, `scheme` `spatial`; `domsp_raw` for the unpooled species |
| plot-index maps, disjoint test | 0.47, 0.48 (0.47 unpooled), −0.28, −0.19 | `family_comparison.csv` | as above, `scheme` `test` |
| the matching variables, plot index against the majority class, cross-validation and test | forest type 0.46 and 0.50 against 0.31 and 0.25; diameter class 0.47 and 0.48 against 0.47 and 0.47; crown density 0.82 and 0.79 against 0.84 and 0.81, at the 50 % cut 0.91 against 0.92 | `family_comparison.csv` | `fld_V_VG`, `fld_V_DC`, `fld_V_DN`, `fld_V_DN2`, `predset` `full`, `stratum` `full`, columns `point` and `nir` |
| direct maps by stand origin, cross-validation | coniferous fraction 0.74 natural, 0.61 planted; biomass 0.28, 0.37; Shannon 0.34, 0.22 | `stratified_accuracy_ftm.csv` | `model` `pooled`, `eval_stratum`, `scheme` `spatial` |
| plot-index coniferous fraction by stand origin, cross-validation | 0.52 natural, 0.28 planted | `family_comparison.csv` | `predset` `full`, `stratum` `natural` / `planted` |
| accuracy after aggregation to coarser resolution, cross-validation | 1 km: 0.78, 0.36, 0.39; 10 km: 0.85, 0.47, 0.55 | `support_scale_ftm.csv` | `support_km` 1 and 10, `min_plots` 1, column `r2` |
| permutation control up to 10 km | −0.64 to −0.10 | `support_scale_ftm.csv` | `min_plots` 1, column `ctl_mean` |
| reference noise: all 100 m units against those holding at least two plots | 0.74 to 0.80, 0.32 to 0.37, 0.34 to 0.41 | `support_scale_ftm.csv` | `support_km` 0.1, `min_plots` 1 and 2 |
| interval coverage at 50, 80 and 95 % | coniferous fraction 0.53, 0.87, 0.96; biomass 0.52, 0.83, 0.96; Shannon 0.50, 0.81, 0.96 | `uncertainty_coverage_ftm.csv` | `stratum` `all` |
| error against the uncertainty layer | monotone over the deciles for coniferous fraction and biomass, not for Shannon; species accuracy 0.24 to 0.93 | `uncertainty_calibration_ftm.csv` | `mean_abs_err` by `bin`; `accuracy` for `domsp_full` |
| accuracy by survey year | biomass 0.25 to 0.35 | `temporal_sensitivity_ftm.csv` | `Bio_a_ha` |
| variance of the map over the variance of the plots | direct 0.66, 0.32, 0.31; plot index 1.06, 1.04, 1.10 | `attribute_consistency_ftm.csv` | `metric` `variance_ratio` |
| overlap of the joint attribute distribution with the plots | direct 0.70, plot index 0.99 | `attribute_consistency_ftm.csv` | `metric` `bhattacharyya_joint` |
| mapped cells; reference plots assigned at least one cell | 61,587,083; 15,999 of 15,999 | `map_summary_plotuse_ftm.csv` | `mapped_cells`, `plots_assigned_at_least_one_cell` |
| classes trained; classes never assigned | 37; 14 | `dominant_species_table_ftm.csv` | rows; empty `users_acc` |
| cluster support, direct maps | coniferous fraction 0.72 to 0.76, biomass 0.31 to 0.35, Shannon 0.32 to 0.41 | `cluster_support_ftm.csv` | `route` `direct`, `support` `1_plot` against `3_cluster` |
| semivariogram ranges of the plot attributes | 6 to 15 km (biomass 10.6, basal area 6.4, dominant species' share 14.9) | `variogram_range.csv` | column `range_raw_km` |
| environmental representativeness | Kolmogorov–Smirnov distance up to 0.12; sample-to-area ratio 0.75 to 1.41 across the deciles | `env_repr_quantiles.csv`, `env_repr_pc1_ratio.csv` | `KS_D`; `repr_ratio` |
