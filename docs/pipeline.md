# Pipeline

The pipeline has eight stages, from inventory records to the deposited rasters. The filenames in
`pipeline/` carry their stage number, so the directory listing is the run order.

R steps run in RStudio with `source(...)`. Python steps run from a shell at the repository root.
Google Earth Engine is used for two exports only: the Sentinel-2 seasonal composite (`03a`) and the
SRTM elevation resampled onto the product grid (`02f`).

Everything is written under `outputs/nfi8/` and `data/interim/`, neither of which is in this
repository. Stage the source data under `data/raw/` first; the README lists what is needed and where
to obtain it.

## Setup

```r
renv::restore()                   # install the versions pinned in renv.lock
```

## 01 · Reference data

| step | script |
|---|---|
| extract the inventory tree records | `01a_extract_nfi_records.R` |
| build the per-plot reference table | `01b_nfi_to_points.R` |
| recover each plot's survey year | `01c_recover_survey_year.R` |
| build the species to forest-type lookup | `01d_species_forest_type_lookup.py` |

For plots surveyed from 2022 on, `01b` takes the stand table's coordinate where it differs from the
register's by more than 50 m. The reference table holds 17,099 plots.

## 02 · Covariate rasters

`02a_forest_type_map_to_rasters.R` (stand attributes) · `02b_forest_origin_to_raster.R` (planted or
natural) · `02c_climate_terrain_to_rasters.R` (climate and terrain) · `02d_forest_mask.R` (the
stocked-forest extent) · `02e_landcover_to_raster.R` (independent land cover, used in stage 06) ·
`02f_srtm_to_grid.py` (Earth Engine: SRTM elevation exported once onto the product grid)

`02c` runs twice: once for the climate layers, which define the grid `02f` exports onto, and again
after the SRTM export has been downloaded into `data/raw/dem_srtm_grid/`, for the terrain layers.
Elevation is therefore resampled once, and slope and aspect are derived on the product grid.

`02e` is the exception to the run order. It belongs here by subject, but it takes its grid template
from a delivered raster, so run it after stage 05. Nothing before stage 06 reads what it writes.

Stages 02a and 02e read tens of gigabytes of vectors. Both are resumable, and both are meant to be
run one edition at a time.

**Manual step.** `02d` writes a mask that nothing else in the pipeline reads. Upload it to Earth
Engine as an Image asset and set `EE_FTM_ASSET` to the asset id. `03a` masks the Sentinel-2
composites with this asset, so that they cover the same extent as the predictions.

## 03 · Predictors

```
03a_sentinel2_composites.py         Earth Engine: the seasonal composites, area-averaged onto the product grid
03b_align_sentinel2.R               put the downloaded tiles on the covariate grid
03c_split_sentinel2_bands.R         one file per band, which the raster package reads quickly in stage 05
03d_extract_predictors_at_plots.R   climate, terrain and stand attributes at the plots
03e_sentinel2_at_plots.R            the 42 Sentinel-2 bands at the plots, read from the aligned raster
```

`03e` reads the plot features from the raster `03b` aligned, so the models train on the raster they
predict onto. `s2_composite.py` holds the one definition of the composite, which `03a` imports. It
reads nothing and writes nothing. `R/paths.R` holds the same 42 band names for the R side.

`03a` writes to Drive, and downloading is manual: the raster tiles go into `data/raw/gee_s2med/` for
`03b`, or into the folder named in `S2_RAW_DIR`.

## 04 · Folds and targets

`04a_spatial_folds.R` freezes the split every later stage reuses: 50 km spatial blocks, K = 10, plus
a spatially disjoint test set of about 20 % of plots (13,923 train, 3,176 test). Blocks, folds and
the test split are assigned by inventory cluster, so the subplots of one cluster never straddle a
block edge or the split. `04b_build_targets.R` builds the four targets from basal area weighted by
each tree's expansion factor.

## 05 · The delivered products

**Direct prediction** — one random forest per attribute, plus a per-cell uncertainty layer. Species
dominant on fewer than 20 reference plots are pooled into the other class, code 36, which leaves 37
classes with codes 0 to 36.

```r
source(here::here("pipeline", "05a_direct_prediction.R"))
```

**Plot-index imputation** — one matched inventory plot per cell, then decoded into attributes. The
reference library holds the 15,999 plots with a complete predictor record and a biomass value.

```r
source(here::here("pipeline", "05b_plot_index_imputation.R"))
source(here::here("pipeline", "05c_index_to_attributes.R"))
```

The index raster is written as a 16-bit unsigned integer with 0 as nodata, so GIS software treats it as
categorical by default.

Keep `IMPUTE_CORES` at 16 or below. Larger SOCK clusters can stall.

Then build the lookup, the deposited table that gives each integer in the index raster its meaning.

```r
source(here::here("pipeline", "05d_build_plot_lookup.R"))
```

## 06 · Post-processing

`06a_landcover_crosstab.R` → `06b_apply_landcover_filter.R` (in place, both families) →
`06d_finalize_rasters.py` → `06e_postprocess.R`

`05a` writes the published class codes, so there is no recoding step.

`06b` must be re-run whenever stage 05 rewrites a raster, because the filter is lost on overwrite.

## 07 · Validation

```
07a_predictor_ladder.R            07b_compare_families.R            07c_imputation_validation.R
07d_extra_metrics.R               07e_forest_type_map_comparison.R  07f_stratified_accuracy.R
07g_uncertainty.R                 07h_uncertainty_coverage.py       07i_temporal_sensitivity.py
07j_attribute_consistency.py      07k_species_reliability.py        07l_map_summaries.py
07m_verify_index_alignment.R      07n_verify_biomass.R              07o_support_scale.R
07p_cluster_support.R             07q_density_weighting.py          07r_block_size_test.R
07s_environmental_representativeness.R
```

`07a` and `07b` score the two product families on the same folds. `07a` is run once per
predictor set (`FTM_ONLY` for the delivered set, `BASE_ONLY` for the sets without the forest-type
map, `RICH_ONLY` for the richer Sentinel-2 features), and each run writes its own table. `07f`,
`07g`, `07o` and `07p` take the same `FTM_ONLY` switch and write `_ftm`-suffixed tables, `07d` writes
`_wholeftm` under it, and `07h`/`07i` read whichever set the `PREDSET` environment variable names
(default `ftm`). Bootstrap intervals resample inventory clusters. The tables behind the paper's
numbers are copied to `validation/` and described in
[`validation/README.md`](../validation/README.md).

## 08 · Data records and deposit

```
python pipeline/08a_verify_raster_integrity.py          # must exit 0 before anything is published
python pipeline/08c_build_data_records.py
python pipeline/08d_stage_deposit.py --dry-run           # every check, nothing written
DEPOSIT_DOI=<DOI> VIEWER_URL=<url> python pipeline/08d_stage_deposit.py
python pipeline/08e_publish_earthengine.py               # optional: the public Earth Engine copy
```

`08a` parses the tile tables of every delivered raster. `--decode` also reads every block through
GDAL, and `--inputs` extends the scan to the predictor rasters under `data/interim/`. `08d` stages
the deposit into `outputs/nfi8/_deposit_v2/`: it checks that every source belongs to the current
run, copies and verifies the rasters and tables, generates the README from the validation tables and
writes the checksums. The dataset DOI and the viewer URL come from the environment, and
a real run stops without them.

`08e` is optional and needs two things this repository does not carry: an Earth-Engine-enabled
Cloud project in `EE_PROJECT` and a staging bucket in `GCS_BUCKET`. `app/` holds the viewer over
the collection `08e` publishes.

## Figures

[`figures/README.md`](../figures/README.md) lists the figures generated from code and how to run
them.

## Shared helpers

`R/paths.R` holds the layout, the reference configuration and the product file naming.
`R/eval_utils.R` holds the accuracy metrics, the cluster-level blocks and folds and the cluster
bootstrap interval, and `R/run_utils.R` the run sizing and timing.
