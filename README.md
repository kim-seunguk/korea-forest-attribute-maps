# Korea forest attribute maps — analysis code

Code behind the nationwide 30 m dataset of four forest attributes for South Korea, and the analyses
reported in the accompanying paper.

**Data** <https://doi.org/10.6084/m9.figshare.34558446> (CC BY 4.0) ·
**Code** MIT (this repository)

**Viewer** <https://ee-korea-forest-maps.projects.earthengine.app/view/korea-forest-maps> — the maps in a browser, two at a time under a wipe.

## What it produces

Four attributes — coniferous fraction, dominant tree species (37 classes), aboveground biomass, and
tree-species diversity — mapped two ways from one predictor set (Sentinel-2 seasonal composites,
climate, topography and the national forest-type map), at 30 m on EPSG:5179:

- **direct prediction** — one random-forest model per attribute, plus a per-cell uncertainty layer
- **plot-index imputation** — `yaImpute` nearest-neighbor, assigning each forest cell one matched
  inventory plot so its attributes are a combination observed together in the field

Reference data is the 8th Korean National Forest Inventory (2021–2025). Accuracy is estimated by
50 km spatial-block cross-validation and a spatially disjoint test set, both with the subplots of an
inventory cluster kept together; coniferous fraction and dominant species are meant to be read cell
by cell, biomass and tree-species diversity at 1 km or coarser.

## Running the pipeline

Filenames in `pipeline/` carry their stage number, so the directory listing is the run order:
reference data, covariates, predictors, folds, products, post-processing, validation,
deposit. [`docs/pipeline.md`](docs/pipeline.md) walks through it with the exact calls.

R steps run in RStudio, and Python steps from a shell at the repository root.

```
pipeline/    the numbered stages, 01 to 08       R/         shared helpers
figures/     the three figures made from code    data/      small lookup tables
app/         the Earth Engine viewer             docs/      pipeline.md, the run order
scripts/     package restore
validation/  the tables behind the reported accuracy
```

Scripts resolve their inputs from the repository root, through `here()` or by walking up from their
own file, so run them where they sit.

`pipeline/08d_stage_deposit.py` copies the delivered rasters and tables into the deposit under the
deposit's file names (its `RENAME` table).

### Input data

This repository holds code, not data. It expects the following under `data/raw/`, which is not tracked:

| input | source |
|---|---|
| National Forest Inventory plot records | Korea Forest Service, via the public NFI database at <https://kfss.forest.go.kr/> |
| Forest-type map (1:5,000) | Korea Forest Service, Forest Geospatial Information System |
| Land cover (Level 3, 1:5,000) | Ministry of Climate, Energy and Environment, EGIS |
| Climate (MK-PRISM, GK2A) | Korea Meteorological Administration |
| Topography | NASA SRTM, SRTMGL1 v3, exported onto the product grid via Google Earth Engine |
| Sentinel-2 | European Space Agency Copernicus, via Google Earth Engine |

The lookup tables in `data/` are tracked here and read by the pipeline: `01d` writes
`species_forest_type_lookup.csv` from the curated species table, and the others hold the 37-class
species legend and its colors, the recode of the forest-type map's codes and the species taxonomy.

**Inventory plot coordinates are confidential and are in neither this repository nor the deposited
dataset.** `.gitignore` blocks any `*_coords.csv` from being committed.

Set `EE_PROJECT` to your own Earth-Engine-enabled Cloud project before running either Earth Engine
export.

### Package versions

R 4.6.0, with the packages pinned in `renv.lock`. Restore them with `renv::restore()`.
`scripts/init_renv.R` produced the lockfile. Re-running it re-pins against whatever is installed,
so use it only to update the pins.

Python 3.11 with numpy 2.4.4, pandas 3.0.3, openpyxl 3.1.5, rasterio 1.4.4 and pyproj 3.7.2, plus
earthengine-api 1.7.32 for the two Earth Engine exports and the optional Earth Engine publisher.

## Checking the outputs

A GDAL statistics write can silently corrupt a GeoTIFF: part of the tile index in the
full-resolution image is overwritten while the overviews survive, so a GIS opens the file normally
and a header-only check reports it as sound. Run `pipeline/08a_verify_raster_integrity.py` before
trusting or publishing any raster, and do not conclude a raster is undamaged from a sampled read.

## Citation

Kim, S., Shin, J. H., Han, H., and Choe, H.: Nationwide 30 m maps of forest composition, biomass, and
diversity in South Korea (2021–2025) from direct prediction and plot-index imputation
[Data set], figshare, https://doi.org/10.6084/m9.figshare.34558446, 2026.

## Acknowledgments

OpenAI ChatGPT and Anthropic Claude Code were used to assist with code refinement and debugging.
