# Korea forest attribute maps — analysis code

Code behind the nationwide 30 m dataset of four forest attributes for South Korea, and the analyses
reported in the accompanying data descriptor.

**Data** <https://doi.org/10.5281/zenodo.21701424> (CC BY 4.0) · **Code** MIT (this repository)

**Viewer** <https://ee-korea-forest-maps.projects.earthengine.app/view/korea-forest-maps> — the maps in a browser, two at a time under a wipe.

## What it produces

Four attributes — coniferous fraction, dominant tree species, aboveground biomass, and tree-species
diversity — mapped two ways, each under two predictor sets, at 30 m on EPSG:5179:

- **direct prediction** — one random-forest model per attribute, plus a per-cell uncertainty layer
- **plot-index imputation** — `yaImpute` nearest-neighbour, assigning each forest cell one matched
  inventory plot so its attributes are a combination observed together in the field

Reference data is the 8th Korean National Forest Inventory (2021–2025).

## Running the pipeline

`pipeline/` is the run order. Every filename carries its stage number, so the directory listing is
the sequence: reference data, covariates, predictors, folds, products, post-processing, validation,
deposit. [`docs/pipeline.md`](docs/pipeline.md) walks through it with the exact calls.

R steps run in RStudio; Python steps from a shell at the repository root. Google Earth Engine is
used for two steps only, both Sentinel-2 compositing.

```
pipeline/    the numbered stages, 01 to 08       R/         shared helpers
figures/     the three figures made from code    data/      small lookup tables
app/         the Earth Engine viewer             docs/      pipeline.md, the run order
scripts/     package restore
```

Scripts resolve their inputs from the repository root, either through `here()` or by walking up from
their own file, so run them where they sit rather than copying one out on its own.

### Input data

This repository holds code, not data. It expects the following under `data/raw/`, which is not tracked:

| input | source |
|---|---|
| National Forest Inventory plot records | Korea Forest Service, via the public NFI database at <https://kfss.forest.go.kr/> |
| Forest-type map (1:5,000) | Korea Forest Service, Forest Geospatial Information System |
| Land cover (Level 3, 1:5,000) | Ministry of Climate, Energy and Environment, EGIS |
| Climate (MK-PRISM, GK2A) | Korea Meteorological Administration |
| Topography | NASA SRTM, SRTMGL1 v3 |
| Sentinel-2 | European Space Agency Copernicus, via Google Earth Engine |

None of these are ours to redistribute. Each must be obtained from the provider named above.

**Inventory plot coordinates are confidential and are in neither this repository nor the deposited
dataset.** `.gitignore` blocks any `*_coords.csv` from being committed.

Set `EE_PROJECT` to your own Earth-Engine-enabled Cloud project before running either compositing
step.

### Package versions

R packages are pinned in `renv.lock`. Restore them with `renv::restore()`. `scripts/init_renv.R` is
the script that produced the lockfile; re-running it re-pins against whatever is installed, so use it
only when you mean to update the pins.

The Python steps need `numpy`, `pandas`, `openpyxl`, `rasterio` and `pyproj`, plus `earthengine-api`
for the two compositing steps.

## Checking the outputs

A GDAL statistics write can silently corrupt a GeoTIFF: part of the tile index in the
full-resolution image is overwritten while the overviews survive, so a GIS opens the file normally
and a header-only check reports it as sound. Run `pipeline/08a_verify_raster_integrity.py` before
trusting or publishing any raster, and do not conclude one is undamaged from a sampled read — a
corrupted file here passed a 164-block sample.

## Citation

Kim, S., Shin, J. H., Han, H., and Choe, H.: Nationwide 30 m maps of forest composition, biomass, and
diversity in South Korea (2021–2025) from direct prediction and plot-index imputation [Data set], Zenodo,
https://doi.org/10.5281/zenodo.21701424, 2026.

## Acknowledgments

Google Gemini, OpenAI ChatGPT, and Anthropic Claude Code were used to assist with code refinement
and debugging.
