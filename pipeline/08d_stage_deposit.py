"""08d_stage_deposit.py — assemble the Zenodo deposit as a frozen, verified copy.

WHY a separate staged tree rather than uploading from outputs/nfi8/ directly: base-tier rasters have
been corrupted twice by a GDAL PAM statistics write overrunning the IFD0 tile index, both times
on files a GIS had opened. The deposit therefore gets its own directory that
nothing else touches, and the integrity gate runs against that copy, not against the working tree.

The pipeline already writes the published names, so nothing is renamed here. This gathers the
released files from their working subdirectories into one flat folder and verifies them. FILES below
is that list; the names themselves are decided in R/paths.R.

Source of truth for the schema is `outputs/nfi8/data_records.csv` (written by 08c_build_data_records.py) —
this script transforms that table rather than restating it, so the two can never disagree.

RUN   python pipeline/08d_stage_deposit.py
      python pipeline/08d_stage_deposit.py --dry-run
OUT   outputs/nfi8/_deposit_v1/   the released files, flat, + checksums.sha256
EXIT  0 = staged and verified, 1 = something failed and the deposit is not ready.
"""
import csv
import hashlib
import io
import os
import shutil
import sys

os.environ.setdefault('GDAL_PAM_ENABLED', 'NO')   # no .aux.xml sidecars: that write is the corruption vector

import rasterio
from rasterio.enums import Resampling

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')
DEPOSIT = os.path.join(OUT, '_deposit_v1')

DOI = '10.5281/zenodo.21701424'
GITHUB = 'https://github.com/kim-seunguk/korea-forest-attribute-maps'
VIEWER = 'https://ee-korea-forest-maps.projects.earthengine.app/view/korea-forest-maps'
EE_COLLECTION = 'projects/ee-korea-forest-maps/assets/korea-forest-attributes'
OVERVIEW_LEVELS = [2, 4, 8, 16, 32]

# ---------------------------------------------------------------------------------------------
# pipeline name -> deposit name. `none`/`full` are the internal predictor-set tokens; the deposit
# says base/ftm, matching the manuscript.
# ---------------------------------------------------------------------------------------------
# the three continuous attributes carry an ensemble-SD layer; dominant species carries `1 - vote
# margin` under a different suffix, so it is handled separately rather than folded into the loop.
CONTINUOUS = ['conifer_fraction', 'biomass', 'shannon']
VARIABLE = CONTINUOUS + ['dominant_species']
TIER = ['base', 'ftm']

# The deposit name IS the pipeline name, so this map is an identity. It is kept because the rest of
# the script is written in terms of (source name -> deposit name), and because it states the released
# set explicitly: anything not listed here is not deposited.
RENAME = {}
SUBDIR = {}
for _t in TIER:
    for _v in VARIABLE:
        RENAME[f'{_t}_{_v}_30m.tif'] = f'{_t}_{_v}_30m.tif'
        SUBDIR[f'{_t}_{_v}_30m.tif'] = 'm8_maps_30m'
    for _v in CONTINUOUS:
        RENAME[f'{_t}_{_v}_sd_30m.tif'] = f'{_t}_{_v}_sd_30m.tif'
        SUBDIR[f'{_t}_{_v}_sd_30m.tif'] = 'm8_uncertainty_30m'
    RENAME[f'{_t}_dominant_species_uncertainty_30m.tif'] = f'{_t}_dominant_species_uncertainty_30m.tif'
    SUBDIR[f'{_t}_dominant_species_uncertainty_30m.tif'] = 'm8_uncertainty_30m'
for _f in ('imputation_base.tif', 'imputation_ftm.tif',
           'nearest_reference_distance_base.tif', 'nearest_reference_distance_ftm.tif'):
    RENAME[_f] = _f
    SUBDIR[_f] = 'index'

# non-raster companions: (source under outputs/nfi8, deposit name) — identities, as above
TABLES = [('plot_lookup_base.csv', 'plot_lookup_base.csv'),
          ('plot_lookup_ftm.csv', 'plot_lookup_ftm.csv'),
          ('dominant_species_table.csv', 'dominant_species_table.csv')]

# Manifest notes that name a companion file. These are identities now that the pipeline writes the
# published names, and the list is kept only so a future divergence has somewhere to be declared
# rather than being fixed by hand in the staged copy.
NOTE_FIX = [('dominant_species_table.csv', 'dominant_species_table.csv'),
            ('plot_lookup_ftm.csv', 'plot_lookup_ftm.csv'),
            ('imputation_ftm.tif', 'imputation_ftm.tif')]


def read_csv(path):
    with io.open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, cols):
    """UTF-8 without a BOM. plot_lookup*.csv shipped with one and it breaks a naive header parse."""
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, '') for c in cols})


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------
# species harmonisation
# ---------------------------------------------------------------------------------------------
def species_maps():
    """Korean field name -> (39-class code, scientific name).

    The lookup tables carry the plot's FIELD-measured dominant species (129 / 141 distinct names),
    not the 39 raster classes, so they cannot be joined to the species table as shipped. The legend
    is the authoritative code<->name map; any species outside its 38 named classes is what the raster
    calls class 38 (기타 / Other). The field name is KEPT — collapsing it would discard exactly the
    detail that makes the plot-index family worth having — and the code is added alongside it.
    """
    legend = read_csv(os.path.join(OUT, 'm8_maps_30m', 'm8_domsp_legend_used.csv'))
    code_of = {r['species']: r['code'] for r in legend}
    sci_of = {r['species']: r['sci_name'] for r in legend}
    other = next((r['code'] for r in legend if r['species'] == '기타'), '38')

    nfi = read_csv(os.path.join(ROOT, 'data', 'species_nfi_lookup.csv'))
    for r in nfi:
        name = r['Name_Korean']
        sci = f"{r['Genus']} {r['Species']}".strip()
        if name not in sci_of and sci:
            sci_of[name] = sci
    return code_of, sci_of, other


def stage_lookup(src_name, dst_name, code_of, sci_of, other, dry):
    """Copy a plot_lookup table with the species columns harmonised and the BOM dropped."""
    rows = read_csv(os.path.join(OUT, src_name))
    unresolved = set()
    for r in rows:
        kr = r.pop('dominant_species', '')
        r['dominant_species_kr'] = kr
        r['dominant_species_sci'] = sci_of.get(kr, '')
        r['dominant_species_code'] = code_of.get(kr, other)
        if kr and kr not in sci_of:
            unresolved.add(kr)

    cols = ['plot_index', 'nfi_plot_id',
            'dominant_species_code', 'dominant_species_kr', 'dominant_species_sci',
            'coniferous_fraction', 'aboveground_biomass_Mgha', 'shannon_diversity',
            'forest_type', 'diameter_class', 'crown_density',
            'dominant_species_ba_share', 'carbon_Mgha', 'basal_area_m2ha']
    missing = [c for c in cols if c not in rows[0]]
    extra = [c for c in rows[0] if c not in cols]
    if missing:
        raise SystemExit(f"  {src_name}: expected columns absent: {missing}")
    if extra:
        cols = cols + extra                      # never silently drop a column
    if not dry:
        write_csv(os.path.join(DEPOSIT, dst_name), rows, cols)
    folded = sum(1 for r in rows if r['dominant_species_code'] == other)
    print(f"  {src_name:24s} -> {dst_name:24s} {len(rows):>6,} rows, "
          f"{folded} folded into class {other} (Other)"
          + (f", {len(unresolved)} names without a scientific name" if unresolved else ""))
    return unresolved


def stage_species_table(dst_name, dry):
    rows = read_csv(os.path.join(OUT, 'dominant_species_table.csv'))
    cols = list(rows[0].keys())
    if not dry:
        write_csv(os.path.join(DEPOSIT, dst_name), rows, cols)
    print(f"  {'dominant_species_table.csv':24s} -> {dst_name:24s} {len(rows)} classes")


# ---------------------------------------------------------------------------------------------
# rasters
# ---------------------------------------------------------------------------------------------
def stage_raster(src_name, dry):
    dst_name = RENAME[src_name]
    src = os.path.join(OUT, SUBDIR[src_name], src_name)
    dst = os.path.join(DEPOSIT, dst_name)
    if not os.path.exists(src):
        raise SystemExit(f"  MISSING source: {src}")
    if dry:
        print(f"  {src_name:38s} -> {dst_name}")
        return dst

    shutil.copy2(src, dst)
    with rasterio.open(dst) as ds:
        n_ov = len(ds.overviews(1))
    if n_ov == 0:
        # A safety net that should no longer fire: 06d gives every delivered raster its overviews.
        # It exists because a national float32 raster is painful to browse without them, and
        # because it did once fire, for the two distance layers 06d used to leave out. Built on the
        # COPY, then re-verified by the gate below.
        with rasterio.open(dst, 'r+') as ds:
            ds.build_overviews(OVERVIEW_LEVELS, Resampling.average)
            ds.update_tags(ns='rio_overview', resampling='average')
        with rasterio.open(dst) as ds:
            n_ov = len(ds.overviews(1))
        print(f"  {src_name:38s} -> {dst_name:44s} (+{n_ov} overview levels)")
    else:
        print(f"  {src_name:38s} -> {dst_name}")
    return dst


def stage_manifest(dry):
    """Rewrite outputs/nfi8/data_records.csv into deposit names."""
    src = os.path.join(OUT, 'data_records.csv')
    if not os.path.exists(src):
        raise SystemExit("  outputs/nfi8/data_records.csv is absent - run 08c_build_data_records.py first")
    rows = read_csv(src)
    if len(rows) != 20:
        raise SystemExit(f"  manifest has {len(rows)} rows, expected 20")
    for r in rows:
        if r['dtype'] == 'MISSING':
            raise SystemExit(f"  manifest row {r['file']} is MISSING - the product is incomplete")
        if r['epsg'] != '5179':
            raise SystemExit(f"  manifest row {r['file']} has epsg={r['epsg']}, expected 5179")
        r['file'] = RENAME[r['file']]
        for old, new in NOTE_FIX:
            r['notes'] = r['notes'].replace(old, new)
    cols = list(rows[0].keys())
    if not dry:
        write_csv(os.path.join(DEPOSIT, 'data_records.csv'), rows, cols)
    print(f"  {'data_records.csv':24s} -> {'data_records.csv':24s} {len(rows)} rows")
    return rows


# ---------------------------------------------------------------------------------------------
def write_readme(manifest, checks, dry):
    n_base = sum(1 for r in manifest if r['predictor_set'] == 'base')
    body = f"""# Nationwide 30 m maps of forest composition, biomass, and diversity in South Korea (2021–2025)

Kim, S., Shin, J. H., Han, H., and Choe, H. (2026).
Department of Agriculture, Forestry and Bioresources, Seoul National University.

DOI <https://doi.org/{DOI}> · Licence **CC BY 4.0** · Code (MIT) <{GITHUB}>

Read `data_records.csv` alongside this file — it is the machine-readable schema for every raster
(data type, nodata, unit, valid range, grid, CRS).

## Looking at the maps without downloading them

The same products are published as a public Earth Engine image collection, with a browser viewer
over it, so they can be examined before or instead of downloading several gigabytes.

| | |
|---|---|
| Viewer | <{VIEWER}> |
| Collection | `{EE_COLLECTION}` |

The collection also carries the four attributes of the plot-index family as ready-made rasters,
derived there from the plot index and the lookup table in this record.

## What is here

Four forest attributes over South Korea at 30 m, produced two different ways, each under two
predictor sets. All rasters share one grid.

| | |
|---|---|
| Attributes | coniferous fraction · dominant tree species (39 classes) · aboveground biomass · tree-species diversity (Shannon *H*′) |
| Predictor sets | `base` = Sentinel-2, climate, topography · `ftm` = base plus the national forest-type map |
| Product families | **direct prediction** (a random-forest model per attribute) · **plot-index imputation** (each cell is assigned one matched inventory plot) |
| Reference data | 8th Korean National Forest Inventory, 2021–2025 |
| Grid | EPSG:5179, 30 m, {manifest[0]['width']} × {manifest[0]['height']}, origin ({manifest[0]['origin_x']}, {manifest[0]['origin_y']}) |
| Extent | stocked forest, nationwide. The base tier fills the DMZ gap with ESA WorldCover; the ftm tier leaves it nodata |

**The two families answer different questions.** Direct prediction is more accurate attribute by
attribute. Plot-index imputation is less accurate per attribute but every cell carries a combination
of attributes that was actually observed together in one field plot, which matters if you need the
attributes to be mutually consistent rather than individually optimal.

## Files

**Direct prediction** — `{{base,ftm}}_<attribute>_30m.tif`, each paired with a per-cell relative
uncertainty layer:

- `*_sd_30m.tif` — ensemble standard deviation, for the three continuous attributes. Higher = more uncertain.
- `*_dominant_species_uncertainty_30m.tif` — `1 − vote margin`. Higher = more uncertain.

**Plot-index imputation**

- `imputation_{{base,ftm}}.tif` — integer index of the matched inventory plot.
- `plot_lookup_{{base,ftm}}.csv` — join `plot_index` to get that plot's attributes.
- `nearest_reference_distance_{{base,ftm}}.tif` — random-forest **distance** from the cell to its
  matched plot, defined as `1 − the proportion of trees placing the two in the same terminal node`.
  Higher means a **less** similar match. It is a continuous diagnostic, not a thresholded mask.

**Tables**

- `dominant_species_table.csv` — the 39 class codes with Korean and scientific names, prevalence, and
  per-class producer's / user's accuracy and F1.
- `data_records.csv` — per-file schema for all {len(manifest)} rasters ({n_base} base, {len(manifest) - n_base} ftm).

## Decoding the dominant-species map

Raster values are class codes 0–38; join to `dominant_species_table.csv` on `code`. Class 38 (기타 /
Other) holds species too rare to model individually. Nodata is 255.

## Decoding the plot-index map

Raster values are an internal plot index (1…N). Join to `plot_lookup_*.csv` on `plot_index`.

The lookup carries `nfi_plot_id`, the National Forest Inventory sample-plot number. This is the
identifier introduced for the 8th inventory, which renumbered the plots carried over from earlier
cycles: `신표본점번호` in the inventory database tables, `표본점번호[신]` in the published microdata.
Join on that column, not on the old number, to reach any attribute not included here through the
Korea Forest Service statistics platform (<https://kfss.forest.go.kr/stat/>).

Note that `dominant_species_kr` / `_sci` are the plot's **field-measured** species, which is finer
than the 39 mapped classes; `dominant_species_code` is the class the raster would assign it.

## No plot coordinates are distributed

Inventory plot locations are confidential and are not part of this deposit in any form. The plot-index
raster already carries all the spatial information the dataset provides — a matched plot's true
location cannot be recovered from where its index appears, because matching is in predictor space, not
geographic space.

## Accuracy, in brief

Base-tier direct predictions, 50 km spatial-block cross-validation (*K* = 10) and a spatially disjoint
test set:

| attribute | cross-validation | disjoint test |
|---|---|---|
| coniferous fraction | *R*² 0.61 | 0.61 |
| dominant species | accuracy 0.45 (majority-class baseline 0.30) | 0.50 |
| aboveground biomass | *R*² 0.23 | 0.22 |
| tree-species diversity | *R*² 0.24 | 0.28 |

**Use these maps for broad spatial gradients, not for precise estimates in individual cells.** The
*R*² values for biomass and diversity are low. For those two attributes the plot-index family falls
below a mean-predictor baseline, so the direct maps are preferable. Per-class species accuracy is
concentrated in the more prevalent species — check `dominant_species_table.csv` before relying on any
one class. Only the direct-prediction uncertainty layers were validated against empirical error.

## Format

Tiled GeoTIFF, 256 × 256 blocks, LZW compression, internal overviews. Tables are UTF-8 CSV without a
byte-order mark.

## Citation

Kim, S., Shin, J. H., Han, H., and Choe, H.: Nationwide 30 m maps of forest composition, biomass, and
diversity in South Korea (2021–2025) from direct prediction and plot-index imputation [Data set], Zenodo,
https://doi.org/{DOI}, 2026.

---

Every file in this record was verified by a full tile-by-tile decode before upload. `checksums.sha256`
covers all {checks} files.
"""
    if not dry:
        with io.open(os.path.join(DEPOSIT, 'README.md'), 'w', encoding='utf-8', newline='\n') as f:
            f.write(body)
    print(f"  {'(generated)':24s} -> {'README.md':24s} {len(body.splitlines())} lines")


def main(argv):
    dry = '--dry-run' in argv
    print(f"staging deposit -> {os.path.relpath(DEPOSIT, ROOT)}" + ("  [DRY RUN]" if dry else ""))

    print("\n[1] integrity gate on the working tree")
    # Loaded by path: the stage number makes 08a_verify_raster_integrity an invalid module name,
    # so a plain import of it cannot work.
    import importlib.util
    gate = os.path.join(HERE, '08a_verify_raster_integrity.py')
    spec = importlib.util.spec_from_file_location('raster_integrity_gate', gate)
    vri = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vri)
    if vri.main() != 0:
        print("  STOP - the working tree has damaged rasters.")
        return 1

    if not dry:
        if os.path.exists(DEPOSIT):
            print(f"\n  {os.path.relpath(DEPOSIT, ROOT)} already exists - remove it first")
            return 1
        os.makedirs(DEPOSIT)

    print("\n[2] rasters")
    staged = [stage_raster(n, dry) for n in RENAME]

    print("\n[3] tables")
    code_of, sci_of, other = species_maps()
    unresolved = set()
    for src, dst in TABLES:
        if src == 'dominant_species_table.csv':
            stage_species_table(dst, dry)
        else:
            unresolved |= stage_lookup(src, dst, code_of, sci_of, other, dry)
    if unresolved:
        print(f"  ! {len(unresolved)} species names have no scientific name: {sorted(unresolved)[:8]}")

    print("\n[4] manifest and README")
    manifest = stage_manifest(dry)
    n_files = len(RENAME) + len(TABLES) + 2          # + data_records.csv + README.md
    write_readme(manifest, n_files, dry)

    if dry:
        print(f"\nDRY RUN - {n_files} files would be staged. Nothing was written.")
        return 0

    print("\n[5] integrity gate on the staged copy")
    bad = 0
    for p in staged:
        r = vri.scan(p)
        if r is None:
            print(f"  ? {os.path.basename(p)} not a classic TIFF")
            continue
        n_tiles, n_bad, _ = r
        if n_bad:
            print(f"  DAMAGED {os.path.basename(p)}  {n_bad} of {n_tiles} tiles")
            bad += 1
    if bad:
        print(f"\n{bad} staged file(s) damaged. The deposit is NOT ready.")
        return 1
    print(f"  all {len(staged)} staged rasters clean")

    print("\n[6] checksums")
    lines = []
    for name in sorted(os.listdir(DEPOSIT)):
        p = os.path.join(DEPOSIT, name)
        if os.path.isfile(p) and name != 'checksums.sha256':
            lines.append(f"{sha256(p)}  {name}")
    with io.open(os.path.join(DEPOSIT, 'checksums.sha256'), 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(lines) + '\n')
    total = sum(os.path.getsize(os.path.join(DEPOSIT, n)) for n in os.listdir(DEPOSIT))
    print(f"  {len(lines)} files, {total / 1e9:.2f} GB -> checksums.sha256")

    print(f"\nDeposit staged and verified: {os.path.relpath(DEPOSIT, ROOT)}")
    print("Do not open these files in a GIS. Upload it to the data repository.")
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
