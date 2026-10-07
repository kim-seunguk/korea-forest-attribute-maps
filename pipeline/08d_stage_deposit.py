"""08d_stage_deposit.py — assemble version 2 of the data deposit as a frozen, verified copy.

Version 2 replaces every file of version 1 (10.5281/zenodo.21701424): one predictor set, 37
dominant-species classes, corrected plot positions, rebuilt inputs and validation by inventory
cluster. The staged folder is uploaded as a new version of the same record, so the version DOI the
repository reserves for that draft is passed in through NFIMAP_DEPOSIT_DOI and written into the
README; version 1 stays available under its own DOI.

WHY a separate staged tree rather than uploading from outputs/nfi8/ directly: rasters of this dataset
have been corrupted by a GDAL PAM statistics write overrunning the IFD0 tile index, on files a GIS had
opened. The deposit therefore gets its own directory that nothing else touches, and the integrity
gate runs against that copy, not against the working tree.

The pipeline writes the published raster names (naming is decided once, in R/paths.R). Only the tables
whose working names carry the tier suffix are renamed on the way in.

SOURCES
    m8_maps_30m/ftm_<attribute>_30m.tif · m8_uncertainty_30m/ftm_*_{sd,uncertainty}_30m.tif
    index/imputation_ftm.tif · index/nearest_reference_distance_ftm.tif
    data_records.csv                08c_build_data_records.py      -> data_records.csv
    plot_lookup_ftm.csv             05d_build_plot_lookup.R        -> plot_lookup_ftm.csv (+ species columns)
    dominant_species_table_ftm.csv  07k_species_reliability.py     -> dominant_species_table.csv

CURRENCY. Tables from an earlier run can lie beside the current rasters in outputs/nfi8/, and nothing
else stops a stale one from shipping. Before anything is copied:
    - the lookup must be a plain renaming of models/train_reference_ftm.csv: the same rows in the
      same order, the published columns, the same values (numbers compared as numbers, so the
      writer's number formatting cannot fail or pass the check on its own);
    - each raster must be newer than the file its run wrote first (the trained legend, the
      reference table);
    - each table that describes the rasters must be newer than the newest raster;
    - the species table must list exactly the classes the product was trained on;
    - the lookup must number its plots 1..N and carry no blank biomass (donors without biomass
      are kept out of the reference library).
After the copy, every staged raster is decoded block by block, its band metadata is checked, and the
staged index must hold only integers the lookup can decode.

RUN   python pipeline/08d_stage_deposit.py --dry-run          (checks; writes nothing)
      NFIMAP_DEPOSIT_DOI=10.5281/zenodo.<n> NFIMAP_VIEWER_URL=<url> python pipeline/08d_stage_deposit.py
      The DOI is the version DOI the repository reserves for the new-version draft. The README
      template holds ⟨…⟩ placeholders for the viewer URL and for the numbers the validation stages
      produce; a real run refuses to stage until every one is filled in write_readme().
OUT   outputs/nfi8/_deposit_v2/   (10 rasters, 3 tables, README.md, checksums.sha256)
EXIT  0 = staged and verified (or, with --dry-run, every check passed), 1 = not ready.
"""
import csv
import hashlib
import io
import os
import re
import shutil
import sys

os.environ.setdefault('GDAL_PAM_ENABLED', 'NO')   # no .aux.xml sidecars: that write is the corruption vector

import numpy as np
import rasterio
from rasterio.enums import Resampling

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')
DEPOSIT = os.path.join(OUT, '_deposit_v2')

DOI = os.environ.get('NFIMAP_DEPOSIT_DOI', '')
DOI_V1 = '10.5281/zenodo.21701424'
VIEWER = os.environ.get('NFIMAP_VIEWER_URL', '⟨viewer URL⟩')
GITHUB = 'https://github.com/kim-seunguk/korea-forest-attribute-maps'
EE_COLLECTION = 'projects/ee-korea-forest-maps/assets/korea-forest-attributes-v2'
OVERVIEW_LEVELS = [2, 4, 8, 16, 32]

# ---------------------------------------------------------------------------------------------
# the released set. Anything not listed here is not deposited.
# ---------------------------------------------------------------------------------------------
# the three continuous attributes carry an ensemble-SD layer; dominant species carries `1 - vote
# margin` under a different suffix, so it is handled separately rather than folded into the loop.
CONTINUOUS = ['conifer_fraction', 'biomass', 'shannon']
VARIABLE = CONTINUOUS + ['dominant_species']
TIER = ['ftm']

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
    for _f in (f'imputation_{_t}.tif', f'nearest_reference_distance_{_t}.tif'):
        RENAME[_f] = _f
        SUBDIR[_f] = 'index'

# non-raster companions: (source under outputs/nfi8, deposit name)
LOOKUP = ('plot_lookup_ftm.csv', 'plot_lookup_ftm.csv')
SPECIES = ('dominant_species_table_ftm.csv', 'dominant_species_table.csv')
MANIFEST = ('data_records.csv', 'data_records.csv')
LEGEND_USED = os.path.join('m8_maps_30m', 'm8_domsp_legend_used.csv')   # written first by 05a_direct_prediction.R
REFERENCE = os.path.join('models', 'train_reference_ftm.csv')             # written first by 05b_plot_index_imputation.R

# published lookup column <- reference-library column, in the published order: the renaming 05d applies
LOOKUP_COLUMNS = [('plot_index', 'plot_idx'), ('nfi_plot_id', 'Plot_ID'),
                  ('coniferous_fraction', 'coniferfrac_full'), ('dominant_species', 'dom_species'),
                  ('aboveground_biomass_Mgha', 'Bio_a_ha'), ('shannon_diversity', 'shannon_full'),
                  ('forest_type', 'V_VG'), ('diameter_class', 'V_DC'), ('crown_density', 'V_DN'),
                  ('dominant_species_ba_share', 'dom_ba_share'), ('carbon_Mgha', 'Carbon_a_ha'),
                  ('basal_area_m2ha', 'BAPH')]

# Manifest notes that name a companion file by its working name.
NOTE_FIX = [('dominant_species_table_ftm.csv', 'dominant_species_table.csv')]


def read_csv(path):
    with io.open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, cols):
    """UTF-8 without a BOM. The working lookup carries one, and it breaks a naive header parse."""
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


def src_path(name):
    return os.path.join(OUT, SUBDIR[name], name)


def same_cell(a, b):
    """Two CSV cells hold the same value: equal as numbers when both parse, else equal as text."""
    a, b = a.strip(), b.strip()
    if a == b or {a, b} <= {'', 'NA'}:      # an empty cell and R's NA are the same missing value
        return True
    try:
        return float(a) == float(b)
    except ValueError:
        return False


# ---------------------------------------------------------------------------------------------
# currency: is every source from the current run?
# ---------------------------------------------------------------------------------------------
def check_currency():
    """Return a list of problems; empty means every source belongs to the current run."""
    need = [os.path.join(SUBDIR[n], n) for n in RENAME] + [LOOKUP[0], SPECIES[0], MANIFEST[0],
                                                          LEGEND_USED, REFERENCE]
    missing = [p for p in need if not os.path.exists(os.path.join(OUT, p))]
    if missing:
        return [f'missing {p}' for p in missing]

    def mt(rel):
        return os.path.getmtime(os.path.join(OUT, rel))

    problems = []
    for n in RENAME:
        first = REFERENCE if SUBDIR[n] == 'index' else LEGEND_USED
        if mt(os.path.join(SUBDIR[n], n)) < mt(first):
            problems.append(f'{n} is older than {first}: not written by the run that wrote {first}')
    newest = max(mt(os.path.join(SUBDIR[n], n)) for n in RENAME)
    for t in (MANIFEST[0], SPECIES[0]):
        if mt(t) < newest:
            problems.append(f'{t} is older than the newest raster: rebuild it after the rasters')
    if mt(LOOKUP[0]) < mt(REFERENCE):
        problems.append(f'{LOOKUP[0]} is older than {REFERENCE}')

    # The lookup is 05d's renaming of the reference library, rebuilt here from that library and
    # compared row by row in library order, so a lookup built from an older library, or edited by
    # hand, cannot pass. Cells are compared as parsed values: a number is a number however it was
    # formatted, and anything that does not parse is compared as text.
    rows = read_csv(os.path.join(OUT, LOOKUP[0]))
    ref = read_csv(os.path.join(OUT, REFERENCE))
    want = [c for c, _ in LOOKUP_COLUMNS]
    if list(rows[0].keys()) != want:
        problems.append(f'{LOOKUP[0]}: columns {list(rows[0].keys())} are not the published {want}')
    elif len(rows) != len(ref):
        problems.append(f'{LOOKUP[0]} has {len(rows)} rows, {REFERENCE} has {len(ref)}')
    else:
        differ = sum(1 for a, b in zip(rows, ref)
                     if any(not same_cell(a[c], b[s]) for c, s in LOOKUP_COLUMNS))
        if differ:
            problems.append(f'{LOOKUP[0]}: {differ} rows differ from {REFERENCE}: '
                            f'run 05d_build_plot_lookup.R again')

    idx = sorted(int(r['plot_index']) for r in rows)
    if idx != list(range(1, len(rows) + 1)):
        problems.append(f'{LOOKUP[0]}: plot_index is not 1..{len(rows)}')
    blank = sum(1 for r in rows if r['aboveground_biomass_Mgha'].strip() in ('', 'NA'))
    if blank:
        problems.append(f'{LOOKUP[0]}: {blank} plots without biomass (they are kept out of the library)')

    trained = {r['code'] for r in read_csv(os.path.join(OUT, LEGEND_USED))}
    listed = {r['code'] for r in read_csv(os.path.join(OUT, SPECIES[0]))}
    if trained != listed:
        problems.append(f'{SPECIES[0]} lists codes {sorted(listed - trained)} not trained and misses '
                        f'{sorted(trained - listed)}')
    return problems


# ---------------------------------------------------------------------------------------------
# species harmonisation
# ---------------------------------------------------------------------------------------------
def species_maps():
    """Korean field name -> (raster class code, scientific name).

    The lookup carries the plot's FIELD-measured dominant species (well over a hundred distinct names),
    not the raster classes, so it cannot be joined to the species table as shipped. The trained legend
    is the authoritative code<->name map; any species outside its named classes is what the raster calls
    class 38 (기타 / Other). The field name is KEPT, since collapsing it would discard exactly the detail
    that makes the plot-index family worth having, and the code is added alongside it.
    """
    legend = read_csv(os.path.join(OUT, LEGEND_USED))
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


def stage_lookup(code_of, sci_of, other, dry):
    """Copy the plot lookup with the species columns harmonised and the BOM dropped."""
    src_name, dst_name = LOOKUP
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
    print(f"  {src_name:30s} -> {dst_name:28s} {len(rows):>6,} rows, "
          f"{folded} folded into class {other} (Other)"
          + (f", {len(unresolved)} names without a scientific name" if unresolved else ""))
    return unresolved, len(rows)


def stage_species_table(dry):
    src_name, dst_name = SPECIES
    rows = read_csv(os.path.join(OUT, src_name))
    if not dry:
        write_csv(os.path.join(DEPOSIT, dst_name), rows, list(rows[0].keys()))
    print(f"  {src_name:30s} -> {dst_name:28s} {len(rows)} classes")
    return len(rows)


# ---------------------------------------------------------------------------------------------
# rasters
# ---------------------------------------------------------------------------------------------
def overview_resampling(name):
    """Same rule as 06d_finalize_rasters.py, which builds the overviews upstream: uncertainty and
    distance layers are continuous (checked first, since their names also contain the attribute), the
    species classes take the mode, and the plot index takes the nearest value (a mode of plot ids
    means nothing, an average invents plots)."""
    low = name.lower()
    if '_sd_' in low or 'uncertainty' in low or 'distance' in low:
        return Resampling.average, 'average'
    if 'dominant_species' in low:
        return Resampling.mode, 'mode'
    if low.startswith('imputation'):
        return Resampling.nearest, 'nearest'
    return Resampling.average, 'average'


def band_statistics(path):
    """min, max, mean, sd and valid share over every block, read-only; nodata = the declared value or NaN."""
    n = s1 = s2 = 0.0
    lo = hi = None
    with rasterio.open(path) as ds:
        nd = ds.nodata
        total = ds.width * ds.height
        for _, w in ds.block_windows(1):
            a = ds.read(1, window=w).astype('float64')
            m = np.isfinite(a)
            if nd is not None and np.isfinite(nd):
                m &= a != nd
            v = a[m]
            if v.size:
                n += v.size; s1 += float(v.sum()); s2 += float((v * v).sum())
                lo = float(v.min()) if lo is None else min(lo, float(v.min()))
                hi = float(v.max()) if hi is None else max(hi, float(v.max()))
    if n == 0:
        raise SystemExit(f"  {os.path.basename(path)} has no valid cell")
    mean = s1 / n
    sd = max(0.0, s2 / n - mean * mean) ** 0.5
    return dict(STATISTICS_MINIMUM=repr(lo), STATISTICS_MAXIMUM=repr(hi), STATISTICS_MEAN=repr(mean),
                STATISTICS_STDDEV=repr(sd), STATISTICS_VALID_PERCENT=f"{100.0 * n / total:.3f}")


def stage_raster(src_name, dry, meta):
    """Copy one raster, then fix its metadata on the COPY: true band statistics in place of the -9999
    placeholders the writers leave, the variable as band description and the unit as band unit, and
    overviews where the master has none. Every write is to the copy, which step [6] then decodes in full."""
    dst_name = RENAME[src_name]
    src = src_path(src_name)
    dst = os.path.join(DEPOSIT, dst_name)
    if dry:
        print(f"  {src_name:42s} -> {dst_name}")
        return dst

    shutil.copy2(src, dst)
    stats = band_statistics(dst)
    row = meta[dst_name]
    with rasterio.open(dst, 'r+') as ds:
        ds.update_tags(1, **stats)
        ds.set_band_description(1, row['variable'])
        ds.set_band_unit(1, row['unit'])
        n_ov = len(ds.overviews(1))
        tag = 'kept'
        if n_ov == 0:
            rs, tag = overview_resampling(dst_name)
            ds.build_overviews(OVERVIEW_LEVELS, rs)
            ds.update_tags(ns='rio_overview', resampling=tag)
    with rasterio.open(dst) as ds:
        n_ov = len(ds.overviews(1))
    print(f"  {src_name:42s} -> {dst_name:44s} ({n_ov} overview levels, {tag}; mean {float(stats['STATISTICS_MEAN']):.4g}, "
          f"valid {stats['STATISTICS_VALID_PERCENT']}%)")
    return dst


def check_metadata(path, meta):
    """The staged copy must carry a description, a unit, true statistics and no PAM sidecar."""
    problems = []
    row = meta[os.path.basename(path)]
    with rasterio.open(path) as ds:
        tags = ds.tags(1)
        if (ds.descriptions[0] or '') != row['variable']:
            problems.append(f"band description {ds.descriptions[0]!r} != {row['variable']!r}")
        if (ds.units[0] or '') != row['unit']:
            problems.append(f"band unit {ds.units[0]!r} != {row['unit']!r}")
        for k in ('STATISTICS_MEAN', 'STATISTICS_STDDEV'):
            if tags.get(k, '-9999') == '-9999':
                problems.append(f"{k} missing or the -9999 placeholder")
        if ds.nodata is None and ds.dtypes[0] != 'float32':
            problems.append("no nodata declared")
    if os.path.exists(path + '.aux.xml'):
        problems.append(".aux.xml sidecar present")
    return problems


def check_index(path, n_rows):
    """Every valid cell of the staged index must be an integer the lookup decodes (1..n_rows)."""
    lo, hi, nonint, seen = None, None, 0, set()
    with rasterio.open(path) as ds:
        nd = ds.nodata
        for _, w in ds.block_windows(1):
            a = ds.read(1, window=w)
            v = a[np.isfinite(a)]
            if nd is not None and np.isfinite(nd):
                v = v[v != nd]
            if v.size == 0:
                continue
            nonint += int(np.count_nonzero(v != np.round(v)))
            u = np.unique(v)
            seen.update(u.tolist())
            lo = u[0] if lo is None else min(lo, u[0])
            hi = u[-1] if hi is None else max(hi, u[-1])
    problems = []
    if lo is None:
        return ['the index raster has no valid cell']
    if nonint:
        problems.append(f'{nonint:,} cells hold a non-integer index')
    if lo < 1 or hi > n_rows:
        problems.append(f'index values {lo:g}..{hi:g} fall outside the lookup 1..{n_rows}')
    print(f"  index decodes: {len(seen):,} distinct plots of {n_rows:,}, values {lo:g}..{hi:g}")
    return problems


def stage_manifest(dry):
    """Copy data_records.csv with the companion files named as the deposit names them."""
    src = os.path.join(OUT, MANIFEST[0])
    rows = read_csv(src)
    if len(rows) != len(RENAME):
        raise SystemExit(f"  {MANIFEST[0]} has {len(rows)} rows, expected {len(RENAME)}")
    for r in rows:
        if r['dtype'] == 'MISSING':
            raise SystemExit(f"  manifest row {r['file']} is MISSING: the product is incomplete")
        if r['epsg'] != '5179':
            raise SystemExit(f"  manifest row {r['file']} has epsg={r['epsg']}, expected 5179")
        r['file'] = RENAME[r['file']]
        for old, new in NOTE_FIX:
            r['notes'] = r['notes'].replace(old, new)
    if not dry:
        write_csv(os.path.join(DEPOSIT, MANIFEST[1]), rows, list(rows[0].keys()))
    print(f"  {MANIFEST[0]:30s} -> {MANIFEST[1]:28s} {len(rows)} rows")
    return rows


def raster_format(paths):
    """The Format line of the README, read from the files instead of asserted."""
    blocks, comps = set(), set()
    for p in paths:
        with rasterio.open(p) as ds:
            blocks.add(f"{ds.block_shapes[0][1]} × {ds.block_shapes[0][0]}")
            comps.add(str(ds.compression.value if ds.compression else 'no'))
    return f"Tiled GeoTIFF, {' / '.join(sorted(blocks))} blocks, {' / '.join(sorted(comps))} compression"


# ---------------------------------------------------------------------------------------------
def numbers():
    """The accuracy figures the README quotes, read from the validation tables, never typed.
    Keys match the ⟨…⟩ placeholders in the template; a placeholder with no key here stays and blocks a real run."""
    def pick(rows, val, **kw):
        hit = [r for r in rows if all(str(r.get(k)) == str(v) for k, v in kw.items())]
        if len(hit) != 1:
            raise SystemExit(f"  numbers(): {len(hit)} rows for {kw}")
        return float(hit[0][val])
    m4a = read_csv(os.path.join(OUT, 'm4a_ladder_ftmonly.csv'))      # 07a_predictor_ladder.R, FTM_ONLY
    m11 = read_csv(os.path.join(OUT, 'm11_support_scale_ftm.csv'))   # 07o_support_scale.R, FTM_ONLY
    n = {}
    for tg, key in (('coniferfrac_full', 'conifer'), ('domsp_full', 'species'), ('Bio_a_ha', 'biomass'), ('shannon_full', 'shannon')):
        for sc, lab in (('spatial', 'cv'), ('test', 'test')):
            n[f'{lab} {key}'] = f"{pick(m4a, 'point', target=tg, predset='full', scheme=sc):.2f}"
    n['nir'] = f"{pick(m4a, 'nir', target='domsp_full', predset='full', scheme='spatial'):.2f}"
    for sup, lab in (('0', 'l30'), ('1', 'l1k'), ('10', 'l10k')):
        for tg, key in (('coniferfrac_full', 'conifer'), ('Bio_a_ha', 'biomass'), ('shannon_full', 'shannon')):
            n[f'{lab} {key}'] = f"{pick(m11, 'r2', target=tg, predset='full', support_km=sup, min_plots='1'):.2f}"
    return n


def write_readme(manifest, n_plots, n_classes, fmt, n_files, dry):
    m0 = manifest[0]
    doi = DOI or '⟨version 2 DOI⟩'
    lines = f"""# Nationwide 30 m maps of forest composition, biomass, and diversity in South Korea (2021–2025)

Kim, S., Shin, J. H., Han, H., and Choe, H. (2026).
Department of Agriculture, Forestry and Bioresources, Seoul National University.

Version 2 · DOI <https://doi.org/{doi}> · Licence **CC BY 4.0** · Code (MIT) <{GITHUB}>

Read `data_records.csv` alongside this file: it is the machine-readable schema for every raster (data
type, nodata, unit, valid range, grid, CRS).

## What changed from version 1

Version 1 (<https://doi.org/{DOI_V1}>) stays available. Version 2 replaces every file.

- **One predictor set.** Version 1 delivered every product twice, without and with the national
  forest-type map (`base`, `ftm`). Version 2 delivers the `ftm` set only.
- **{n_classes} dominant-species classes.** Two species fell below 20 reference plots and joined class 38
  (Other). Their codes, 7 and 17, are unused; every other code keeps its version 1 meaning.
- **Corrected plot positions.** For 191 reference plots first surveyed in 2025, the source coordinates
  disagreed with the plot number by 0.2–1.1 km. Version 2 takes their coordinates from the same
  survey's stand records, which agree with it.
- **Rebuilt inputs.** The Sentinel-2 composites are averaged into each 30 m cell of the product grid
  (version 1 resampled them twice), elevation is resampled once from the native SRTM grid, and the
  climate grids are aligned on cell centres with a coastal fill.
- **Validation by inventory cluster.** Cross-validation blocks, folds and bootstrap intervals treat the
  subplots of one inventory cluster as a single unit.
- ⚠ **Plot indices are not comparable between versions.** The reference library changed, so an index
  value names a different plot in each version. Decode each version's raster with its own lookup.
- **A recommended minimum support per layer.** `data_records.csv` gained the column
  `recommended_min_support_km`: 0.03 (the cell) for coniferous fraction, dominant species and the plot
  index; 1 (one kilometre) for aboveground biomass and tree-species diversity, whose cell-level accuracy
  is low (see "Accuracy, in brief").

## Looking at the maps without downloading them

The maps are published as a public Earth Engine image collection, with a browser viewer over it, so
they can be examined before or instead of downloading several gigabytes. The viewer shows version 2 of
the maps; the version 1 collection stays online under its own name.

| | |
|---|---|
| Viewer | <{VIEWER}> |
| Collection | `{EE_COLLECTION}` |

The collection also carries the four attributes of the plot-index family as ready-made rasters,
derived there from the plot index and the lookup table in this record.

## What is here

| | |
|---|---|
| Attributes | coniferous fraction · dominant tree species ({n_classes} classes) · aboveground biomass · tree-species diversity (Shannon *H*′) |
| Predictors | Sentinel-2 seasonal composites (2021–2025), climate, topography, and the national forest-type map (2025) |
| Product families | **direct prediction** (a random-forest model per attribute) · **plot-index imputation** (each cell is assigned one matched inventory plot) |
| Reference data | 8th Korean National Forest Inventory, 2021–2025; {n_plots:,} plots in the imputation reference library |
| Grid | EPSG:5179, {m0['width']} × {m0['height']} cells of {m0['pixel_m']} × {m0['pixel_y_m']} m, origin ({m0['origin_x']}, {m0['origin_y']}) |
| Extent | stocked forest in the 2025 national forest-type map, nationwide. The strip along the DMZ, which the map does not cover, is nodata, and cells that the national land-cover map (2024 imagery) shows as built-up, water or barren are removed |

**The two families answer different questions.** Direct prediction is more accurate attribute by
attribute. Plot-index imputation is less accurate per attribute, but every cell carries a combination
of attributes that was actually observed together in one field plot, which matters if you need the
attributes to be mutually consistent rather than individually optimal.

## Files

**Direct prediction**: `ftm_<attribute>_30m.tif` for `conifer_fraction`, `dominant_species`, `biomass`
and `shannon`, each paired with a per-cell relative uncertainty layer:

- `ftm_<attribute>_sd_30m.tif`: ensemble standard deviation, for the three continuous attributes. Higher = more uncertain.
- `ftm_dominant_species_uncertainty_30m.tif`: `1 − vote margin`. Higher = more uncertain.

**Plot-index imputation**

- `imputation_ftm.tif`: index of the matched inventory plot.
- `plot_lookup_ftm.csv`: join `plot_index` to get that plot's attributes.
- `nearest_reference_distance_ftm.tif`: random-forest **distance** from the cell to its matched plot,
  defined as `1 − the proportion of trees placing the two in the same terminal node`. Higher means a
  **less** similar match. It is a continuous diagnostic, not a thresholded mask.

**Tables**

- `dominant_species_table.csv`: the {n_classes} class codes with Korean and scientific names, prevalence,
  and per-class producer's and user's accuracy and F1. Accuracy cells are empty for classes the
  classifier never assigned: user's accuracy is undefined there, not zero.
- `data_records.csv`: per-file schema for the {len(manifest)} rasters (data type, nodata, unit, valid
  range, grid, recommended minimum support). The three tables are described in this README only.

## Decoding the dominant-species map

Raster values are class codes 0–38, with 7 and 17 unused. Join to `dominant_species_table.csv` on
`code`. Class 38 (기타 / Other) holds species too rare to model individually. Nodata is 255.

## Decoding the plot-index map

Raster values are plot indices 1…{n_plots:,} stored as float32. Join to `plot_lookup_ftm.csv` on
`plot_index`.

⚠ **Never resample this raster with an interpolating method.** Bilinear or cubic averages two plot
indices into a third plot that was never matched. Use nearest neighbour or a majority rule.

| column | content |
|---|---|
| `plot_index` | the raster value |
| `nfi_plot_id` | official inventory plot number |
| `dominant_species_code` | the class the dominant-species map would give this plot (38 = Other) |
| `dominant_species_kr`, `dominant_species_sci` | the plot's field-recorded dominant species, finer than the mapped classes |
| `coniferous_fraction` | conifer share of basal area per hectare, 0–1 |
| `aboveground_biomass_Mgha` | aboveground biomass, Mg ha⁻¹ |
| `shannon_diversity` | Shannon *H*′ of the species' shares of basal area per hectare |
| `forest_type` | forest type: 10 other conifers or conifer mixture · 11 *Pinus thunbergii* · 12 *Larix kaempferi* · 13 *Pinus rigida* · 14 *Pinus densiflora* · 15 *Pinus koraiensis* · 16 *Chamaecyparis obtusa* · 20 other oaks or oak mixture · 21 *Quercus variabilis* · 22 *Quercus acutissima* · 23 *Quercus mongolica* · 30 *Castanea crenata* · 40 evergreen broadleaves · 50 other broadleaves or broadleaf mixture · 60 mixed conifer–broadleaf |
| `diameter_class` | representative breast-height diameter in cm: 3 (below 6 cm), 12 (6–18 cm), 24 (18–30 cm), 36 (30 cm and above) |
| `crown_density` | representative canopy cover in %: 25 (50 % or less), 60 (51–75 %), 85 (76 % or more) |
| `dominant_species_ba_share` | the dominant species' share of basal area per hectare, 0–1 |
| `carbon_Mgha` | aboveground carbon, Mg ha⁻¹ |
| `basal_area_m2ha` | basal area, m² ha⁻¹ |

`forest_type`, `diameter_class` and `crown_density` are the plot's field records, not values read from
the forest-type map, and they are what the imputation matches on. `diameter_class = 36` is a diameter,
not the 36th class.

## Plot locations

Exact plot coordinates are not distributed. The lookup identifies each plot by its inventory number,
which encodes the approximate position of the plot cluster. The field measurements of the 8th
inventory are to be released by the Korea Forest Service at <https://kfss.forest.go.kr/stat/>, which
already distributes the microdata of earlier inventory cycles.

## Accuracy, in brief

Direct predictions, spatial-block cross-validation (50 km blocks, *K* = 10, inventory clusters kept
whole) and a spatially disjoint test set:

| attribute | cross-validation | disjoint test |
|---|---|---|
| coniferous fraction | *R*² ⟨cv conifer⟩ | ⟨test conifer⟩ |
| dominant species | accuracy ⟨cv species⟩ (majority-class baseline ⟨nir⟩) | ⟨test species⟩ |
| aboveground biomass | *R*² ⟨cv biomass⟩ | ⟨test biomass⟩ |
| tree-species diversity | *R*² ⟨cv shannon⟩ | ⟨test shannon⟩ |

**The four attributes are not equally reliable at a single cell.** Coniferous fraction and dominant
species can be read cell by cell at the accuracy above. Aboveground biomass and tree-species diversity
should be aggregated to at least 1 km and read as gradients; accuracy rises with aggregation (*R*²,
cross-validation):

| aggregation | coniferous fraction | biomass | Shannon |
|---|---|---|---|
| single cell (30 m) | ⟨l30 conifer⟩ | ⟨l30 biomass⟩ | ⟨l30 shannon⟩ |
| 1 km | ⟨l1k conifer⟩ | ⟨l1k biomass⟩ | ⟨l1k shannon⟩ |
| 10 km | ⟨l10k conifer⟩ | ⟨l10k biomass⟩ | ⟨l10k shannon⟩ |

For biomass and diversity the plot-index family falls below a mean-predictor baseline, so the direct
maps are preferable for those two attributes. Per-class species accuracy is concentrated in the more
prevalent species; check `dominant_species_table.csv` before relying on any one class. The uncertainty
layers' intervals hold their nominal coverage on the held-out plots as a whole; they rank cells by
reliability and are not cell-level error bounds. The nearest-reference distance layer was not
evaluated against observed error.

## Using the maps

- **Aggregating the plot index:** do not average the matched plots' attributes inside a unit. That
  shrinks their variance, which is the property the index exists to keep; take the modal plot or one
  donor instead.
- **Aggregating diversity:** the mean of predicted Shannon over an area is the mean local (alpha)
  diversity, not the diversity of the area. For the second, pool the matched plots' stem lists through
  the index.
- **Mask vintage:** the extent comes from the 2025 forest-type map and the 2024 land-cover map. Before
  studying a period that ended in a major disturbance, check the mask has not removed the area of
  interest: a recent wildfire scar is mapped as non-stocked.
- **Do not let a GIS write to these files.** Opening them is safe; computing statistics or building
  pyramids writes into the TIFF, and such a write has corrupted rasters of this dataset before. Work
  on copies.

## Format

{fmt}, internal overviews. Tables are UTF-8 CSV without a byte-order mark.

## Citation

Kim, S., Shin, J. H., Han, H., and Choe, H.: Nationwide 30 m maps of forest composition, biomass, and
diversity in South Korea (2021–2025) from direct prediction and plot-index imputation, version 2
[Data set], Zenodo, https://doi.org/{doi}, 2026.

---

Every raster in this record was read block by block in full before upload, and the plot-index raster
was checked against its lookup. `checksums.sha256` covers all {n_files} files.
"""
    n = numbers()
    lines = re.sub(r'⟨([^⟩]+)⟩', lambda m: n.get(m.group(1), m.group(0)), lines)
    holes = sorted(set(re.findall(r'⟨[^⟩]*⟩', lines)))
    if holes and not dry:
        raise SystemExit(f"  README still has {len(holes)} placeholder(s): {holes}")
    if not dry:
        with io.open(os.path.join(DEPOSIT, 'README.md'), 'w', encoding='utf-8', newline='\n') as f:
            f.write(lines)
    print(f"  {'(generated)':30s} -> {'README.md':28s} {len(lines.splitlines())} lines"
          + (f"; {len(holes)} placeholder(s) to fill: {', '.join(holes)}" if holes else ""))


def main(argv):
    dry = '--dry-run' in argv
    print(f"staging deposit -> {os.path.relpath(DEPOSIT, ROOT)}" + ("  [DRY RUN]" if dry else ""))
    if not dry and not DOI:
        print("  set NFIMAP_DEPOSIT_DOI to the version DOI reserved for the new-version draft")
        return 1

    print("\n[1] integrity gate on the working tree")
    # Loaded by path: the stage number makes 08a_verify_raster_integrity an invalid module name,
    # so a plain import of it cannot work.
    import importlib.util
    gate = os.path.join(HERE, '08a_verify_raster_integrity.py')
    spec = importlib.util.spec_from_file_location('raster_integrity_gate', gate)
    vri = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vri)
    if vri.main() != 0:
        print("  STOP: the working tree has damaged rasters.")
        return 1

    print("\n[2] currency of the sources")
    problems = check_currency()
    for p in problems:
        print(f"  ! {p}")
    if problems:
        print(f"  STOP: {len(problems)} source(s) are not from the current run.")
        return 1
    print("  every source belongs to the current run")

    if not dry:
        if os.path.exists(DEPOSIT):
            print(f"\n  {os.path.relpath(DEPOSIT, ROOT)} already exists: remove it first")
            return 1
        os.makedirs(DEPOSIT)

    print("\n[3] manifest (read first: it supplies each raster's variable and unit)")
    manifest = stage_manifest(dry)
    meta = {r['file']: r for r in manifest}

    print("\n[4] rasters")
    staged = [stage_raster(n, dry, meta) for n in RENAME]

    print("\n[5] tables and README")
    code_of, sci_of, other = species_maps()
    unresolved, n_plots = stage_lookup(code_of, sci_of, other, dry)
    if unresolved:
        print(f"  ! {len(unresolved)} species names have no scientific name: {sorted(unresolved)[:8]}")
    n_classes = stage_species_table(dry)
    n_files = len(RENAME) + 3 + 1                 # rasters + lookup, species table, data_records + README
    fmt = raster_format(staged if not dry else [src_path(n) for n in RENAME])
    write_readme(manifest, n_plots, n_classes, fmt, n_files, dry)

    if dry:
        print(f"\nDRY RUN: every check passed; {n_files} files would be staged. Nothing was written.")
        return 0

    print("\n[6] integrity gate on the staged copy (tile tables, every block decoded, metadata)")
    bad = 0
    for p in staged:
        r = vri.scan(p)
        if r is None:
            print(f"  UNPARSED {os.path.basename(p)}")
            bad += 1
            continue
        n_tiles, n_bad, n_unw, _ = r
        failed = vri.decode(p) if not (n_bad or n_unw) else 0
        if n_bad or n_unw or failed:
            print(f"  DAMAGED {os.path.basename(p)}  {n_bad} bad + {n_unw} unwritten of {n_tiles} tiles, "
                  f"{failed} blocks failed to decode")
            bad += 1
        for q in check_metadata(p, meta):
            print(f"  ! {os.path.basename(p)}: {q}")
            bad += 1
    for p in staged:
        if os.path.basename(p).startswith('imputation_'):
            for q in check_index(p, n_plots):
                print(f"  ! {q}")
                bad += 1
    if bad:
        print(f"\n{bad} problem(s) in the staged copy. The deposit is NOT ready.")
        return 1
    print(f"  all {len(staged)} staged rasters clean and fully decoded")

    print("\n[7] checksums")
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
