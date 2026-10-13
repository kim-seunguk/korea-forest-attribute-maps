"""08c_build_data_records.py — the machine-readable file manifest that ships with the dataset.

Opens every released raster and emits one row per file describing it: data type, nodata value,
unit, valid range, pixel size, grid origin, dimensions and CRS. Reading this from the files
themselves means the published table cannot drift from what was actually written, which
hand-typing it could not guarantee.

Categorical layers name their legend file and uncertainty layers state which direction is worse.

Pixel size and origin are written at full precision. The cell is 30.000454545 m wide and 30 m tall
and the origin is 746150.0931818181; rounded to 30 and 746150.1, a grid rebuilt from this table sits
8.45 m off at the eastern edge.

The filenames in this table are the ones the pipeline writes; 08d renames them to the published
names when it stages the deposit.

OUTPUT outputs/nfi8/data_records.csv
RUN    python pipeline/08c_build_data_records.py
"""
import os, csv, io
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')

MD, UN = 'm8_maps_30m', 'm8_uncertainty_30m'
# (relative path, family, predictor set, variable, unit, valid range, notes)
FILES = [
    (f'{MD}/ftm_conifer_fraction_30m.tif', 'direct', 'FTM-added', 'coniferous fraction', '0–1 fraction', '0–1', ''),
    (f'{MD}/ftm_dominant_species_30m.tif', 'direct', 'FTM-added', 'dominant species', 'class code', '0–36', 'class codes, species names and per-class accuracy in dominant_species_table.csv'),
    (f'{MD}/ftm_biomass_30m.tif', 'direct', 'FTM-added', 'aboveground biomass', 'Mg ha-1', '>=0', 'accuracy by level of aggregation in the README'),
    (f'{MD}/ftm_shannon_30m.tif', 'direct', 'FTM-added', 'tree-species diversity (Shannon H′)', 'unitless', '>=0', 'accuracy by level of aggregation in the README'),
    (f'{UN}/ftm_conifer_fraction_sd_30m.tif', 'direct uncertainty', 'FTM-added', 'coniferous fraction SD', '0–1 fraction', '>=0', 'ensemble SD; higher = more uncertain'),
    (f'{UN}/ftm_biomass_sd_30m.tif', 'direct uncertainty', 'FTM-added', 'biomass SD', 'Mg ha-1', '>=0', 'ensemble SD; higher = more uncertain'),
    (f'{UN}/ftm_shannon_sd_30m.tif', 'direct uncertainty', 'FTM-added', 'Shannon SD', 'unitless', '>=0', 'ensemble SD; higher = more uncertain'),
    (f'{UN}/ftm_dominant_species_uncertainty_30m.tif', 'direct uncertainty', 'FTM-added', 'dominant-species uncertainty', 'index', '0–1', '1 - vote margin; higher = more uncertain'),
    ('index/imputation_ftm.tif', 'plot-index', 'FTM-added', 'matched NFI plot index', 'plot index', '1..N', 'integer index; join to plot_lookup_ftm.csv for the matched plot attributes; resample nearest or modal only; do not read biomass or diversity off it'),
    ('index/nearest_reference_distance_ftm.tif', 'plot-index nearest-reference distance', 'FTM-added', 'random-forest distance to the matched reference plot', 'index', '0–1', '1 - terminal-node proximity; higher = less similar match; pairs with imputation_ftm.tif'),
]

rows = []
for rel, fam, ps, var, unit, vrange, note in FILES:
    p = os.path.join(OUT, rel)
    if not os.path.exists(p):
        rows.append(dict(file=os.path.basename(rel), family=fam, variable=var,
                         dtype='MISSING', nodata='', unit=unit, valid_range=vrange, pixel_m='', pixel_y_m='',
                         width='', height='', epsg='', origin_x='', origin_y='', notes=note))
        continue
    with rasterio.open(p) as s:
        t = s.transform
        rows.append(dict(file=os.path.basename(rel), family=fam, variable=var,
                         dtype=s.dtypes[0], nodata=('' if s.nodata is None else s.nodata),
                         unit=unit, valid_range=vrange, pixel_m=repr(t.a), pixel_y_m=repr(-t.e),
                         width=s.width, height=s.height, epsg=(s.crs.to_epsg() if s.crs else ''),
                         origin_x=repr(t.c), origin_y=repr(t.f), notes=note))

dst = os.path.join(OUT, 'data_records.csv')
cols = ['file', 'family', 'variable', 'dtype', 'nodata', 'unit', 'valid_range',
        'pixel_m', 'pixel_y_m', 'width', 'height', 'epsg', 'origin_x', 'origin_y', 'notes']
with io.open(dst, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in rows: w.writerow(r)
print('wrote', os.path.relpath(dst, ROOT), f'({len(rows)} files)')
for r in rows:
    print(f"  {r['file']:42s} {str(r['dtype']):8s} nodata={str(r['nodata']):>5s} {r['pixel_m']}m "
          f"{r['width']}x{r['height']} epsg={r['epsg']}  {r['family']}")
