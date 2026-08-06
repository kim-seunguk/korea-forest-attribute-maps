"""08c_build_data_records.py — the machine-readable file manifest that ships with the dataset.

Opens every released raster and emits one row per file describing it: data type, nodata value,
unit, valid range, pixel size, grid origin, dimensions and CRS. Reading this from the files
themselves means the published table cannot drift from what was actually written, which
hand-typing it could not guarantee.

Categorical layers name their legend file, and uncertainty layers state which direction is worse.

The filenames in this table are the ones the pipeline writes, which are also the ones the dataset
is published under; nothing is renamed downstream.

OUTPUT outputs/nfi8/data_records.csv
RUN    python pipeline/08c_build_data_records.py
"""
import os, csv, io
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')

# (relative path, family, predictor set, variable, unit, valid range, notes/legend/direction)
FILES = []
for tier, ps in (('base', 'base'), ('ftm', 'FTM-added')):
    md = 'm8_maps_30m'; un = 'm8_uncertainty_30m'
    FILES += [
        (f'{md}/{tier}_conifer_fraction_30m.tif', 'direct', ps, 'coniferous fraction', '0–1 fraction', '0–1', ''),
        (f'{md}/{tier}_dominant_species_30m.tif', 'direct', ps, 'dominant species', 'class code', '0–38', 'class codes, species names and per-class accuracy in dominant_species_table.csv'),
        (f'{md}/{tier}_biomass_30m.tif', 'direct', ps, 'aboveground biomass', 'Mg ha-1', '>=0', ''),
        (f'{md}/{tier}_shannon_30m.tif', 'direct', ps, 'tree-species diversity (Shannon H′)', 'unitless', '>=0', ''),
        (f'{un}/{tier}_conifer_fraction_sd_30m.tif', 'direct uncertainty', ps, 'coniferous fraction SD', '0–1 fraction', '>=0', 'ensemble SD; higher = more uncertain'),
        (f'{un}/{tier}_biomass_sd_30m.tif', 'direct uncertainty', ps, 'biomass SD', 'Mg ha-1', '>=0', 'ensemble SD; higher = more uncertain'),
        (f'{un}/{tier}_shannon_sd_30m.tif', 'direct uncertainty', ps, 'Shannon SD', 'unitless', '>=0', 'ensemble SD; higher = more uncertain'),
        (f'{un}/{tier}_dominant_species_uncertainty_30m.tif', 'direct uncertainty', ps, 'dominant-species uncertainty', 'index', '0–1', '1 - vote margin; higher = more uncertain'),
    ]
FILES += [
    ('index/imputation_base.tif', 'plot-index', 'base', 'matched NFI plot index', 'plot index', '1..N', 'integer index; join to plot_lookup_base.csv for the matched plot attributes'),
    ('index/imputation_ftm.tif', 'plot-index', 'FTM-added', 'matched NFI plot index', 'plot index', '1..N', 'integer index; join to plot_lookup_ftm.csv for the matched plot attributes'),
    ('index/nearest_reference_distance_base.tif', 'plot-index nearest-reference distance', 'base', 'random-forest distance to the matched reference plot', 'index', '0–1', '1 - terminal-node proximity; higher = less similar match; pairs with imputation_base.tif'),
    ('index/nearest_reference_distance_ftm.tif', 'plot-index nearest-reference distance', 'FTM-added', 'random-forest distance to the matched reference plot', 'index', '0–1', '1 - terminal-node proximity; higher = less similar match; pairs with imputation_ftm.tif'),
]

rows = []
for rel, fam, ps, var, unit, vrange, note in FILES:
    p = os.path.join(OUT, rel)
    if not os.path.exists(p):
        rows.append(dict(file=os.path.basename(rel), family=fam, predictor_set=ps, variable=var,
                         dtype='MISSING', nodata='', unit=unit, valid_range=vrange, pixel_m='',
                         width='', height='', epsg='', origin_x='', origin_y='', notes=note))
        continue
    with rasterio.open(p) as s:
        t = s.transform
        rows.append(dict(file=os.path.basename(rel), family=fam, predictor_set=ps, variable=var,
                         dtype=s.dtypes[0], nodata=('' if s.nodata is None else s.nodata),
                         unit=unit, valid_range=vrange, pixel_m=int(round(t.a)),
                         width=s.width, height=s.height, epsg=(s.crs.to_epsg() if s.crs else ''),
                         origin_x=round(t.c, 1), origin_y=round(t.f, 1), notes=note))

dst = os.path.join(OUT, 'data_records.csv')
cols = ['file', 'family', 'predictor_set', 'variable', 'dtype', 'nodata', 'unit', 'valid_range',
        'pixel_m', 'width', 'height', 'epsg', 'origin_x', 'origin_y', 'notes']
with io.open(dst, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in rows: w.writerow(r)
print('wrote', os.path.relpath(dst, ROOT), f'({len(rows)} files)')
for r in rows:
    print(f"  {r['file']:38s} {str(r['dtype']):8s} nodata={str(r['nodata']):>5s} {r['pixel_m']}m "
          f"{r['width']}x{r['height']} epsg={r['epsg']}  {r['family']}/{r['predictor_set']}")
