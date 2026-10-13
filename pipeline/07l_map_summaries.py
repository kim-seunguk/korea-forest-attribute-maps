"""07l_map_summaries.py — summary statistics of the delivered products.

Elevation-band means, dominant-species shares by elevation band, uncertainty by class of predicted
value, and how the plot-index map uses its reference library. The 200 m band edges, the value-class
widths and the treatment of nodata are set here.

All summaries are of the delivered rasters. Every layer sits on the
covariate grid with 256 x 256 tiles, so the pass reads one tile stripe at a time. Reading single rows
instead would decode each stripe once per row. Each layer is read through its declared nodata, so a
class map's nodata code and a float layer's NaN are both treated as absent.

⚠ These are descriptions of the mapped products. Elevation is itself a predictor, so an elevation
gradient in the maps is not independent evidence about forests.

INPUT   outputs/nfi8/m8_maps_30m/ftm_{conifer_fraction,biomass,shannon,dominant_species}_30m.tif
        outputs/nfi8/m8_uncertainty_30m/ftm_{conifer_fraction_sd,biomass_sd,shannon_sd,dominant_species_uncertainty}_30m.tif
        outputs/nfi8/index/imputation_ftm.tif · data/interim/covar_5179/T_ELE.tif
        data/species_codes.csv
OUTPUT  outputs/nfi8/m8_map_summary_elevation_ftm.csv    (band means + cell counts)
        outputs/nfi8/m8_map_summary_domsp_elev_ftm.csv   (dominant-species share per band)
        outputs/nfi8/m8_map_summary_uncertainty_ftm.csv  (uncertainty by class of predicted value)
        outputs/nfi8/m8_map_summary_plotuse_ftm.csv      (plot-index reference-library use)
RUN     python pipeline/07l_map_summaries.py
"""
import csv
import io
import os

import numpy as np
import rasterio

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')
MD, UN = os.path.join(OUT, 'm8_maps_30m'), os.path.join(OUT, 'm8_uncertainty_30m')
TIER = 'ftm'                                         # the delivered file prefix (R/paths.R, FILE_TIER)
SUF = '_' + TIER

CONT = {
    'coniferfrac': os.path.join(MD, f'{TIER}_conifer_fraction_30m.tif'),
    'biomass':     os.path.join(MD, f'{TIER}_biomass_30m.tif'),
    'shannon':     os.path.join(MD, f'{TIER}_shannon_30m.tif'),
    'coniferfrac_sd': os.path.join(UN, f'{TIER}_conifer_fraction_sd_30m.tif'),
    'biomass_sd':     os.path.join(UN, f'{TIER}_biomass_sd_30m.tif'),
    'shannon_sd':     os.path.join(UN, f'{TIER}_shannon_sd_30m.tif'),
    'domsp_uncertainty': os.path.join(UN, f'{TIER}_dominant_species_uncertainty_30m.tif'),
}
NAMES = list(CONT)
ELE_EDGES = [200, 400, 600, 800]                     # 200 m bands; digitize gives 5 classes
ELE_LAB = ['<200', '200-400', '400-600', '600-800', '>800']
STRIPE = 256

src = {k: rasterio.open(v) for k, v in CONT.items()}
domsp = rasterio.open(os.path.join(MD, f'{TIER}_dominant_species_30m.tif'))
index = rasterio.open(os.path.join(OUT, 'index', f'imputation_{TIER}.tif'))
elev = rasterio.open(os.path.join(ROOT, 'data', 'interim', 'covar_5179', 'T_ELE.tif'))
H, W = src['coniferfrac'].shape
for r in list(src.values()) + [domsp, index, elev]:
    assert r.shape == (H, W), 'grid mismatch: %s' % r.name

acc = {}                 # (band, layer) -> [sum, n]
domsp_elev = {}          # band -> {code: cells}
cf_class = {}            # coniferous-fraction decile -> [sum sd, n]
bio_class = {}           # 50 Mg/ha biomass class -> [sum sd, n]
plot_use = {}            # plot index -> cells


def add(key, v):
    s = acc.setdefault(key, [0.0, 0])
    g = np.isfinite(v)
    s[0] += float(v[g].sum())
    s[1] += int(g.sum())


for r0 in range(0, H, STRIPE):
    win = rasterio.windows.Window(0, r0, W, min(STRIPE, H - r0))
    cf = src['coniferfrac'].read(1, window=win, masked=True).filled(np.nan)   # declared nodata -> NaN
    ok = np.isfinite(cf)                             # the forest footprint, identical in every layer
    if not ok.any():
        continue
    band = np.digitize(elev.read(1, window=win), ELE_EDGES)
    vals = {}
    for name in NAMES:
        v = src[name].read(1, window=win, masked=True).astype('float32').filled(np.nan)
        v[~ok] = np.nan
        vals[name] = v
        for b in range(len(ELE_LAB)):
            m = ok & (band == b)
            if m.any():
                add((b, name), v[m])
    # uncertainty against the predicted value: coniferous-fraction deciles, 50 Mg/ha biomass classes
    for key, q, store in ((None, np.clip(np.nan_to_num(cf * 10, nan=-1), -1, 9).astype('int16'), cf_class),
                          (None, np.clip(np.nan_to_num(vals['biomass'] / 50, nan=-1), -1, 5).astype('int16'), bio_class)):
        sd = vals['coniferfrac_sd'] if store is cf_class else vals['biomass_sd']
        for b in range(10 if store is cf_class else 6):
            m = ok & (q == b)
            if m.any():
                s = store.setdefault(b, [0.0, 0])
                s[0] += float(np.nansum(sd[m]))
                s[1] += int(m.sum())
    ds = domsp.read(1, window=win)
    for b in range(len(ELE_LAB)):
        m = ok & (band == b) & (ds != 255)
        if m.any():
            c = domsp_elev.setdefault(b, {})
            u, n = np.unique(ds[m], return_counts=True)
            for uu, nn in zip(u.tolist(), n.tolist()):
                c[uu] = c.get(uu, 0) + nn
    ix = index.read(1, window=win, masked=True).astype('float64').filled(np.nan)
    mi = np.isfinite(ix)
    if mi.any():
        u, n = np.unique(ix[mi].astype('int32'), return_counts=True)
        for uu, nn in zip(u.tolist(), n.tolist()):
            plot_use[uu] = plot_use.get(uu, 0) + nn
    if (r0 // STRIPE) % 20 == 0:
        print('  stripe %d/%d' % (r0 // STRIPE, H // STRIPE), flush=True)


def w(name, header, rows):
    name = name.replace('.csv', SUF + '.csv')
    p = os.path.join(OUT, name)
    with io.open(p, 'w', encoding='utf-8', newline='') as f:
        cw = csv.writer(f)
        cw.writerow(header)
        cw.writerows(rows)
    print('wrote %s (%d rows)' % (name, len(rows)))


def mean(b, n):
    s = acc.get((b, n))
    return round(s[0] / s[1], 4) if s and s[1] else ''


w('m8_map_summary_elevation.csv', ['elevation_band_m'] + NAMES + ['cells'],
  [[ELE_LAB[b]] + [mean(b, n) for n in NAMES] + [acc.get((b, 'coniferfrac'), [0, 0])[1]]
   for b in range(len(ELE_LAB))])

leg = {int(r['code']): r['sci_name'] for r in csv.DictReader(
    io.open(os.path.join(ROOT, 'data', 'species_codes.csv'),
            encoding='utf-8-sig'))}
rows = []
for b in range(len(ELE_LAB)):
    c = domsp_elev.get(b, {})
    tot = sum(c.values()) or 1
    for code, n in sorted(c.items(), key=lambda x: -x[1]):
        rows.append([ELE_LAB[b], code, leg.get(code, ''), n, round(100.0 * n / tot, 2)])
w('m8_map_summary_domsp_elev.csv', ['elevation_band_m', 'code', 'sci_name', 'cells', 'share_pct'], rows)

rows = [['coniferfrac', '%.1f-%.1f' % (b / 10, (b + 1) / 10), round(s[0] / s[1], 4), s[1]]
        for b, s in sorted(cf_class.items())]
# the top class is open-ended: clip() folds every value of 250 Mg/ha and above into it
rows += [['biomass_Mgha', ('%d-%d' % (b * 50, (b + 1) * 50)) if b < 5 else '>=250', round(s[0] / s[1], 2), s[1]]
         for b, s in sorted(bio_class.items())]
w('m8_map_summary_uncertainty.csv', ['attribute', 'predicted_value_class', 'mean_uncertainty', 'cells'], rows)

sizes = np.array(sorted(plot_use.values()))
tot = int(sizes.sum())
w('m8_map_summary_plotuse.csv', ['statistic', 'value'], [
    ['plots_assigned_at_least_one_cell', len(plot_use)],
    ['mapped_cells', tot],
    ['cells_per_plot_min', int(sizes.min())],
    ['cells_per_plot_median', int(np.median(sizes))],
    ['cells_per_plot_q25', int(np.percentile(sizes, 25))],
    ['cells_per_plot_q75', int(np.percentile(sizes, 75))],
    ['cells_per_plot_max', int(sizes.max())],
    ['ha_per_plot_min', round(float(sizes.min()) * 0.09, 2)],
    ['ha_per_plot_max', round(float(sizes.max()) * 0.09, 1)],
    ['ha_per_plot_median', round(float(np.median(sizes)) * 0.09, 1)],
    ['area_share_pct_of_most_used_tenth', round(100.0 * sizes[int(0.9 * len(sizes)):].sum() / tot, 1)],
])
print('\nplots used = %d | median %d cells (%.0f ha) | most-used tenth = %.1f %% of area'
      % (len(plot_use), np.median(sizes), np.median(sizes) * 0.09,
         100.0 * sizes[int(0.9 * len(sizes)):].sum() / tot))
