"""07i_temporal_sensitivity.py — accuracy by the year a plot was measured.

The inventory visits a fifth of its plots each year, so the reference data spans 2021-2025 while the
predictors are one five-year Sentinel-2 composite and a single forest-type map edition. A plot
measured at the start of the period is therefore matched to imagery that partly postdates it. This
reports accuracy separately for each survey year.

This script reuses the per-plot spatial out-of-fold predictions that 07g writes and refits no model.
PREDSET picks the 07g run: ftm (the default) reads the delivered set's *_ftm files, and base
reads the unsuffixed files of the comparison run without the forest-type map.

INPUT  data/interim/nfi8/plot_year.csv               (Plot_ID, Survey_year)
       outputs/nfi8/uncertainty_points<SUF>.csv      (target, Plot_ID, pred, obs, ...)
OUTPUT outputs/nfi8/temporal_sensitivity<SUF>.csv    (target, survey_year, metric, value, n)
RUN    python pipeline/07i_temporal_sensitivity.py
"""
import os, csv, io
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')
PREDSET = os.environ.get('PREDSET', 'ftm')
if PREDSET not in ('base', 'ftm'):
    raise SystemExit('PREDSET must be "ftm" (the delivered set) or "base" (the run without the forest-type map)')
SUF = '_ftm' if PREDSET == 'ftm' else ''

year = {r['Plot_ID']: r['Survey_year']
        for r in csv.DictReader(io.open(os.path.join(ROOT, 'data', 'interim', 'nfi8', 'plot_year.csv'), encoding='utf-8-sig'))}
pts = list(csv.DictReader(io.open(os.path.join(OUT, f'uncertainty_points{SUF}.csv'), encoding='utf-8-sig')))

CONT = {'coniferfrac_full', 'Bio_a_ha', 'shannon_full'}           # R2 ; domsp -> accuracy
rows = []
def r2(o, p):
    o, p = np.asarray(o, float), np.asarray(p, float)
    ss = np.sum((o - o.mean()) ** 2)
    return np.nan if ss == 0 else 1 - np.sum((o - p) ** 2) / ss

for T in ['coniferfrac_full', 'Bio_a_ha', 'shannon_full', 'domsp_full']:
    sub = [x for x in pts if x['target'] == T]
    byyr = {}
    for x in sub:
        y = year.get(x['Plot_ID'])
        if y is None: continue
        byyr.setdefault(y, []).append(x)
    for y in sorted(byyr):
        g = byyr[y]
        if T in CONT:
            o = [float(x['obs']) for x in g]; p = [float(x['pred']) for x in g]
            rows.append(dict(target=T, survey_year=y, metric='R2', value=round(r2(o, p), 4), n=len(g)))
        else:
            acc = np.mean([x['obs'] == x['pred'] for x in g])
            rows.append(dict(target=T, survey_year=y, metric='accuracy', value=round(float(acc), 4), n=len(g)))

dst = os.path.join(OUT, f'temporal_sensitivity{SUF}.csv')
with io.open(dst, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['target', 'survey_year', 'metric', 'value', 'n']); w.writeheader()
    for r in rows: w.writerow(r)
print('wrote', os.path.relpath(dst, ROOT))
for T in ['coniferfrac_full', 'Bio_a_ha', 'shannon_full', 'domsp_full']:
    tr = [r for r in rows if r['target'] == T]
    print(f"  {T:16s} " + "  ".join(f"{r['survey_year']}:{r['value']:+.3f}(n{r['n']})" for r in tr))
