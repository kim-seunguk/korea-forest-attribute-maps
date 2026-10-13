"""07h_uncertainty_coverage.py — how much of the truth the uncertainty layers actually cover.

The released uncertainty layer is the standard deviation across the random-forest ensemble. Whether
that quantity ranks locations sensibly is one question; whether a nominal 80 % interval built from
it contains 80 % of the observations is a different and stricter one. This measures the second.

Interval = pred +/- z * ensemble SD, treating the SD as Gaussian, z = 0.6745 / 1.2816 / 1.9600 for
50 / 80 / 95 %. Coverage is the fraction of held-out plots falling inside it. Mean interval width is
reported alongside, since a wide interval can reach nominal coverage without being informative.
Biomass is additionally split into thirds, twice: by the observed value, because coverage that holds
on average can still fail at the top of the range, and by the predicted value, because thirds of the
observed value build regression to the mean into the result (the outer thirds look under-covered
whatever the interval does) and a user only ever knows the prediction.

Reuses the per-plot spatially out-of-fold predictions 07g writes; no model is refitted.
NFIMAP_TIER picks the 07g run: ftm (the default) reads the delivered set's *_ftm files, base reads
the unsuffixed files of the comparison run without the forest-type map.

INPUT  outputs/nfi8/m7_uncertainty_points<SUF>.csv  (target, Plot_ID, pred, obs, uncertainty, err)
OUTPUT outputs/nfi8/m7_coverage<SUF>.csv            (target, level, nominal, coverage, mean_width, n, stratum)
RUN    python pipeline/07h_uncertainty_coverage.py
"""
import os, csv, io
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
TIER = os.environ.get('NFIMAP_TIER', 'ftm')
if TIER not in ('base', 'ftm'):
    raise SystemExit('NFIMAP_TIER must be "ftm" (the delivered set) or "base" (the run without the forest-type map)')
SUF = '_ftm' if TIER == 'ftm' else ''
SRC = os.path.join(ROOT, 'outputs', 'nfi8', f'm7_uncertainty_points{SUF}.csv')
DST = os.path.join(ROOT, 'outputs', 'nfi8', f'm7_coverage{SUF}.csv')
CONT = ['coniferfrac_full', 'Bio_a_ha', 'shannon_full']            # continuous only; domsp uses vote margin
Z = {50: 0.6744898, 80: 1.2815516, 95: 1.9599640}                 # two-sided normal quantiles

rows = list(csv.DictReader(io.open(SRC, encoding='utf-8-sig')))
by = {t: {'pred': [], 'obs': [], 'sd': []} for t in CONT}
for r in rows:
    t = r['target']
    if t not in by: continue
    try:
        p, o, u = float(r['pred']), float(r['obs']), float(r['uncertainty'])
    except (ValueError, TypeError):
        continue
    if not (np.isfinite(p) and np.isfinite(o) and np.isfinite(u)): continue
    by[t]['pred'].append(p); by[t]['obs'].append(o); by[t]['sd'].append(u)

out = []
def cover(pred, obs, sd, level, stratum='all'):
    z = Z[level]; ae = np.abs(obs - pred)
    return dict(target=T, level=f'{level}%', nominal=level / 100,
                coverage=round(float(np.mean(ae <= z * sd)), 4),
                mean_width=round(float(np.mean(2 * z * sd)), 4),
                n=len(pred), stratum=stratum)

for T in CONT:
    pred = np.array(by[T]['pred']); obs = np.array(by[T]['obs']); sd = np.array(by[T]['sd'])
    for lv in (50, 80, 95):
        out.append(cover(pred, obs, sd, lv))
    if T == 'Bio_a_ha':
        # thirds of the observed value: exposes a failure at the high end
        q = np.quantile(obs, [1/3, 2/3])
        terc = np.digitize(obs, q)                                 # 0=low,1=mid,2=high
        for k, lab in ((0, 'low_third'), (1, 'mid_third'), (2, 'high_third')):
            m = terc == k
            out.append(cover(pred[m], obs[m], sd[m], 80, stratum=lab))
        # thirds of the predicted value: the split a user can actually make
        qp = np.quantile(pred, [1/3, 2/3]); tp = np.digitize(pred, qp)
        for k, lab in ((0, 'pred_low_third'), (1, 'pred_mid_third'), (2, 'pred_high_third')):
            m = tp == k
            out.append(cover(pred[m], obs[m], sd[m], 80, stratum=lab))

with io.open(DST, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['target', 'level', 'nominal', 'coverage', 'mean_width', 'n', 'stratum'])
    w.writeheader()
    for r in out: w.writerow(r)

print('wrote', os.path.relpath(DST, ROOT))
for r in out:
    print(f"  {r['target']:16s} {r['level']:>3s} nominal={r['nominal']:.2f}  coverage={r['coverage']:.3f}  "
          f"width={r['mean_width']:>8.3f}  n={r['n']:5d}  {r['stratum']}")
