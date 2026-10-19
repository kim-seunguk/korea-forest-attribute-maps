"""07q_density_weighting.py — accuracy with every plot weighted by the area it represents.

The inventory is not a uniform 4 km sample: some regions have supplementary clusters on 1-2 km
grids, so when every plot counts equally, as it does in every other validation stage, a plot on a
1 km grid counts sixteen times as much per unit of forest area as a plot on the 4 km grid does.
This reports the same out-of-fold predictions with each plot weighted by the square of its local
grid spacing, next to the unweighted figure, as a sensitivity. It fits nothing.

Local spacing = the median distance from a cluster's center to its three nearest cluster centers,
snapped to the grid it sits on (< 1.4 km -> 1 km, < 2.8 km -> 2 km, otherwise 4 km); the median of
three keeps a 4 km cluster beside a densified area at 4 km. Weight = spacing^2, km^2 per cluster.
The share of plots on the 1-2 km grids is printed.

INPUT   outputs/nfi8/cv_folds.csv                         (Plot_ID, UTMK_X, UTMK_Y)
        outputs/nfi8/predictor_ladder_points_ftmonly.csv  (direct route, spatial-CV out-of-fold; 07a with FTM_ONLY)
        outputs/nfi8/family_comparison_points.csv         (plot-index route, spatial-CV out-of-fold; 07b, predset full)
OUTPUT  outputs/nfi8/density_weighted.csv                 (route, target, metric, n, unweighted, weighted)
        outputs/nfi8/plot_spacing.csv                     (Plot_ID, cluster, spacing_km, weight): it carries the
                                                           plot identifier, so it stays in outputs/ and is not published
RUN     python pipeline/07q_density_weighting.py
"""
import csv
import io
import math
import os
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')


def rd(name):
    return list(csv.DictReader(io.open(os.path.join(OUT, name), encoding='utf-8-sig')))


# ---- local grid spacing per cluster ------------------------------------------------------------------
folds = rd('cv_folds.csv')
cxy = defaultdict(list)
for r in folds:
    cxy[r['Plot_ID'][:12]].append((float(r['UTMK_X']), float(r['UTMK_Y'])))   # cluster = digits 1-12 of the plot id
cl = sorted(cxy)
C = np.array([np.mean(cxy[c], axis=0) for c in cl])
cell = 5000.0                                              # spatial hash for the nearest-neighbor search
buckets = defaultdict(list)
for i, (x, y) in enumerate(C):
    buckets[(int(x // cell), int(y // cell))].append(i)
nn = np.full(len(C), np.inf)
for i, (x, y) in enumerate(C):
    bx, by = int(x // cell), int(y // cell)
    d = sorted(math.hypot(C[j][0] - x, C[j][1] - y)
               for dx in (-1, 0, 1) for dy in (-1, 0, 1) for j in buckets.get((bx + dx, by + dy), ()) if j != i)
    if d:
        nn[i] = float(np.median(d[:3]))
spacing = np.where(nn < 1400, 1.0, np.where(nn < 2800, 2.0, 4.0))   # isolated clusters (no neighbor within 5 km) -> 4
sp_of = dict(zip(cl, spacing))
w_of = {c: s * s for c, s in sp_of.items()}
plots = [r['Plot_ID'] for r in folds]
with io.open(os.path.join(OUT, 'plot_spacing.csv'), 'w', encoding='utf-8', newline='') as f:
    wr = csv.writer(f); wr.writerow(['Plot_ID', 'cluster', 'spacing_km', 'weight'])
    for p in plots:
        wr.writerow([p, p[:12], sp_of[p[:12]], w_of[p[:12]]])
cnt = Counter(sp_of[p[:12]] for p in plots)
print('plots by local grid spacing (km):', dict(sorted(cnt.items())),
      '| share on 1-2 km grids: %.1f %%' % (100.0 * (cnt[1.0] + cnt[2.0]) / len(plots)))


# ---- weighted metrics ---------------------------------------------------------------------------------
def r2(o, p, w):
    mo = np.sum(w * o) / np.sum(w)
    return 1.0 - np.sum(w * (o - p) ** 2) / np.sum(w * (o - mo) ** 2)


def acc(o, p, w):
    return float(np.sum(w * (o == p)) / np.sum(w))


rows = []
for route, fname, predset in (('direct', 'predictor_ladder_points_ftmonly.csv', 'full'),
                              ('imputation', 'family_comparison_points.csv', 'full')):
    path = os.path.join(OUT, fname)
    if not os.path.exists(path):
        print('missing', fname, '- skipped'); continue
    pts = [r for r in rd(fname) if r['predset'] == predset]
    for tg in sorted({r['target'] for r in pts}):
        g = [r for r in pts if r['target'] == tg and r['obs'] not in ('', 'NA') and r['pred'] not in ('', 'NA')]
        if not g:
            continue
        w = np.array([w_of[r['Plot_ID'][:12]] for r in g], dtype=float)
        if tg.startswith('domsp') or tg.startswith('fld_'):
            o = np.array([r['obs'] for r in g]); p = np.array([r['pred'] for r in g])
            rows.append((route, tg, 'accuracy', len(g), acc(o, p, np.ones_like(w)), acc(o, p, w)))
        else:
            o = np.array([float(r['obs']) for r in g]); p = np.array([float(r['pred']) for r in g])
            rows.append((route, tg, 'r2', len(g), r2(o, p, np.ones_like(w)), r2(o, p, w)))
with io.open(os.path.join(OUT, 'density_weighted.csv'), 'w', encoding='utf-8', newline='') as f:
    wr = csv.writer(f); wr.writerow(['route', 'target', 'metric', 'n', 'unweighted', 'weighted'])
    for r in rows:
        wr.writerow([r[0], r[1], r[2], r[3], round(r[4], 4), round(r[5], 4)])
for r in rows:
    print('%-10s %-18s %-8s n=%5d  unweighted %.4f  weighted %.4f' % r)
