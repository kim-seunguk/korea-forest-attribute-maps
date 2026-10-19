"""07j_attribute_consistency.py — measure the multivariate consistency the plot-index product offers.

Each cell of the plot-index family holds a combination of attributes observed together on one
plot, while every per-attribute accuracy reported for it is worse than the direct family's. The
script computes three measures of this trade-off:

  (1) VARIANCE RETENTION   Var(map)/Var(plots) per continuous attribute.
      A direct regression maps the conditional mean, so it must under-disperse; a k = 1 transfer copies
      a real value, so it should sit near 1.

  (2) PLAUSIBILITY DISTANCE (bin-free)   For each map cell, the distance in standardized attribute
      space to the nearest REAL plot of the same dominant species. The plot-index map is 0 by
      construction, so the informative number is the DIRECT map's distance, benchmarked against how
      far apart real plots of the same species naturally sit (each plot's nearest neighbor among its
      own species). A direct cell further from any real plot than real plots are from each other is an
      attribute combination the inventory never observed.

  (3) MARGINAL + JOINT DISTRIBUTIONAL FIDELITY   Bhattacharyya coefficient between the map and the
      plot population, for the dominant-species marginal (the class-proportion comparison of
      Riley et al.) and for the joint dominant-species x composition x biomass x diversity table.
      The NFI is a systematic probability sample, so its distribution is the reference.

07q measures the effect of the inventory's denser supplementary grids.

Two restrictions apply. The reference population is the plots the maps cover: the
maps stop at the forest-type map's stocked extent, and the plots outside it are a different
population from those inside, so they are left out before anything is compared. And both families
are sampled at the same random cells, since they share one grid; a cell is kept only where both
hold a value.

INPUT   outputs/nfi8/targets.csv                                 (the plot population and its four attributes)
        data/interim/nfi8/rd_points.csv                          (plot positions, for the mapped-extent restriction)
        outputs/nfi8/maps_30m/ftm_*_30m.tif                      (direct family)
        outputs/nfi8/attributes_ftm/spec30m_*.tif                (plot-index family)
        outputs/nfi8/maps_30m/species_legend_used.csv            (the classes the species map was trained on)
OUTPUT  outputs/nfi8/attribute_consistency_ftm.csv               (tidy: metric, family, attribute, value)
RUN     python pipeline/07j_attribute_consistency.py
"""
import os, csv, io
import numpy as np
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT  = os.path.join(ROOT, 'outputs', 'nfi8')
PREDSET = 'ftm'                                # the delivered file prefix (R/paths.R, FILE_PREDSET)
DIRECT = {'coniferfrac': f'maps_30m/{PREDSET}_conifer_fraction_30m.tif',
          'Bio_a_ha':    f'maps_30m/{PREDSET}_biomass_30m.tif',
          'shannon':     f'maps_30m/{PREDSET}_shannon_30m.tif',
          'domsp':       f'maps_30m/{PREDSET}_dominant_species_30m.tif'}
IMPUTE = {'coniferfrac': f'attributes_{PREDSET}/spec30m_coniferfrac.tif',
          'Bio_a_ha':    f'attributes_{PREDSET}/spec30m_Bio_a_ha.tif',
          'shannon':     f'attributes_{PREDSET}/spec30m_shannon.tif',
          'domsp':       f'attributes_{PREDSET}/spec30m_dom_species.tif'}
NSAMPLE = 200_000
SEED    = 7214

# ---- dominant-species harmonisation -------------------------------------------------------------
# The three sources code species differently and must be put on ONE basis before anything is compared:
#   plot table         Korean name, with the two regional races of Pinus densiflora still separate
#   direct raster      integer code on the trained legend (pines merged, rare species pooled to 기타)
#   imputation raster  integer on the inventory's own species list, as 05c wrote it (no pooling at all)
# Canonical basis = the legend the species map was trained on, species_legend_used.csv from 05a:
# merge the regional pines, pool anything outside its named species into 기타. data/species_codes.csv
# is read when 05a has not run, but it lists every legend code, including any the rare-species rule
# pooled, so the trained legend is the one to compare with.
def read_legend(path, kc='code', sc='species'):
    # utf-8-sig: these legends are written by R with a BOM, which would otherwise become part of the
    # first column NAME ('﻿code') and make the lookup fail
    return {int(r[kc]): r[sc] for r in csv.DictReader(io.open(path, encoding='utf-8-sig'))}
LEG_USED = os.path.join(OUT, 'maps_30m', 'species_legend_used.csv')
LEG = read_legend(LEG_USED if os.path.exists(LEG_USED) else os.path.join(ROOT, 'data', 'species_codes.csv'))
LEGIMP = read_legend(os.path.join(OUT, f'attributes_{PREDSET}', 'spec30m_dom_species_legend.csv'))
PINE = {'강원지방소나무': '소나무', '중부지방소나무': '소나무'}
NAME2CODE = {v: k for k, v in LEG.items()}
def canon(name):
    name = PINE.get(name, name)
    return NAME2CODE.get(name, NAME2CODE['기타'])                                     # unseen -> 기타

# ---- reference: the NFI plot population, restricted to the plots the maps cover ------------------
rows = list(csv.DictReader(io.open(os.path.join(OUT, 'targets.csv'), encoding='utf-8-sig')))
# the plot table names the whole-stand variants; this script works in the short names throughout
COLUMN = {'coniferfrac': 'coniferfrac_full', 'shannon': 'shannon_full',
          'Bio_a_ha': 'Bio_a_ha', 'domsp': 'domsp_full'}
xy = {r['Plot_ID']: (float(r['UTMK_X']), float(r['UTMK_Y'])) for r in csv.DictReader(
    io.open(os.path.join(ROOT, 'data', 'interim', 'nfi8', 'rd_points.csv'), encoding='utf-8-sig'))}
rows = [r for r in rows if r['Plot_ID'] in xy]
with rasterio.open(os.path.join(OUT, DIRECT['coniferfrac'])) as s:
    at_plot = np.array([x[0] for x in s.sample([xy[r['Plot_ID']] for r in rows])], dtype=float)
    if s.nodata is not None and np.isfinite(s.nodata): at_plot[at_plot == s.nodata] = np.nan
rows = [r for r, v in zip(rows, at_plot) if np.isfinite(v)]
print(f'reference plots on a mapped cell: {len(rows)}')
def col(name, cast=float):
    v = []
    for r in rows:
        x = r.get(COLUMN[name], '')
        v.append(np.nan if x in ('', 'NA', None) else cast(x))
    return np.asarray(v, dtype=float)
P = {k: col(k) for k in ('coniferfrac', 'shannon', 'Bio_a_ha')}
P['domsp'] = np.array([canon(r[COLUMN['domsp']]) if r.get(COLUMN['domsp']) else np.nan
                       for r in rows], dtype=float)
ok = np.all([np.isfinite(P[k]) for k in P], axis=0)
P = {k: v[ok] for k, v in P.items()}
print(f'reference plots: {ok.sum()} of {len(rows)}')

# ---- sample the two map families at the SAME random cells --------------------------------------
rng = np.random.default_rng(SEED)
def read_at(paths, pts):
    vals = {}
    for k, rel in paths.items():
        with rasterio.open(os.path.join(OUT, rel)) as s:
            nod = s.nodata
            v = np.array([x[0] for x in s.sample(pts)], dtype=float)
            if nod is not None and np.isfinite(nod): v[v == nod] = np.nan
            vals[k] = v
    return vals
with rasterio.open(os.path.join(OUT, DIRECT['coniferfrac'])) as s:
    H, W, tr = s.height, s.width, s.transform
rr = rng.integers(0, H, NSAMPLE * 8); cc = rng.integers(0, W, NSAMPLE * 8)   # oversample: most of the box is nodata
pts = list(zip(tr.c + (cc + 0.5) * tr.a, tr.f + (rr + 0.5) * tr.e))
print('sampling both families at the same cells ...')
Dv, Iv = read_at(DIRECT, pts), read_at(IMPUTE, pts)
keep = np.all([np.isfinite(Dv[k]) for k in Dv] + [np.isfinite(Iv[k]) for k in Iv], axis=0)
keep &= (Dv['domsp'] >= 0) & (Iv['domsp'] >= 0)
idx = np.flatnonzero(keep)[:NSAMPLE]
D = {k: v[idx] for k, v in Dv.items()}; I = {k: v[idx] for k, v in Iv.items()}
# put the imputation raster's own species coding onto the trained legend
I['domsp'] = np.array([canon(LEGIMP.get(int(c), '기타')) for c in I['domsp']], dtype=float)
print(f'  direct n={len(D["coniferfrac"]):,}   imputation n={len(I["coniferfrac"]):,}')

res = []
def add(metric, family, attribute, value, note=''):
    res.append(dict(metric=metric, family=family, attribute=attribute, value=value, note=note))

# ---- (1) variance retention --------------------------------------------------------------------
for a in ('coniferfrac', 'Bio_a_ha', 'shannon'):
    vp = np.nanvar(P[a])
    for fam, M in (('direct', D), ('imputation', I)):
        add('variance_ratio', fam, a, float(np.nanvar(M[a]) / vp),
            'Var(map)/Var(plots); 1 = field variance retained, <1 = smoothed toward the mean')

# ---- (2) plausibility distance ------------------------------------------------------------------
# Standardize on the PLOT population so all three sets share one scale.
# Biomass is log1p-transformed because it is right-skewed.
def feat(d):
    return np.column_stack([d['coniferfrac'], np.log1p(np.clip(d['Bio_a_ha'], 0, None)), d['shannon']])
mu, sd = feat(P).mean(0), feat(P).std(0)
Zp, Zd, Zi = [(feat(x) - mu) / sd for x in (P, D, I)]

def nn_dist(Q, qsp, R, rsp, exclude_self=False):
    """Distance from each row of Q to the nearest row of R carrying the SAME dominant species."""
    out = np.full(len(Q), np.nan)
    for sp in np.unique(qsp):
        qi = np.flatnonzero(qsp == sp); ri = np.flatnonzero(rsp == sp)
        if len(ri) == 0: continue                       # species absent from the plot set -> undefined
        A, B = Q[qi], R[ri]
        step = 2000
        for s0 in range(0, len(A), step):
            blk = A[s0:s0 + step]
            d = np.sqrt(((blk[:, None, :] - B[None, :, :]) ** 2).sum(-1))
            if exclude_self:
                d[d == 0] = np.inf                      # a plot's nearest neighbor is not itself
            out[qi[s0:s0 + step]] = d.min(1)
    return out

print('plausibility distance: plot-to-plot baseline ...')
d_base = nn_dist(Zp, P['domsp'], Zp, P['domsp'], exclude_self=True)
print('plausibility distance: direct map ...')
sub = rng.choice(len(Zd), size=min(20000, len(Zd)), replace=False)      # 20k is ample for quantiles
d_dir = nn_dist(Zd[sub], D['domsp'][sub], Zp, P['domsp'])

for nm, d in (('plot_to_plot_baseline', d_base), ('direct', d_dir)):
    dd = d[np.isfinite(d)]
    for q, lab in ((50, 'median'), (90, 'p90')):
        add('plausibility_distance_' + lab, nm, 'joint(coniferfrac,logBio,shannon)',
            float(np.percentile(dd, q)), 'standardized units; nearest real plot of the same species')
add('plausibility_distance_median', 'imputation', 'joint(coniferfrac,logBio,shannon)', 0.0,
    'zero by construction: every cell is a real plot')
frac = float(np.mean(d_dir[np.isfinite(d_dir)] > np.nanpercentile(d_base, 90)))
add('beyond_p90_of_plot_spacing', 'direct', 'joint(coniferfrac,logBio,shannon)', frac,
    'share of direct-map cells further from any real same-species plot than 90% of real plots are from each other')

# ---- (3) distributional fidelity (Bhattacharyya) -------------------------------------------------
def bhat(p, q):
    p = p / p.sum(); q = q / q.sum()
    return float(np.sqrt(p * q).sum())          # 1 = identical, 0 = disjoint

sp_levels = np.unique(np.concatenate([P['domsp'], D['domsp'], I['domsp']]))
def sp_hist(v): return np.array([(v == s).sum() for s in sp_levels], dtype=float) + 1e-12
for fam, M in (('direct', D), ('imputation', I)):
    add('bhattacharyya_domsp_marginal', fam, 'dominant species',
        bhat(sp_hist(P['domsp']), sp_hist(M['domsp'])), 'vs the NFI plot population')

# composition bins, with separate bins for the near-pure ends (<0.01 and >=0.99)
CF = [0, .01, .25, .5, .75, .99, 1.001]
BQ = np.nanpercentile(P['Bio_a_ha'], [0, 33, 67, 100]); BQ[-1] += 1
SQ = np.nanpercentile(P['shannon'],  [0, 33, 67, 100]); SQ[-1] += 1
def joint_hist(d):
    i = np.digitize(d['coniferfrac'], CF[1:-1])
    j = np.digitize(d['Bio_a_ha'],    BQ[1:-1])
    k = np.digitize(d['shannon'],     SQ[1:-1])
    s = np.searchsorted(sp_levels, d['domsp'])
    code = ((s * 6 + i) * 3 + j) * 3 + k
    return np.bincount(code, minlength=len(sp_levels) * 6 * 3 * 3).astype(float) + 1e-12
hp = joint_hist(P)
for fam, M in (('direct', D), ('imputation', I)):
    add('bhattacharyya_joint', fam, 'domsp x composition x biomass x diversity',
        bhat(hp, joint_hist(M)), 'joint table vs the NFI plot population')

# ---- write -------------------------------------------------------------------------------------
dst = os.path.join(OUT, f'attribute_consistency_{PREDSET}.csv')
with io.open(dst, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['metric', 'family', 'attribute', 'value', 'note']); w.writeheader()
    for r in res: w.writerow(r)
print('\nwrote', os.path.relpath(dst, ROOT))
for r in res:
    print(f"  {r['metric']:34s} {r['family']:22s} {r['attribute']:38s} {r['value']:.4f}")
