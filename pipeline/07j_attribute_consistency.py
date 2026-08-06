"""07j_attribute_consistency.py — measure the multivariate consistency the plot-index product offers.

The plot-index family is described as preserving the realistic co-occurrence of species, sizes and
ages, while every per-attribute accuracy reported for it is worse than the direct family's. That
trade-off is usually asserted rather than measured, here and in the source literature. Riley et al.
(2022) note that coherence between separately imputed attributes is not guaranteed and could produce
pixel combinations that are unlikely or impossible, and call the comparison between the two
approaches a next step rather than a result. Having both families over the same country on the
same grid makes the measurement possible. Three of them, each covering one half of the trade-off:

  (1) VARIANCE RETENTION   Var(map)/Var(plots) per continuous attribute.
      A direct regression maps the conditional mean, so it must under-disperse; a k = 1 transfer copies
      a real value, so it should sit near 1. This measures the "retaining more variability on the
      landscape" claim that motivates k = 1 (Riley et al. 2016).

  (2) PLAUSIBILITY DISTANCE (bin-free)   For each map cell, the distance in standardized attribute
      space to the nearest REAL plot of the same dominant species. The plot-index map is 0 by
      construction, which is what the measure is for. The
      informative number is the DIRECT map's distance, benchmarked against how far apart real plots of
      the same species naturally sit (each plot's nearest neighbour among its own species). A direct
      cell further from any real plot than real plots are from each other is an attribute combination
      the inventory never observed.
      Bin-free on purpose: an "impossible combination" rate over binned attributes is dominated by bin
      width and by the sparsity of 17,099 plots spread over a 4-way cross-tabulation, so it would
      measure our binning rather than the maps.

  (3) MARGINAL + JOINT DISTRIBUTIONAL FIDELITY   Bhattacharyya coefficient between the map and the
      plot population, for the dominant-species marginal (Riley's V4 class-proportion comparison) and
      for the joint dominant-species x composition x biomass x diversity table. The NFI is a systematic
      probability sample, so its distribution is the reference — the same inferential move Riley makes.

Comparing an area-weighted map against a plot sample is only fair because the inventory grid is
systematic. That assumption is what makes (3) meaningful, and it is the one to check first if the
result looks wrong.

INPUT   outputs/nfi8/m2_targets.csv                             (the plot population and its four attributes)
        outputs/nfi8/m8_maps_30m/base_*_30m.tif              (direct family, base tier)
        outputs/nfi8/attributes_base/spec30m_*.tif              (plot-index family, base tier)
OUTPUT  outputs/nfi8/m9_consistency.csv                         (tidy: metric, family, attribute, value)
RUN     python pipeline/07j_attribute_consistency.py
"""
import os, csv, io
import numpy as np
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT  = os.path.join(ROOT, 'outputs', 'nfi8')
DIRECT = {'coniferfrac': 'm8_maps_30m/base_conifer_fraction_30m.tif',
          'Bio_a_ha':    'm8_maps_30m/base_biomass_30m.tif',
          'shannon':     'm8_maps_30m/base_shannon_30m.tif',
          'domsp':       'm8_maps_30m/base_dominant_species_30m.tif'}
IMPUTE = {'coniferfrac': 'attributes_base/spec30m_coniferfrac.tif',
          'Bio_a_ha':    'attributes_base/spec30m_Bio_a_ha.tif',
          'shannon':     'attributes_base/spec30m_shannon.tif',
          'domsp':       'attributes_base/spec30m_dom_species.tif'}
NSAMPLE = 200_000
SEED    = 7214

# ---- dominant-species harmonisation -------------------------------------------------------------
# The three sources code species differently and must be put on ONE basis before anything is compared:
#   plot table         Korean name, with the two regional races of Pinus densiflora still separate
#   direct raster      integer 0-38 on the 39-class legend (pines already merged, rare pooled to 기타)
#   imputation raster  integer on the FULL 141-species NFI legend (no pooling at all)
# Canonical basis = the 39-class legend used by the paper: merge the regional pines, pool anything
# outside the 38 named species into 기타.
def read_legend(path, kc='code', sc='species'):
    # utf-8-sig: these legends are written by R with a BOM, which would otherwise become part of the
    # first column NAME ('﻿code') and make the lookup fail
    return {int(r[kc]): r[sc] for r in csv.DictReader(io.open(path, encoding='utf-8-sig'))}
LEG39 = read_legend(os.path.join(ROOT, 'data', 'species_codes.csv'))                       # code -> name (39)
LEGIMP = read_legend(os.path.join(OUT, 'attributes_base', 'spec30m_dom_species_legend.csv'))
PINE = {'강원지방소나무': '소나무', '중부지방소나무': '소나무'}
NAME2CODE = {v: k for k, v in LEG39.items()}
def canon(name):
    name = PINE.get(name, name)
    return NAME2CODE.get(name, NAME2CODE['기타'])                                     # unseen -> 기타

# ---- reference: the NFI plot population -------------------------------------------------------
rows = list(csv.DictReader(io.open(os.path.join(OUT, 'm2_targets.csv'), encoding='utf-8-sig')))
# the plot table names the whole-stand variants; this script works in the short names throughout
COLUMN = {'coniferfrac': 'coniferfrac_full', 'shannon': 'shannon_full',
          'Bio_a_ha': 'Bio_a_ha', 'domsp': 'domsp_full'}
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
def sample_family(paths, n):
    """Draw n random in-raster cells from the family's own grid, return dict of arrays (finite rows only)."""
    ref = os.path.join(OUT, paths['coniferfrac'])
    with rasterio.open(ref) as s:
        H, W, tr = s.height, s.width, s.transform
        # oversample: most of the bounding box is non-forest nodata
        take, got = {}, 0
        rr = rng.integers(0, H, n * 8); cc = rng.integers(0, W, n * 8)
        xs = tr.c + (cc + 0.5) * tr.a; ys = tr.f + (rr + 0.5) * tr.e
    pts = list(zip(xs, ys))
    vals = {}
    for k, rel in paths.items():
        with rasterio.open(os.path.join(OUT, rel)) as s:
            nod = s.nodata
            v = np.array([x[0] for x in s.sample(pts)], dtype=float)
            if nod is not None and np.isfinite(nod): v[v == nod] = np.nan
            vals[k] = v
    keep = np.all([np.isfinite(vals[k]) for k in vals], axis=0)
    keep &= vals['domsp'] >= 0
    idx = np.flatnonzero(keep)[:n]
    return {k: v[idx] for k, v in vals.items()}

print('sampling direct family ...');   D = sample_family(DIRECT, NSAMPLE)
print('sampling imputation family ...'); I = sample_family(IMPUTE, NSAMPLE)
# put the imputation raster's 141-species coding onto the canonical 39-class basis
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
# standardize on the PLOT population so all three sets share one scale; log1p biomass (right-skewed)
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
                d[d == 0] = np.inf                      # a plot's nearest neighbour is not itself
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
    'EXACTLY zero by construction — every cell IS a real plot; reported for completeness, not as a result')
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

CF = [0, .01, .25, .5, .75, .99, 1.001]                       # composition: the two pure ends matter
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
dst = os.path.join(OUT, 'm9_consistency.csv')
with io.open(dst, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['metric', 'family', 'attribute', 'value', 'note']); w.writeheader()
    for r in res: w.writerow(r)
print('\nwrote', os.path.relpath(dst, ROOT))
for r in res:
    print(f"  {r['metric']:34s} {r['family']:22s} {r['attribute']:38s} {r['value']:.4f}")
