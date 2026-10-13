"""07k_species_reliability.py — per-class reliability table for the dominant-species map.

A single overall accuracy is dominated by the common classes and says nothing about the ones a user
is most likely to be unsure of. This emits one row per trained class: prevalence, producer's
accuracy, user's accuracy and F1.

The rows are the classes the species map was trained on, m8_domsp_legend_used.csv from 05a: a legend
code that the rare-species rule pooled into the "other" class has no cells on the map and no row
here. The table ships with the dataset, next to the two dominant-species rasters; 08d copies it in
under its published name.

INPUT   outputs/nfi8/m8_maps_30m/m8_domsp_legend_used.csv     (code, code_v1, Korean and scientific name)
        outputs/nfi8/m4f_domsp_perclass_wholeftm.csv          (predset, scheme, species, n_obs, PA, UA, F1)
        The "_wholeftm" suffix names the SAMPLE, not the predictor set: every plot with a complete
        record under the delivered set, which 07d writes with FTM_ONLY. Unsuffixed m4f_* files are
        the common subset of plots complete under both predictor sets.
OUTPUT  outputs/nfi8/dominant_species_table_ftm.csv
RUN     python pipeline/07k_species_reliability.py
"""
import csv
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')

LEGEND = os.path.join(OUT, 'm8_maps_30m', 'm8_domsp_legend_used.csv')
if not os.path.exists(LEGEND):
    raise SystemExit('missing %s: run pipeline/05a_direct_prediction.R first, it writes the trained legend'
                     % os.path.relpath(LEGEND, ROOT))
leg = list(csv.DictReader(io.open(LEGEND, encoding='utf-8-sig')))
code2sci = {int(r['code']): r['sci_name'] for r in leg}
code2kr = {int(r['code']): r['species'] for r in leg}
code2v1 = {int(r['code']): r.get('code_v1', '') for r in leg}   # the class's code in the earlier record

perclass = [r for r in csv.DictReader(
    io.open(os.path.join(OUT, 'm4f_domsp_perclass_wholeftm.csv'), encoding='utf-8-sig'))
    if r['predset'] == 'full' and r['scheme'] == 'spatial']
by_kr = {r['species']: r for r in perclass}
total_obs = sum(int(r['n_obs']) for r in perclass)

dst = os.path.join(OUT, 'dominant_species_table_ftm.csv')
with io.open(dst, 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['code', 'code_v1', 'species_kr', 'sci_name', 'n_obs', 'prevalence',
                'producers_acc', 'users_acc', 'f1'])
    for code in sorted(code2sci):
        kr = code2kr[code]
        r = by_kr.get(kr, {})
        n = int(r['n_obs']) if r.get('n_obs') else 0
        w.writerow([code, code2v1[code], kr, code2sci[code], n,
                    round(n / total_obs, 4) if total_obs else '',
                    r.get('producers_acc', ''), r.get('users_acc', ''), r.get('f1', '')])

print('%d trained classes, %d evaluation plots -> %s'
      % (len(code2sci), total_obs, os.path.relpath(dst, ROOT)))
