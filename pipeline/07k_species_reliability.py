"""07k_species_reliability.py — per-class reliability table for the dominant-species map.

A single overall accuracy is dominated by the common classes and says nothing about the ones a user
is most likely to be unsure of. This emits one row per species: prevalence, producer's accuracy,
user's accuracy and F1, so the map can be trusted class by class rather than as a whole.

The table ships with the dataset, next to the two dominant-species rasters.

INPUT   data/species_codes.csv                                (code, Korean name, scientific name)
        outputs/nfi8/m4f_domsp_perclass_wholebase.csv         (predset, scheme, species, n_obs, PA, UA, F1)
        The "_wholebase" suffix names the SAMPLE, not the predictor set. Unsuffixed m4f_* files are
        the common subset of plots that have a complete record under both predictor sets; _wholebase
        is every plot with a complete base predictor record, which is the reported standard.
OUTPUT  outputs/nfi8/dominant_species_table.csv
RUN     python pipeline/07k_species_reliability.py
"""
import csv
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')

leg = list(csv.DictReader(io.open(os.path.join(ROOT, 'data', 'species_codes.csv'), encoding='utf-8-sig')))
code2sci = {int(r['code']): r['sci_name'] for r in leg}
code2kr = {int(r['code']): r['species'] for r in leg}

perclass = [r for r in csv.DictReader(
    io.open(os.path.join(OUT, 'm4f_domsp_perclass_wholebase.csv'), encoding='utf-8-sig'))
    if r['predset'] == 'base' and r['scheme'] == 'spatial']
by_kr = {r['species']: r for r in perclass}
total_obs = sum(int(r['n_obs']) for r in perclass)

dst = os.path.join(OUT, 'dominant_species_table.csv')
with io.open(dst, 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['code', 'species_kr', 'sci_name', 'n_obs', 'prevalence',
                'producers_acc', 'users_acc', 'f1'])
    for code in sorted(code2sci):
        kr = code2kr[code]
        r = by_kr.get(kr, {})
        n = int(r['n_obs']) if r.get('n_obs') else 0
        w.writerow([code, kr, code2sci[code], n,
                    round(n / total_obs, 4) if total_obs else '',
                    r.get('producers_acc', ''), r.get('users_acc', ''), r.get('f1', '')])

print('%d classes, %d evaluation plots -> %s'
      % (len(code2sci), total_obs, os.path.relpath(dst, ROOT)))
