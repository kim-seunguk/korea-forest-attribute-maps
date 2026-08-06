"""06c_recode_species_codes.py  — recode the delivered dominant-species maps to the sci-name-ordered codes.

Class codes come out of training in Korean-name order, which scatters congeners across the legend.
The published maps use scientific-name order instead, so that a genus occupies a contiguous block
of codes and a reader scanning the legend sees related species together. data/species_code_remap.csv
holds the permutation.

It is applied to the two dominant-species rasters only. The uncertainty layers beside them are
continuous and carry no class codes.

SAFETY: writes a temp file, then RE-READS it independently and cross-tabulates the written histogram
against the expected (permuted) source histogram; only on an exact match does it atomically replace
the original (os.replace). Any mismatch aborts with the temp left for inspection and the original
untouched. NaN (off-forest) is preserved structurally (only finite class pixels are remapped).

Run:  python pipeline/06c_recode_species_codes.py            # both tiers
      python pipeline/06c_recode_species_codes.py base_dominant_species_30m.tif   # one tier only

! The recode is NOT idempotent. The table is a permutation, so applying it twice re-permutes an
  already-recoded map into nonsense, and the verification above cannot catch that: it confirms the
  permutation was applied, not that it was warranted. What does catch it is a tag written into the
  raster on the way out; a file carrying it is skipped. Do not remove that guard.
"""
import os
import sys
import numpy as np
import pandas as pd
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
LUT_CSV = os.path.join(ROOT, 'data', 'species_code_remap.csv')
MAPS = os.path.join(ROOT, 'outputs', 'nfi8', 'm8_maps_30m')
FILES = sys.argv[1:] if len(sys.argv) > 1 else ['base_dominant_species_30m.tif', 'ftm_dominant_species_30m.tif']
NCLASS = 39
RECODE_TAG = 'sci_name_order_v1'   # stamped on the map as it is recoded; a second run sees it and skips
#   (the recode is a permutation, so applying it twice garbles the map). A freshly written raster
#   carries no tag, so a rebuilt map is recoded exactly once.

lut_df = pd.read_csv(LUT_CSV)
LUT = {int(o): int(n) for o, n in zip(lut_df.old_code, lut_df.new_code)}
assert sorted(LUT.keys()) == list(range(NCLASS)), 'LUT keys must be 0..38'
assert sorted(LUT.values()) == list(range(NCLASS)), 'LUT values must be a permutation of 0..38'
lut_arr = np.full(NCLASS, -1, dtype=np.int64)
for o, n in LUT.items():
    lut_arr[o] = n

def hist_of(path):
    """finite-value histogram (per integer code 0..38) and NaN count; asserts clean integer codes."""
    h = np.zeros(NCLASS, dtype=np.int64); nan = 0
    with rasterio.open(path) as s:
        for _, win in s.block_windows(1):
            a = s.read(1, window=win); fin = np.isfinite(a)
            nan += int((~fin).sum())
            v = a[fin]; iv = np.rint(v).astype(np.int64)
            if v.size and (((iv < 0) | (iv >= NCLASS)).any() or (np.abs(v - iv) > 1e-6).any()):
                raise ValueError(f'{os.path.basename(path)}: non-integer or out-of-range class value')
            if v.size:
                h += np.bincount(iv, minlength=NCLASS)
    return h, nan

n_recoded = 0
for fn in FILES:
    src_path = os.path.join(MAPS, fn); tmp_path = src_path + '.recoded.tmp'
    # GUARD (recurrence-proof): a map already stamped RECODE_TAG must NOT be recoded again (permutation garbles it)
    with rasterio.open(src_path) as s0:
        if s0.tags().get('DOMSP_RECODE') == RECODE_TAG:
            print(f'{fn}: ALREADY recoded ({RECODE_TAG}) -> SKIP (recode is not idempotent)'); continue
    # pass 0: source histogram
    src_h, src_nan = hist_of(src_path)
    # pass 1: write remapped temp (NaN preserved; only finite class pixels changed); stamp the sentinel tag
    with rasterio.open(src_path) as src:
        prof = src.profile
        with rasterio.open(tmp_path, 'w', **prof) as dst:
            for _, win in src.block_windows(1):
                a = src.read(1, window=win); fin = np.isfinite(a)
                out = a.copy()
                out[fin] = lut_arr[np.rint(a[fin]).astype(np.int64)].astype(a.dtype)
                dst.write(out, 1, window=win)
            dst.update_tags(DOMSP_RECODE=RECODE_TAG)   # mark recoded -> a re-run detects this and skips
    # pass 2: independent re-read of the WRITTEN file + cross-tab vs expected permutation
    new_h, new_nan = hist_of(tmp_path)
    expected = np.zeros(NCLASS, dtype=np.int64)
    for o in range(NCLASS):
        expected[LUT[o]] = src_h[o]
    if new_nan != src_nan or not np.array_equal(new_h, expected):
        raise SystemExit(f'ABORT {fn}: verification failed (nan {src_nan}->{new_nan}); temp kept, original untouched')
    os.replace(tmp_path, src_path)  # atomic
    n_recoded += 1
    print(f'{fn}: OK  finite={src_h.sum():,}  nan={src_nan:,}  classes {int((src_h>0).sum())}  -> swapped')

print('DONE - recoded %d of %d dominant-species map(s) to sci-name order' % (n_recoded, len(FILES)))
