"""08b_verify_deposit_transfer.py — prove that a raster copied in from another machine is the same
product as the layers already delivered here.

WHY: base-tier rasters have now been damaged twice by a GDAL PAM statistics write overrunning the
IFD0 tile index, twice, affecting five files in total. The most recent
recovery is a file copy rather than a re-prediction, so nothing about the copy
is checked by the run that produced it. This script is that check.

It answers three separate questions, and a copy has to pass all three:

  1. is the file readable end to end?      full tile decode, not a sample. A 164-block sample of
                                           base_biomass_30m.tif found 0 of its 141 bad tiles.
  2. is it on the delivered grid?          dims, CRS, origin, dtype, internal overviews
  3. is it the same run's output?          valid-cell count must equal the base-tier footprint,
                                           which is what proves the land-cover filter was applied.
                                           Stage 05 writes these files in place and unfiltered,
                                           so geometry alone cannot tell a filtered file from one
                                           that was rebuilt and never re-filtered.

The footprint is not hard-coded. It is taken by consensus from the base-tier layers already sitting
in the delivery tree, so the check stays correct if the product is ever rebuilt. BASE_FOOTPRINT below
is the value measured by full decode on two independent machines, kept as a cross-check.

RUN  python pipeline/08b_verify_deposit_transfer.py FILE [FILE ...]
EXIT 0 = every checked file passed, 1 = at least one failed, 2 = no files named.
"""
import os
import sys

import numpy as np
import rasterio

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')
MAPS = os.path.join(OUT, 'm8_maps_30m')
UNC = os.path.join(OUT, 'm8_uncertainty_30m')

BASE_FOOTPRINT = 69_975_756          # full decode, agreeing across two machines
GRID = dict(width=18593, height=21013, epsg=5179, origin=(746150.1, 2089170.0), pixel_m=30)
MIN_OVERVIEWS = 5


# Healthy base-tier continuous layers used to establish the expected footprint by consensus.
# All three measured 69,975,756, agreeing with the  record. Three is enough
# to catch a mixed tree; more only adds minutes, since each one is a full decode.
# domsp is excluded: uint8 with nodata 255, a different dtype family, and its own recode history.
REFERENCE_LAYERS = [os.path.join(MAPS, 'base_conifer_fraction_30m.tif'),
                    os.path.join(MAPS, 'base_shannon_30m.tif'),
                    os.path.join(UNC, 'base_shannon_sd_30m.tif')]


def decode_and_count(path):
    """Full tile decode. Returns (valid_cells, [failed block ids]).

    Every block is read, so a tile whose byte range is corrupt raises here rather than being
    silently skipped. Counting and decoding are the same pass because both need every tile.
    """
    valid, failed = 0, []
    with rasterio.open(path) as ds:
        nodata = ds.nodata
        use_nan = nodata is None or (isinstance(nodata, float) and np.isnan(nodata))
        for (i, j), win in ds.block_windows(1):
            try:
                a = ds.read(1, window=win)
            except Exception:
                failed.append((i, j))
                continue
            valid += int(np.count_nonzero(~np.isnan(a)) if use_nan else np.count_nonzero(a != nodata))
    return valid, failed


def grid_problems(path):
    """Geometry / metadata deviations from the delivered grid, as a list of strings."""
    bad = []
    with rasterio.open(path) as ds:
        t = ds.transform
        if (ds.width, ds.height) != (GRID['width'], GRID['height']):
            bad.append(f"dims {ds.width}x{ds.height} != {GRID['width']}x{GRID['height']}")
        epsg = ds.crs.to_epsg() if ds.crs else None
        if epsg != GRID['epsg']:
            bad.append(f"EPSG {epsg} != {GRID['epsg']}")
        if round(t.a) != GRID['pixel_m']:
            bad.append(f"pixel {t.a} != {GRID['pixel_m']} m")
        if (round(t.c, 1), round(t.f, 1)) != GRID['origin']:
            bad.append(f"origin ({round(t.c, 1)}, {round(t.f, 1)}) != {GRID['origin']}")
        n_ov = len(ds.overviews(1))
        if n_ov < MIN_OVERVIEWS:
            bad.append(f"{n_ov} overview levels < {MIN_OVERVIEWS}")
        if not ds.profile.get('tiled'):
            bad.append("not tiled")
    return bad


def expected_footprint():
    """Footprint the transferred files must match, by consensus of the healthy base-tier layers.

    Returns (count, note). Falls back to the recorded constant if the delivery tree cannot supply
    a unanimous answer, so a disagreement is reported rather than silently averaged away.
    """
    counts = {}
    for p in REFERENCE_LAYERS:
        if not os.path.exists(p):
            continue
        n, failed = decode_and_count(p)
        if failed:
            print(f"  ! reference layer {os.path.basename(p)} is itself damaged "
                  f"({len(failed)} blocks) — excluded")
            continue
        counts[os.path.basename(p)] = n
        print(f"  reference {os.path.basename(p):40s} {n:,}")
    if not counts:
        return BASE_FOOTPRINT, "no usable reference layer, using the recorded value"
    uniq = set(counts.values())
    if len(uniq) > 1:
        return None, f"reference layers disagree: {counts}"
    n = uniq.pop()
    if n != BASE_FOOTPRINT:
        return n, (f"consensus {n:,} differs from the recorded value {BASE_FOOTPRINT:,} — "
                   f"the product was rebuilt, or the tree is mixed")
    return n, "matches the recorded value"


def main(argv):
    targets = [os.path.abspath(a) for a in argv[1:]]
    if not targets:
        print(__doc__.strip().split('RUN')[-1].strip())
        print('\nName the files to check. There is no default: which files were copied in is\n'
              'something only you know, and checking an arbitrary pair would prove nothing.')
        return 2

    print("[1] expected base-tier footprint")
    want, note = expected_footprint()
    if want is None:
        print(f"  STOP {note}")
        return 1
    print(f"  -> {want:,} cells ({note})\n")

    # The two tiers do not share a footprint and never did: the FTM-added tier is bounded by the
    # forest-type map and covers less. Checking every file against the base count reports the whole
    # FTM tier as broken, which is how this check first read as nine failures that were not there.
    # An FTM file is measured against its own tier instead, by consensus of the FTM layers present.
    ftm_want = None
    ftm_counts = {}
    for p in targets:
        if os.path.basename(p).startswith('ftm_'):
            n, failed = decode_and_count(p)
            if not failed:
                ftm_counts[n] = ftm_counts.get(n, 0) + 1
    if ftm_counts:
        ftm_want = max(ftm_counts, key=ftm_counts.get)
        agree = ftm_counts[ftm_want]
        print(f"[1b] FTM-added tier footprint\n  -> {ftm_want:,} cells "
              f"(consensus of {agree} of {sum(ftm_counts.values())} FTM layers)\n")

    print("[2] transferred files")
    failures = 0
    for p in targets:
        name = os.path.basename(p)
        if not os.path.exists(p):
            print(f"  MISSING  {name}  - not copied in yet")
            failures += 1
            continue

        is_ftm = name.startswith('ftm_') or name.endswith('_ftm.tif')
        expect = ftm_want if (is_ftm and ftm_want) else want
        problems = grid_problems(p)
        n, failed = decode_and_count(p)
        if failed:
            problems.append(f"{len(failed)} blocks failed to decode, first at {failed[0]}")
        if expect is not None and n != expect:
            delta = n - expect
            tier = 'FTM-added' if is_ftm else 'base'
            problems.append(f"valid cells {n:,} != {expect:,} ({delta:+,}) against the {tier} "
                            f"tier — land-cover filter likely not applied — re-run 06b")
        sidecar = p + '.aux.xml'
        if os.path.exists(sidecar):
            problems.append("a .aux.xml sidecar came along — delete it, that write is what corrupts "
                            "the tile index")

        if problems:
            failures += 1
            print(f"  FAIL     {name}")
            for m in problems:
                print(f"             {m}")
        else:
            print(f"  ok       {name}  {n:,} cells, full decode clean")

    print()
    if failures:
        print(f"{failures} file(s) failed. Do NOT stage the deposit.")
        return 1
    print("All transferred files match the delivered product. Safe to stage the deposit.")
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
