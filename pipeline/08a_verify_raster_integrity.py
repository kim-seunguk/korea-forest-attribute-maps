"""08a_verify_raster_integrity.py — pre-deposit integrity check for every delivered GeoTIFF.

WHY: a tile index can be damaged while the overviews stay intact. A GIS still opens such a file,
and a header-only check (08c_build_data_records.py reads the profile, never the pixels) or a
sampled read misses the damage.

This script parses the TIFF tile tables directly, for every IFD, and flags any tile whose byte range
cannot exist inside the file. By default it decodes no pixels, so it runs in seconds over the whole
product. Three cases a tile-table parse alone would miss or misreport are covered as well:

  - BigTIFF files (8-byte offsets, LONG8 tags) are parsed like classic TIFF. Several writers use
    BIGTIFF=IF_SAFER, so a layer can change format once it passes 4 GB, and the predictor rasters
    under data/interim/ are all BigTIFF.
  - a tile with offset 0 and byte count 0 was never written. None of this pipeline's writers produce
    sparse files, so such a tile marks a file caught mid-write, and it is reported as damage.
  - --decode reads every block of IFD0 through GDAL, read-only and with PAM disabled: the full decode
    the deposit README promises. --inputs also scans the predictor rasters under data/interim/.

RUN  python pipeline/08a_verify_raster_integrity.py [--inputs] [--decode]
EXIT 0 = all clean, 1 = at least one damaged or unparseable file, or nothing found to check (so it
     can gate a deposit script).
"""
import glob
import os
import struct
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')
DECODE = '--decode' in sys.argv
INPUTS = '--inputs' in sys.argv
# PATTERNS covers the delivered rasters and the staged deposit. A deposit gate has to read the
# uploaded copies as well as the masters they were copied from. The deposit glob is empty until 08d
# has run.
PATTERNS = ['maps_30m/*.tif', 'uncertainty_30m/*.tif', 'index/*.tif', 'attributes_ftm/*.tif',
            '_deposit_v2/*.tif']
INPUT_PATTERNS = ['data/interim/covar_5179/*.tif', 'data/interim/covar_5179_s2/*.tif',
                  'data/interim/covar_5179_s2/s2_bands/*.tif', 'data/interim/fmap_30m/*/*.tif',
                  'data/interim/landcover_mcee_30m/*.tif']
TYPE_SIZE = {1: 1, 3: 2, 4: 4, 16: 8}
FMT = {1: 'B', 3: 'H', 4: 'I', 16: 'Q'}
# A 256x256 float64 tile is ~0.5 MB, so a byte count above MAX_TILE_BYTES marks a corrupt tile table.
MAX_TILE_BYTES = 50_000_000


def scan(path):
    """Return (n_tiles, n_bad, n_unwritten, [(ifd, n_bad), ...]), or None if the file is not a TIFF."""
    size = os.path.getsize(path)
    with open(path, 'rb') as f:
        head = f.read(16)
        if head[:2] not in (b'II', b'MM'):
            return None
        e = '<' if head[:2] == b'II' else '>'
        big = struct.unpack(e + 'H', head[2:4])[0] == 43
        if big:
            offset = struct.unpack(e + 'Q', head[8:16])[0]
            cnt_fmt, cnt_len, ent_len, off_fmt, off_len = 'Q', 8, 20, 'Q', 8
        else:
            offset = struct.unpack(e + 'I', head[4:8])[0]
            cnt_fmt, cnt_len, ent_len, off_fmt, off_len = 'H', 2, 12, 'I', 4

        per_ifd, total, unwritten, tiles = [], 0, 0, 0
        ifd = 0
        while offset and 0 < offset < size:
            f.seek(offset)
            n_entries = struct.unpack(e + cnt_fmt, f.read(cnt_len))[0]
            if not 0 < n_entries <= 500:
                break
            tags = {}
            for _ in range(n_entries):
                entry = f.read(ent_len)
                if big:
                    tag, typ, count = struct.unpack(e + 'HHQ', entry[:12]); val = entry[12:20]
                else:
                    tag, typ, count = struct.unpack(e + 'HHI', entry[:8]); val = entry[8:12]
                tags[tag] = (typ, count, val)
            nxt = struct.unpack(e + off_fmt, f.read(off_len))[0]

            def array(tag):
                typ, count, val = tags[tag]
                if typ not in FMT:
                    return None
                nbytes = count * TYPE_SIZE[typ]
                if nbytes <= len(val):
                    data = val[:nbytes]
                else:
                    ptr = struct.unpack(e + off_fmt, val)[0]
                    if ptr + nbytes > size:
                        return None
                    f.seek(ptr)
                    data = f.read(nbytes)
                return list(struct.unpack(e + str(count) + FMT[typ], data))

            if 324 in tags and 325 in tags:
                offs, counts = array(324), array(325)
                if offs is None or counts is None:
                    per_ifd.append((ifd, 'tile table unreadable')); total += 1
                else:
                    bad = sum(1 for i, c in enumerate(counts)
                              if (c > 0 and (offs[i] == 0 or offs[i] + c > size)) or c > MAX_TILE_BYTES)
                    unw = sum(1 for i, c in enumerate(counts) if c == 0 and offs[i] == 0)
                    tiles += len(offs); total += bad; unwritten += unw
                    if bad or unw:
                        per_ifd.append((ifd, bad + unw))
            offset = nxt
            ifd += 1
    return tiles, total, unwritten, per_ifd


def decode(path):
    """Read every block of IFD0 through GDAL; return the number of blocks that failed."""
    os.environ['GDAL_PAM_ENABLED'] = 'NO'      # read-only: no statistics are written back into the file
    import rasterio
    failed = 0
    with rasterio.open(path, 'r') as ds:
        for b in range(1, ds.count + 1):
            for _, w in ds.block_windows(b):
                try:
                    ds.read(b, window=w)
                except Exception:          # noqa: BLE001
                    failed += 1
    return failed


def main():
    damaged, checked = [], 0
    targets = [(os.path.join(OUT, p), OUT) for p in PATTERNS]
    if INPUTS:
        targets += [(os.path.join(ROOT, p), ROOT) for p in INPUT_PATTERNS]
    for pattern, base in targets:
        for path in sorted(glob.glob(pattern)):
            checked += 1
            rel = os.path.relpath(path, base)
            result = scan(path)
            if result is None:
                damaged.append(rel)
                print('UNPARSED %-58s not a TIFF' % rel)
                continue
            tiles, bad, unw, per_ifd = result
            failed = decode(path) if DECODE and not (bad or unw) else 0
            if bad or unw or failed:
                damaged.append(rel)
                print('DAMAGED  %-58s %d bad + %d unwritten of %d tiles, %d blocks failed to decode, IFDs %s'
                      % (rel, bad, unw, tiles, failed, per_ifd))
            else:
                print('ok       %-58s %d tiles%s' % (rel, tiles, ', all blocks decoded' if DECODE else ''))

    print()
    if not checked:
        print('no rasters found under %s' % os.path.relpath(OUT, ROOT))
        print('Nothing was verified. Run this where the delivered rasters are.')
        return 1
    if damaged:
        print('%d damaged file(s). Overviews may still be intact, so a GIS can open these normally '
              'and a header-only manifest will not notice.' % len(damaged))
        print('Re-write the affected layers before depositing.')
        return 1
    print('all files clean')
    return 0


if __name__ == '__main__':
    sys.exit(main())
