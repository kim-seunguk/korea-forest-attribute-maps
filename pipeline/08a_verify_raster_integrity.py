"""08a_verify_raster_integrity.py — pre-deposit integrity check for every delivered GeoTIFF.

WHY: base-tier rasters have been found with corrupted TileOffsets / TileByteCounts
in IFD0 while their overviews stayed intact, so ArcGIS opened them normally and a header-only check
(08c_build_data_records.py reads the profile, never the pixels) reported them as fine. Sampling a few
windows is not enough either — one of the three passed a three-window read test.

This script parses the TIFF tile tables directly, for every IFD, and flags any tile whose byte range
cannot exist inside the file. It decodes no pixels, so it runs in seconds over 15 GB.

RUN  python pipeline/08a_verify_raster_integrity.py
EXIT 0 = all clean, 1 = at least one damaged file (so it can gate a deposit script).
"""
import glob
import os
import struct
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = os.path.join(ROOT, 'outputs', 'nfi8')
PATTERNS = ['m8_maps_30m/*.tif', 'm8_uncertainty_30m/*.tif', 'index/*.tif',
            'attributes_ftm/*.tif', 'attributes_base/*.tif']
TYPE_SIZE = {1: 1, 3: 2, 4: 4}
FMT = {1: 'B', 3: 'H', 4: 'I'}
MAX_TILE_BYTES = 50_000_000          # a 256x256 float64 tile is ~0.5 MB; anything near this is garbage


def scan(path):
    """Return (n_tiles, n_bad, [(ifd, n_bad), ...]), or a string saying why it was not checked."""
    size = os.path.getsize(path)
    with open(path, 'rb') as f:
        head = f.read(4)
        if head[:2] not in (b'II', b'MM'):
            return 'not a TIFF'
        endian = '<' if head[:2] == b'II' else '>'
        if struct.unpack(endian + 'H', head[2:4])[0] == 43:
            # BigTIFF. Several writers use BIGTIFF=IF_SAFER, so one can appear once a layer passes
            # 4 GB, and this parser only understands the classic tile tables.
            return 'BigTIFF: this parser reads classic TIFF only'
        f.seek(4)
        offset = struct.unpack(endian + 'I', f.read(4))[0]

        per_ifd, total, tiles = [], 0, 0
        ifd = 0
        while offset and 0 < offset < size:
            f.seek(offset)
            n_entries = struct.unpack(endian + 'H', f.read(2))[0]
            if not 0 < n_entries <= 200:
                break
            tags = {}
            for _ in range(n_entries):
                entry = f.read(12)
                tag, typ, count = struct.unpack(endian + 'HHI', entry[:8])
                tags[tag] = (typ, count, entry[8:12])
            nxt = struct.unpack(endian + 'I', f.read(4))[0]

            def array(tag):
                typ, count, val = tags[tag]
                nbytes = count * TYPE_SIZE.get(typ, 4)
                if nbytes <= 4:
                    data = val[:nbytes]
                else:
                    ptr = struct.unpack(endian + 'I', val)[0]
                    if ptr + nbytes > size:
                        return None
                    f.seek(ptr)
                    data = f.read(nbytes)
                return list(struct.unpack(endian + str(count) + FMT[typ], data))

            if 324 in tags and 325 in tags:
                offs, counts = array(324), array(325)
                if offs and counts:
                    bad = sum(1 for i, c in enumerate(counts)
                              if (c > 0 and (offs[i] == 0 or offs[i] + c > size)) or c > MAX_TILE_BYTES)
                    tiles += len(offs)
                    total += bad
                    if bad:
                        per_ifd.append((ifd, bad))
            offset = nxt
            ifd += 1
    return tiles, total, per_ifd


def main():
    damaged, skipped, checked = [], [], 0
    for pattern in PATTERNS:
        for path in sorted(glob.glob(os.path.join(OUT, pattern))):
            checked += 1
            result = scan(path)
            rel = os.path.relpath(path, OUT)
            if isinstance(result, str):
                skipped.append((rel, result))
                print('SKIPPED  %-58s %s' % (rel, result))
                continue
            tiles, bad, per_ifd = result
            if bad:
                damaged.append((rel, tiles, bad, per_ifd))
                print('DAMAGED  %-58s %d of %d tiles unreadable, IFDs %s'
                      % (rel, bad, tiles, per_ifd))
            else:
                print('ok       %-58s %d tiles' % (rel, tiles))

    print()
    if not checked:
        # Finding nothing is not the same as finding nothing wrong, and this gates a deposit.
        print('no rasters found under %s' % os.path.relpath(OUT, ROOT))
        print('Nothing was verified. Run this where the delivered rasters are.')
        return 1
    if damaged:
        print('%d damaged file(s). Overviews may still be intact, so a GIS can open these normally '
              'and a header-only manifest will not notice.' % len(damaged))
        print('Re-write the affected layers before depositing.')
    if skipped:
        # A file this cannot parse is not a file it found clean. Saying so is the whole point: the
        # damage this catches was invisible to every check that stopped short of the tile tables.
        print('%d file(s) were not checked. They are not known to be clean.' % len(skipped))
    if damaged or skipped:
        return 1
    print('all files clean')
    return 0


if __name__ == '__main__':
    sys.exit(main())
