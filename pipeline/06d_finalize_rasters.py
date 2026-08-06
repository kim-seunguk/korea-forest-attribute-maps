# 06d_finalize_rasters.py — regenerate the int8 Shannon display raster + (re)build overviews.
# =============================================================================================
# Run AFTER 05c_index_to_attributes.R + 06b_apply_landcover_filter.R (the filter drops overviews and the
# compact Shannon raster goes stale). Three jobs, all operating IN PLACE on outputs/nfi8/:
#
#   1) SHANNON int8 display raster  m8_maps_30m/<tier>_shannon_30m_int.tif
#      = round(H' x 20) as INT1S, NoData = -128 (H' in [0, ~1.9] -> [0, ~38], fits int8; 0.05 H' steps).
#        This is a COMPACT display/symbology raster for ArcGIS (NOT a delivered data layer — the accuracy table
#        ships Shannon as float). Stale ArcGIS sidecars (.tfw/.vat.dbf/.vat.cpg/.xml) are removed so the
#        renderer rebuilds them against the new values.
#
#   2) NORMALIZE CLASS MAPS: small-code class rasters (domsp, V_*) exported as float32/NaN -> uint8/nodata=255.
#      Integer nodata is the correct dtype for a discrete class raster AND unlocks safe MODE overviews (step 3).
#
#   3) INTERNAL overviews (pyramids) on every delivered raster the filter touched (levels 2/4/8/16/32):
#        - CLASS maps (domsp, V_*, now uint8/255) -> MODE: fills the footprint like AVERAGE does for the
#          continuous panels, so the national zoom-out matches instead of looking sparse. Safe because the
#          nodata is an INTEGER. With a floating-point NaN nodata, GDAL's mode fills all-NoData
#          blocks such as the ocean with class code 0, which is why step 2 exists.
#        - imputation plot-index (codes > 255, still integer) -> NEAREST (mode is meaningless on plot indices).
#        - CONTINUOUS rasters (coniferous fraction, biomass, carbon, Shannon, BA share, SD/margin) -> AVERAGE.
#      Internal (COG-style) overviews: single file, read by ArcGIS/QGIS, and the delivery format.
#
#   python pipeline/06d_finalize_rasters.py
# =============================================================================================
import os, glob, re, sys, time
import numpy as np
import rasterio
from rasterio.enums import Resampling

VERSION = "nfi8"
# Which tier gets the int8 Shannon display raster. It is derived from that tier's map, so name the
# one you have just rebuilt: DELIVER_TIER=ftm python pipeline/06d_finalize_rasters.py
DELIVER_TIER = os.environ.get("DELIVER_TIER", "base")
SHANNON_SCALE = 20          # H' x 20 -> int8 (0.05 H' resolution; max ~1.9 -> ~38, within [-127,127])
OVR_LEVELS = [2, 4, 8, 16, 32]
# categorical = discrete codes, matched by these tokens in the file name.
#   CLASS_UINT8 = small-code CLASS maps (domsp, V_*) -> normalized to uint8/nodata=255 -> MODE overviews
#                 (mode fills the footprint like average for the continuous panels; with an INTEGER nodata it is
#                  safe, whereas on float32/NaN GDAL mode floods all-NoData blocks (the ocean) with code 0).
#   CATEGORICAL = the rest of the discrete rasters (imputation plot-index, codes > 255) -> NEAREST (mode is
#                 meaningless on arbitrary plot indices; averaging invents non-existent codes).
CLASS_UINT8 = ("domsp", "dom_species", "spec30m_v_vg", "spec30m_v_dc", "spec30m_v_dn")
CATEGORICAL = CLASS_UINT8 + ("imputation",)
# ⚠ CONTINUOUS uncertainty layers whose NAME contains a class token ('domsp' in 'domsp_margin') — they must
#   NEVER be class-treated (uint8/mode/colormap): the vote-margin & SD are 0..1 / continuous. Guard on this.
UNCERTAINTY = ("_sd_", "margin")
def _is_uncertainty(bn):
    return any(t in bn.lower() for t in UNCERTAINTY)

# The names the pipeline publishes, mirroring R/paths.R::product_file and what 05b and 05c write.
# This stage is destructive — it rewrites dtype and rebuilds overviews in place — so it must only
# touch what it produced. A file that merely sits in the same directory is something else's: a
# display variant with a colour table already baked in, an experiment, a copy kept for a figure.
DELIVERED = re.compile(
    r'^(?:'
    r'(?:base|ftm)_(?:conifer_fraction|biomass|shannon|dominant_species)(?:_sd|_uncertainty)?_30m\.tif'
    r'|(?:base|ftm)_shannon_30m_int\.tif'
    r'|imputation_(?:base|ftm)\.tif'
    r'|nearest_reference_distance_(?:base|ftm)\.tif'
    r'|spec30m_(?:V_VG|V_DC|V_DN|dom_species|dom_ba_share'
    r'|Bio_a_ha|Carbon_a_ha|BAPH|shannon|coniferfrac)\.tif'
    r')$')

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT  = os.path.join(ROOT, "outputs", VERSION)


def _display_files():
    files = (sorted(glob.glob(os.path.join(OUT, "m8_maps_30m", "*.tif"))) +
             sorted(glob.glob(os.path.join(OUT, "m8_uncertainty_30m", "*.tif"))) +
             [os.path.join(OUT, "index", "imputation_ftm.tif"),
              os.path.join(OUT, "index", "imputation_base.tif"),
              # delivered alongside the index, so they need the same overviews
              os.path.join(OUT, "index", "nearest_reference_distance_ftm.tif"),
              os.path.join(OUT, "index", "nearest_reference_distance_base.tif")] +
             sorted(glob.glob(os.path.join(OUT, "attributes_ftm", "spec30m_*.tif"))) +
             sorted(glob.glob(os.path.join(OUT, "attributes_base", "spec30m_*.tif"))))
    files = [f for f in files
             if os.path.exists(f) and DELIVERED.match(os.path.basename(f))]
    # Optional scope limiter for repair runs, mirroring LC_ONLY in 06b_apply_landcover_filter.R.
    # Building internal overviews rewrites the file tail, which is where the tile-index damage this
    # dataset has seen appeared, so a run that replaced a few layers should not rewrite the rest.
    # Name the files to limit it to, as arguments or in FINALIZE_ONLY. Arguments are preferable:
    # an environment variable left set will silently narrow a later run.
    #   python pipeline/06d_finalize_rasters.py a.tif b.tif
    only = ",".join(sys.argv[1:]) or os.environ.get("FINALIZE_ONLY", "").strip()
    if only:
        keep = {s.strip() for s in only.split(",") if s.strip()}
        files = [f for f in files if os.path.basename(f) in keep]
        missing = keep - {os.path.basename(f) for f in files}
        if missing:
            raise SystemExit("FINALIZE_ONLY names no such delivered raster: %s" % ", ".join(sorted(missing)))
        print("FINALIZE_ONLY: %d file(s)" % len(files))
    return files


def normalize_class_maps():
    """Small-code CLASS maps (domsp, V_*) -> uint8 / nodata=255, idempotent & in-place. They are exported as
    float32/NaN; integer nodata (a) makes MODE overviews fill the footprint safely (GDAL mode floods NaN blocks
    with 0) and (b) is the right dtype for a discrete class raster (¼ the size). Leaves the imputation
    plot-index (codes > 255) and all continuous maps untouched. DESTRUCTIVE: if you want the float32
    originals, archive them before running this."""
    conv = skip = 0
    for f in _display_files():
        bn = os.path.basename(f).lower()
        if _is_uncertainty(bn):                    # domsp_margin / *_sd are CONTINUOUS despite the class token
            continue
        if not any(tok in bn for tok in CLASS_UINT8):
            continue
        with rasterio.open(f) as s:
            if s.dtypes[0] == "uint8":
                skip += 1; continue
            a = s.read(1); prof = s.profile.copy()
        mx = np.nanmax(a)
        if not np.isfinite(mx) or mx > 254:
            print("  SKIP (codes exceed uint8 or empty):", bn); skip += 1; continue
        b = np.where(np.isnan(a), 255, np.rint(a)).astype("uint8")
        prof.update(dtype="uint8", nodata=255, compress="LZW", predictor=2,
                    tiled=True, blockxsize=256, blockysize=256)
        prof.pop("photometric", None)
        tmp = f + ".tmp.tif"
        with rasterio.open(tmp, "w", **prof) as o:
            o.write(b, 1)
        os.replace(tmp, f)
        print("  -> uint8/255:", bn); conv += 1
    print("class-map normalize: %d converted / %d already-ok" % (conv, skip))


def regenerate_shannon_int():
    src = os.path.join(OUT, "m8_maps_30m", "%s_shannon_30m.tif" % DELIVER_TIER)   # the tier set above
    dst = os.path.join(OUT, "m8_maps_30m", "%s_shannon_30m_int.tif" % DELIVER_TIER)
    if not os.path.exists(src):
        print("  (skip shannon_int: %s not found)" % os.path.basename(src)); return None
    for side in (dst + ".vat.dbf", dst + ".vat.cpg", dst + ".xml", dst + ".aux.xml",
                 dst + ".ovr", dst[:-4] + ".tfw"):
        if os.path.exists(side):
            os.remove(side); print("  removed stale sidecar:", os.path.basename(side))
    with rasterio.open(src) as s:
        prof = s.profile.copy()
        prof.update(driver="GTiff", dtype="int8", nodata=-128, count=1,
                    compress="LZW", tiled=True, blockxsize=256, blockysize=256)
        prof.pop("photometric", None)
        mn, mx = 127, -128
        with rasterio.open(dst, "w", **prof) as d:
            for _, w in s.block_windows(1):
                a = s.read(1, window=w)
                out = np.full(a.shape, -128, dtype="int8")
                m = np.isfinite(a)
                if m.any():
                    sv = np.clip(np.round(a[m] * SHANNON_SCALE), -127, 127).astype("int8")
                    out[m] = sv; mn = min(mn, int(sv.min())); mx = max(mx, int(sv.max()))
                d.write(out, 1, window=w)
    print("  shannon_int: round(H' x %d) int8 NoData=-128  range %d..%d" % (SHANNON_SCALE, mn, mx))
    return dst


def build_overviews():
    files = _display_files()
    print("building internal overviews on %d rasters (levels %s):" % (len(files), OVR_LEVELS))
    ok = skip = fail = 0
    for f in files:
        bn = os.path.basename(f)
        low = bn.lower()
        # CLASS_UINT8 (domsp, V_*, now uint8/255) -> MODE: fills the footprint like average does for the
        #   continuous panels, so the national zoom-out matches (safe: integer nodata, no ocean-0 flood).
        # other CATEGORICAL (imputation plot-index) -> NEAREST. everything else (continuous) -> AVERAGE.
        if _is_uncertainty(low):                     # domsp_margin / *_sd are CONTINUOUS -> average (check FIRST)
            rs, tag = Resampling.average, "avg "
        elif any(tok in low for tok in CLASS_UINT8):
            rs, tag = Resampling.mode, "mode"
        elif any(tok in low for tok in CATEGORICAL):
            rs, tag = Resampling.nearest, "near"
        else:
            rs, tag = Resampling.average, "avg "
        try:
            t0 = time.time()
            with rasterio.Env(COMPRESS_OVERVIEW="LZW", GDAL_TIFF_OVR_BLOCKSIZE="256"):
                with rasterio.open(f, "r+") as ds:
                    ds.build_overviews(OVR_LEVELS, rs)
            print("  OK [%s] %-38s %5.1fs" % (tag, bn, time.time() - t0)); ok += 1
        except rasterio.errors.RasterioIOError as e:
            print("  SKIP (could not open for writing) %s: %s" % (bn, e)); skip += 1
        except Exception as e:                                    # noqa: BLE001
            print("  FAIL %s: %r" % (bn, e)); fail += 1
    print("DONE: %d built / %d skipped / %d failed" % (ok, skip, fail))


if __name__ == "__main__":
    print("== 1. regenerate Shannon int8 display raster ==")
    regenerate_shannon_int()
    print("\n== 2. normalize class maps (domsp, V_*) -> uint8/255 ==")
    normalize_class_maps()
    print("\n== 3. (re)build internal overviews ==")
    build_overviews()
