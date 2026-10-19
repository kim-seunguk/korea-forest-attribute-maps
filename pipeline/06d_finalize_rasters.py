# 06d_finalize_rasters.py — normalize the class maps, embed the species color table and (re)build
# overviews.
# =============================================================================================
# Run AFTER 05c_index_to_attributes.R + 06b_apply_landcover_filter.R (the filter drops overviews).
# Three jobs, all operating IN PLACE on outputs/nfi8/:
#
#   1) NORMALIZE CLASS MAPS: small-code class rasters (dominant species, V_*) that are still float32/NaN
#      -> uint8/nodata=255. 05a already writes the dominant-species map as uint8, so it is skipped here. The
#      decoded attribute rasters from 05c are the ones converted. An integer dtype with an integer nodata
#      suits a discrete class raster AND makes MODE overviews safe (step 3).
#
#   2) COLOR TABLE on the dominant-species map, from data/species_colors.csv — the palette the figures use,
#      the named species in their own color and the rest gray — so a GIS colors the map on open.
#
#   3) INTERNAL overviews (pyramids) on every delivered raster the filter touched (levels 2/4/8/16/32):
#        - CLASS maps (dominant species, V_*, now uint8/255) -> MODE: fills the footprint as AVERAGE does for
#          the continuous layers, so the national zoom-out does not look sparse. MODE is safe here because
#          the nodata is an INTEGER. With a floating-point NaN nodata, GDAL's mode fills all-NoData blocks such
#          as the ocean with class code 0, which is why step 1 exists.
#        - imputation plot-index (codes > 255, still integer) -> NEAREST (mode is meaningless on plot indices).
#        - CONTINUOUS rasters (coniferous fraction, biomass, carbon, Shannon, BA share, SD/uncertainty) -> AVERAGE.
#      The overviews are internal (COG-style), the delivery format: each raster stays a single file, and
#      ArcGIS and QGIS read the overviews from it.
#
#   python pipeline/06d_finalize_rasters.py
# =============================================================================================
import csv
import os, glob, re, sys, time
import numpy as np
import rasterio
from rasterio.enums import Resampling

VERSION = "nfi8"
PREDSET = "ftm"             # the delivered file prefix (R/paths.R::FILE_PREDSET)
OVR_LEVELS = [2, 4, 8, 16, 32]
# Discrete rasters, matched by file-name token (overview rules in the header). CLASS_UINT8 maps
#   become uint8/255 with MODE overviews, the rest of CATEGORICAL (the plot index) gets NEAREST, since
#   averaging invents codes. "dominant_species" is the published name (R/paths.R), "dom_species" the
#   decoded attribute.
CLASS_UINT8 = ("dominant_species", "dom_species", "spec30m_v_vg", "spec30m_v_dc", "spec30m_v_dn")
CATEGORICAL = CLASS_UINT8 + ("imputation",)
# ⚠ Uncertainty layers can carry a class token ('dominant_species_uncertainty') but are continuous
#   (0..1 vote margin, SD), so they are NEVER class-treated (uint8/mode/colormap).
UNCERTAINTY = ("_sd_", "uncertainty")
def _is_uncertainty(bn):
    return any(t in bn.lower() for t in UNCERTAINTY)

# The names the pipeline publishes, mirroring R/paths.R::product_file and what 05b and 05c write.
# This stage is destructive — it rewrites dtype and rebuilds overviews in place — so it touches only
# the files the pipeline writes. Any other file in these directories is left alone.
DELIVERED = re.compile(
    r'^(?:'
    r'ftm_(?:conifer_fraction|biomass|shannon|dominant_species)(?:_sd|_uncertainty)?_30m\.tif'
    r'|imputation_ftm\.tif'
    r'|nearest_reference_distance_ftm\.tif'
    r'|spec30m_(?:V_VG|V_DC|V_DN|dom_species|dom_ba_share'
    r'|Bio_a_ha|Carbon_a_ha|BAPH|shannon|coniferfrac)\.tif'
    r')$')

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT  = os.path.join(ROOT, "outputs", VERSION)
PALETTE = os.path.join(ROOT, "data", "species_colors.csv")


def _display_files():
    files = (sorted(glob.glob(os.path.join(OUT, "maps_30m", "*.tif"))) +
             sorted(glob.glob(os.path.join(OUT, "uncertainty_30m", "*.tif"))) +
             [os.path.join(OUT, "index", "imputation_%s.tif" % PREDSET),
              # delivered alongside the index, so it needs the same overviews
              os.path.join(OUT, "index", "nearest_reference_distance_%s.tif" % PREDSET)] +
             sorted(glob.glob(os.path.join(OUT, "attributes_%s" % PREDSET, "spec30m_*.tif"))))
    files = [f for f in files
             if os.path.exists(f) and DELIVERED.match(os.path.basename(f))]
    # Optional scope limiter, mirroring LC_ONLY in 06b_apply_landcover_filter.R. Building internal
    # overviews rewrites the file tail, where a tile index can be damaged, so a run that replaced a
    # few layers need not rewrite the rest.
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
    """Convert small-code class maps (dominant species, V_*) from float32/NaN to uint8/255, in place and
    idempotent (header step 1). Leaves the plot index and all continuous maps untouched. DESTRUCTIVE:
    archive the float32 originals first if you want them."""
    conv = skip = 0
    for f in _display_files():
        bn = os.path.basename(f).lower()
        if _is_uncertainty(bn):
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


def embed_species_colormap():
    """Bake the species symbology into the uint8 dominant-species GeoTIFF as an internal color table, so a GIS
    colors the map on open with the same palette the figures use. The table in data/species_colors.csv indexes
    by class code and covers every legend code. Codes left unused by the trained legend have no cells.
    Must run AFTER normalize_class_maps() (colour tables need Byte/uint8). Idempotent."""
    if not os.path.exists(PALETTE):
        print("  (skip colormap: %s not found)" % os.path.relpath(PALETTE, ROOT)); return
    cmap = {}
    with open(PALETTE, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            h = r["hex"].lstrip("#")
            cmap[int(r["code"])] = (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16), 255)
    cmap[255] = (0, 0, 0, 0)                                   # nodata -> transparent
    n = 0
    for f in _display_files():
        bn = os.path.basename(f).lower()
        if "dominant_species" not in bn or _is_uncertainty(bn):   # the species MAP only
            continue
        try:
            with rasterio.open(f, "r+") as ds:
                if ds.dtypes[0] != "uint8":
                    print("  colormap SKIP %s: not uint8 (run normalize first)" % os.path.basename(f)); continue
                ds.write_colormap(1, cmap)
            print("  colormap ->", os.path.basename(f)); n += 1
        except Exception as e:                                # noqa: BLE001
            print("  colormap SKIP %s: %r" % (os.path.basename(f), e))
    print("species colormap embedded on %d raster(s)" % n)


def build_overviews():
    files = _display_files()
    print("building internal overviews on %d rasters (levels %s):" % (len(files), OVR_LEVELS))
    ok = skip = fail = 0
    for f in files:
        bn = os.path.basename(f)
        low = bn.lower()
        if _is_uncertainty(low):                     # SD and vote-margin layers are CONTINUOUS -> average (check FIRST)
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
    print("== 1. normalize class maps (dominant species, V_*) -> uint8/255 ==")
    normalize_class_maps()
    print("\n== 2. embed the species colour table ==")
    embed_species_colormap()
    print("\n== 3. (re)build internal overviews ==")
    build_overviews()
