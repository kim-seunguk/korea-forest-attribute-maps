# 03a_sentinel2_composites.py — export the Sentinel-2 seasonal composites as a raster.
# =============================================================================================
# The 42 predictor bands, wall to wall, on the delivered grid. This is what the models predict onto;
# 03e samples the same composites at the plots, which is what they are trained on. Both take the
# image from s2_composite.py so that they cannot drift apart.
#
# The export is masked to the forest extent, which keeps it small enough to move. Earth Engine may
# still split it into several tiles.
#
# Out: Drive/<DRIVE_FOLDER>/*.tif — download into data/raw/gee_s2med/, then run 03b to put it on
#      the covariate grid.
# The mask is forest_base, the forest-type map's forest together with ESA tree cover where that map
# does not reach, so the base products cover the same cells the composite does.
# =============================================================================================
import os
import ee
from s2_composite import BAND_NAMES, seasonal_median

# Your own Earth-Engine-enabled Cloud project, and your own copy of the forest-type mask asset.
# Set EE_PROJECT before running, EE_FTM_ASSET if you uploaded the mask elsewhere, and
# EE_DRIVE_FOLDER if you want the export somewhere other than the default folder.
PROJECT = os.environ.get('EE_PROJECT', '')
if not PROJECT:
    raise SystemExit('set EE_PROJECT to your Earth Engine Cloud project id')
FTM_ASSET = os.environ.get('EE_FTM_ASSET', f'projects/{PROJECT}/assets/ftm2025_storunst_30m')
try:
    ee.Initialize(project=PROJECT)             # persistent credentials when run locally
except Exception:
    ee.Authenticate(); ee.Initialize(project=PROJECT)   # first run: opens an auth link

CRS, SCALE   = 'EPSG:5179', 30                          # the delivered 30 m product grid
DRIVE_FOLDER = os.environ.get('EE_DRIVE_FOLDER', 'nfimap_s2_export')
korea        = ee.Geometry.Rectangle([124.5, 33.0, 131.9, 38.9])   # S. Korea bbox (EPSG:4326); ymax 38.9 ≥ extended covar 38.8

# Create DRIVE_FOLDER in Drive yourself before starting the export, and leave it empty. The export
# splits into many tiles, and if the folder does not exist yet they race to create it, which can end
# with the tiles spread over several folders of the same name.

# ---- forest_base mask = FTM forest ∪ ESA tree cover where the map does not reach -------------------
# Must match the base direct map and the base imputation extent. A bare STORUNST == 1 mask stops at the
# forest-type map's own coverage, which ends around 38.37N and is also patchy elsewhere; forest_base
# fills what the map does not cover with the ESA WorldCover tree class. The northern strip is the most
# visible part of that, but not the largest: it is roughly a fifth of the added area.
# unmask(0, False) is essential: without it the mask stops at the asset footprint.
ftm        = ee.Image(FTM_ASSET)
ftm_forest = ftm.eq(1)                                              # authoritative (stocked) forest
coverage   = ftm.mask()                                            # FTM surveyed extent
sk_land    = ee.Image.constant(1).clipToCollection(               # SK land (NON-simplified GAUL: coast/MDL follow-through)
               ee.FeatureCollection('FAO/GAUL/2015/level0')
                 .filter(ee.Filter.eq('ADM0_NAME', 'Republic of Korea'))).mask()
esa        = ee.ImageCollection('ESA/WorldCover/v200').first().eq(10)   # global tree cover (gap fallback only)
gap        = coverage.unmask(0, False).Not().And(sk_land)          # South Korean land the FTM does not cover
forest_base = ftm_forest.unmask(0, False).Or(esa.And(gap)).selfMask()  # unmask(...,False) extends past the asset footprint

# ---- the composite, masked to the forest extent -------------------------------------------------
img = seasonal_median(korea).updateMask(forest_base).toFloat()

# ---- export ------------------------------------------------------------------------------------------
print('S2 export bands =', img.bandNames().size().getInfo(), '(expect %d)' % len(BAND_NAMES))
ee.batch.Export.image.toDrive(image=img, description='m8_s2med42_30m', folder=DRIVE_FOLDER,
    fileNamePrefix='m8_s2med42_30m', region=korea, scale=SCALE, crs=CRS, maxPixels=1e13).start()
print('STARTED  m8_s2med42_30m  ->  Drive/%s' % DRIVE_FOLDER)
print('Monitor: https://code.earthengine.google.com/tasks  (may split into -0000000000-*.tif tiles)')
