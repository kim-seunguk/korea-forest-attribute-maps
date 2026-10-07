# 03a_sentinel2_composites.py — export the Sentinel-2 seasonal composites onto the product grid.
# =============================================================================================
# The 42 predictor bands, wall to wall, on the delivered grid. Stage 05 predicts onto this raster,
# and once 03b and 03c have put it on disk, 03e reads the training features at the plots from it,
# so training and prediction see one set of numbers. The composite itself is defined in
# s2_composite.py.
#
# The composite has no native projection, so it is given a 10 m one and then averaged into each
# 30 m cell of the product grid (reduceResolution, mean): a cell holds the area mean of the seasonal
# medians inside it, and that mean is the only resampling in the whole chain. The export is written
# with the grid's own transform and dimensions, so the tiles land on the covariate grid cell for
# cell and 03b has nothing to resample.
#
# The export is masked to the forest-type map's stocked forest, which is the extent of the delivered
# products and keeps the export small enough to move. Earth Engine may still split it into tiles.
#
# Out: Drive/<DRIVE_FOLDER>/m8_s2med42_grid30m*.tif — download into data/raw/gee_s2med/, then run
#      03b to mosaic the tiles and label the bands.
#
# RUN  set EE_PROJECT to your own Earth-Engine-enabled Cloud project, then
#        python pipeline/03a_sentinel2_composites.py
#      EE_FTM_ASSET names your copy of the 02d mask asset, EE_DRIVE_FOLDER the Drive folder.
#      EE_DRY_RUN=1 evaluates the image over a small test box on the export grid and starts nothing.
# =============================================================================================
import os
import ee
from s2_composite import BAND_NAMES, seasonal_median

PROJECT = os.environ.get('EE_PROJECT', '')
if not PROJECT:
    raise SystemExit('set EE_PROJECT to your Earth Engine Cloud project id')
FTM_ASSET = os.environ.get('EE_FTM_ASSET', f'projects/{PROJECT}/assets/ftm2025_storunst_30m')
try:
    ee.Initialize(project=PROJECT)             # persistent credentials when run locally
except Exception:
    ee.Authenticate(); ee.Initialize(project=PROJECT)   # first run: opens an auth link

# The product grid: EPSG:5179, the covariate grid of 02c. The transform is column width, 0, west
# edge, 0, negative row height, north edge; the dimensions are columns x rows. 02f carries the same
# three lines, and 02c prints the template they must match.
GRID_CRS       = 'EPSG:5179'
GRID_TRANSFORM = [30.000454545454545, 0, 746150.0931818181, 0, -30.0, 2089170.0]
GRID_DIMS      = '18593x21013'
DRIVE_FOLDER   = os.environ.get('EE_DRIVE_FOLDER', 'nfimap_s2_export')
korea          = ee.Geometry.Rectangle([124.5, 33.0, 131.9, 38.9])   # the scene filter; the grid bounds the export itself

# Create DRIVE_FOLDER in Drive yourself before starting the export, and leave it empty. The export
# splits into many tiles, and if the folder does not exist yet they race to create it, which can end
# with the tiles spread over several folders of the same name.

# ---- the mask: the forest-type map's stocked forest ---------------------------------------------
# The 02d mask asset is 1 on stocked forest. Everything else — 0 inside the asset, and the area
# outside its footprint — is masked out of the composite.
ftm_forest = ee.Image(FTM_ASSET).eq(1)

# ---- the composite, averaged into the product grid cells and masked ----------------------------
s2 = seasonal_median(korea)
s2 = (s2.setDefaultProjection(ee.Projection(GRID_CRS).atScale(10))   # the median has no projection of its own
        .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=64)
        .reproject(crs=GRID_CRS, crsTransform=GRID_TRANSFORM))
img = s2.updateMask(ftm_forest).toFloat()

if os.environ.get('EE_DRY_RUN'):
    box = ee.Geometry.Rectangle([128.48, 37.48, 128.52, 37.52])   # forest near Pyeongchang
    v = img.reduceRegion(ee.Reducer.mean(), box, crs=GRID_CRS, crsTransform=GRID_TRANSFORM, maxPixels=1e7).getInfo()
    print('DRY RUN: %d bands evaluated over the test box; sample:' % len(v), dict(list(v.items())[:6]))
    print('DRY RUN: projection', img.projection().getInfo())
    raise SystemExit(0)

# ---- export ------------------------------------------------------------------------------------------
print('S2 export bands =', img.bandNames().size().getInfo(), '(expect %d)' % len(BAND_NAMES))
name = 'm8_s2med42_grid30m'
ee.batch.Export.image.toDrive(image=img, description=name, folder=DRIVE_FOLDER, fileNamePrefix=name,
    crs=GRID_CRS, crsTransform=GRID_TRANSFORM, dimensions=GRID_DIMS, maxPixels=1e13).start()
print('STARTED  %s  ->  Drive/%s' % (name, DRIVE_FOLDER))
print('Monitor: https://code.earthengine.google.com/tasks  (may split into -0000000000-*.tif tiles)')
