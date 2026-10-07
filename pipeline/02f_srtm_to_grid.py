# 02f_srtm_to_grid.py — export SRTM elevation onto the product grid (Earth Engine).
# =============================================================================================
# The terrain covariates (02c) start from an elevation raster that is already on the product grid.
# This export samples the native 1" SRTM bilinearly, once, directly onto that grid, so 02c derives
# slope and aspect without reprojecting anything. A DEM resampled onto a geographic grid and
# projected afterwards would carry duplicated rows and columns into slope and aspect as stripes.
#
# The grid is the covariate grid 02c builds from the forest-type map: its origin, cell size and
# dimensions are written out below, and 02c checks the downloaded tiles against the template it
# builds. Until the tiles are in place, 02c builds the climate layers and stops at the terrain
# layers, so it runs once before this export and once after the download.
#
# Out: Drive/<DRIVE_FOLDER>/srtm_elev_grid30m*.tif — download into data/raw/dem_srtm_grid/, then run
#      02c, with the four terrain layers in REBUILD if an earlier DEM has already produced them.
#
# RUN  set EE_PROJECT to your own Earth-Engine-enabled Cloud project, then
#        python pipeline/02f_srtm_to_grid.py
#      EE_DRIVE_FOLDER chooses the Drive folder. EE_DRY_RUN=1 evaluates the image over a small test
#      box on the export grid and starts nothing.
# =============================================================================================
import os
import ee

PROJECT = os.environ.get('EE_PROJECT', '')
if not PROJECT:
    raise SystemExit('set EE_PROJECT to your Earth Engine Cloud project id')
try:
    ee.Initialize(project=PROJECT)             # persistent credentials when run locally
except Exception:
    ee.Authenticate(); ee.Initialize(project=PROJECT)   # first run: opens an auth link

# The product grid: EPSG:5179, the covariate grid of 02c. The transform is column width, 0, west
# edge, 0, negative row height, north edge; the dimensions are columns x rows. 03a carries the same
# three lines, and 02c prints the template they must match.
GRID_CRS       = 'EPSG:5179'
GRID_TRANSFORM = [30.000454545454545, 0, 746150.0931818181, 0, -30.0, 2089170.0]
GRID_DIMS      = '18593x21013'
DRIVE_FOLDER   = os.environ.get('EE_DRIVE_FOLDER', 'nfimap_srtm_export')

# Create DRIVE_FOLDER in Drive yourself before starting the export, and leave it empty. If Earth
# Engine splits the export into tiles and the folder does not exist yet, the tiles race to create
# it and can end up spread over several folders of the same name.

dem = ee.Image('USGS/SRTMGL1_003').select('elevation').toFloat().resample('bilinear')

if os.environ.get('EE_DRY_RUN'):
    box = ee.Geometry.Rectangle([128.48, 37.48, 128.52, 37.52])   # forest near Pyeongchang
    v = dem.reduceRegion(ee.Reducer.mean(), box, crs=GRID_CRS, crsTransform=GRID_TRANSFORM, maxPixels=1e7).getInfo()
    print('DRY RUN: mean elevation over the test box', v)
    print('DRY RUN: projection', dem.projection().getInfo())
    raise SystemExit(0)

name = 'srtm_elev_grid30m'
ee.batch.Export.image.toDrive(image=dem, description=name, folder=DRIVE_FOLDER, fileNamePrefix=name,
    crs=GRID_CRS, crsTransform=GRID_TRANSFORM, dimensions=GRID_DIMS, maxPixels=1e13).start()
print('STARTED  %s  ->  Drive/%s' % (name, DRIVE_FOLDER))
print('Monitor: https://code.earthengine.google.com/tasks  (may split into -0000000000-*.tif tiles)')
