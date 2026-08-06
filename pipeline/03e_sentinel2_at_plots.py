# 03e_sentinel2_at_plots.py — sample the Sentinel-2 seasonal composites at the inventory plots.
# =============================================================================================
# Stage 03a exports the seasonal composites as a wall-to-wall raster, which is what the models
# predict onto. This samples the same composites at the plot locations, which is what the models
# are trained on. Both take the image from s2_composite.py, so the two cannot drift apart.
#
# Sampling happens on Earth Engine's servers; only the resulting table is downloaded.
#
# In : a CSV of plot coordinates you export from data/interim/nfi8/rd_points.csv:
#        fwrite(rd[, .(Plot_ID, UTMK_X, UTMK_Y)], "plot_coords.csv")
#      Plot coordinates are confidential and are not distributed with this repository.
# Out: Drive/<DRIVE_FOLDER>/plot_features_s2_NN.csv, one per chunk. Download them into
#        data/raw/gee_plot_features/ and run 03f to assemble the feature table.
#
# RUN  set EE_PROJECT to your own Earth-Engine-enabled Cloud project, then
#        python pipeline/03e_sentinel2_at_plots.py
#      EE_DRIVE_FOLDER_PLOTS chooses the Drive folder; PLOT_COORDS the coordinate CSV.
# =============================================================================================
import os
import ee
import pandas as pd
from pyproj import Transformer
from s2_composite import seasonal_median

PROJECT = os.environ.get('EE_PROJECT', '')
if not PROJECT:
    raise SystemExit('set EE_PROJECT to your Earth Engine Cloud project id')
try:
    ee.Initialize(project=PROJECT)                      # reuses cached credentials
except Exception:
    ee.Authenticate(); ee.Initialize(project=PROJECT)   # first run only: opens a browser

COORDS       = os.environ.get('PLOT_COORDS', 'plot_coords.csv')
DRIVE_FOLDER = os.environ.get('EE_DRIVE_FOLDER_PLOTS', 'nfimap_plot_features')
CRS          = 'EPSG:5179'
CHUNK        = 750                               # plots per export task


def to_feature_collection(df):
    feats = [ee.Feature(ee.Geometry.Point([r.lon, r.lat]), {'Plot_ID': r.Plot_ID})
             for r in df.itertuples()]
    return ee.FeatureCollection(feats)


def main():
    plots = pd.read_csv(COORDS, dtype={'Plot_ID': str})
    tf = Transformer.from_crs(CRS, 'EPSG:4326', always_xy=True)
    plots['lon'], plots['lat'] = tf.transform(plots['UTMK_X'].values, plots['UTMK_Y'].values)
    plots = plots.reset_index(drop=True)

    # Exported in chunks. One request carrying every plot as an inline feature collection is large
    # enough to stall on upload, while CHUNK plots at a time uploads immediately. The compute itself
    # runs server-side either way.
    started = 0
    for start in range(0, len(plots), CHUNK):
        part = plots.iloc[start:start + CHUNK]
        fc = to_feature_collection(part)
        sampled = (seasonal_median(fc)
                   .reduceRegions(collection=fc, reducer=ee.Reducer.mean(),
                                  scale=30, crs=CRS, tileScale=16)
                   .map(lambda f: f.setGeometry(None)))
        name = f'plot_features_s2_{start // CHUNK:02d}'
        ee.batch.Export.table.toDrive(collection=sampled, description=name,
                                      folder=DRIVE_FOLDER, fileNamePrefix=name,
                                      fileFormat='CSV').start()
        started += 1
        print('started', name, '-', len(part), 'plots')

    print(f'{started} export tasks started -> Drive/{DRIVE_FOLDER}/')
    print('When they finish, download them into data/raw/gee_plot_features/ and run 03f.')


if __name__ == '__main__':
    main()
