"""s2_composite.py — the Sentinel-2 seasonal composite, defined once.

Two stages need this image and they must agree exactly. 03a exports it wall to wall, which is what
stage 05 predicts onto; 03e samples it at the inventory plots, which is what the models are trained
on. If the two ever described the bands differently the models would be trained on one thing and
applied to another, and nothing in the outputs would say so. So both build the image from here
rather than each carrying its own copy.

This is a shared definition, not a pipeline stage — it reads nothing and writes nothing.

The band names it produces are season-major: spring's fourteen, then summer's, then fall's. That is
the physical order of the exported raster, and R/… stages label the raster with it. The plot table
that 03e exports is not in that order, because Earth Engine's table export sorts its properties
alphabetically; downstream code selects those columns by name.
"""
import ee

# Surface-reflectance bands, then the indices computed from them.
BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']
INDICES = ['NDVI', 'NDWI', 'NDRE', 'EVI']
S2_BANDS = BANDS + INDICES

# NDWI here is the NIR-SWIR1 moisture index (B8, B11) = NDMI (Wilson & Sader 2002). It is NOT Gao
# (1996) NDWI, which pairs 0.86 um with 1.24 um rather than B11 (1.61 um). The internal label stays
# NDWI to match the existing feature tables and models; it is reported and deposited as NDMI.

SEASONS = {'spring': ('-04-01', '-05-31'),
           'summer': ('-06-01', '-08-31'),
           'fall':   ('-09-01', '-10-31')}
YEARS = [2021, 2022, 2023, 2024, 2025]      # the NFI-8 measurement period

# Cloud Score+ threshold. The monsoon leaves few clear scenes, and a stricter cut empties the summer
# composite over parts of the country.
CS_THRESHOLD = 0.60

BAND_NAMES = [f's2_{season}_{band}' for season in SEASONS for band in S2_BANDS]   # 42, season-major


def add_indices(img):
    """The four indices, computed on the seasonal median rather than per scene, so that how many
    scenes a pixel kept after cloud masking cannot change what an index means between seasons."""
    return img.addBands([
        img.normalizedDifference(['B8', 'B4']).rename('NDVI'),
        img.normalizedDifference(['B8', 'B11']).rename('NDWI'),
        img.normalizedDifference(['B8', 'B5']).rename('NDRE'),
        img.expression('2.5*((N-R)/(N+6*R-7.5*B+1))',
                       {'N': img.select('B8'), 'R': img.select('B4'),
                        'B': img.select('B2')}).rename('EVI')])


def seasonal_median(region):
    """One image of 42 named bands: per season, the median over every cloud-free observation of
    every year in the window. `region` may be a geometry or a feature collection."""
    csp = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
    seasons = []
    for season, (start, end) in SEASONS.items():
        col = ee.ImageCollection([])
        for year in YEARS:
            col = col.merge(
                ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterBounds(region)
                  .filterDate(f'{year}{start}', f'{year}{end}')
                  .linkCollection(csp, ['cs'])
                  .map(lambda img: img.updateMask(img.select('cs').gte(CS_THRESHOLD)))
                  .select(BANDS))
        seasons.append(add_indices(col.median())
                       .rename([f's2_{season}_{band}' for band in S2_BANDS]))
    return ee.Image.cat(seasons)
