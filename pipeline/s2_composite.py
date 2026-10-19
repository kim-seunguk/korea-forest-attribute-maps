"""s2_composite.py — the Sentinel-2 seasonal composite, defined once.

03a exports this image onto the product grid. Stage 05 predicts onto that raster, and 03e samples it
at the plots for training. Nothing else on the Earth Engine side carries its own copy of the band
definitions.

This is a shared definition, not a pipeline stage — it reads nothing and writes nothing.

The band names it produces are season-major: spring's fourteen, then summer's, then fall's. That is
the physical order of the exported raster. R/paths.R holds the same names for the R side, which
labels the raster with them.
"""
import ee

# Surface-reflectance bands, then the indices computed from them.
BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']
INDICES = ['NDVI', 'NDWI', 'NDRE', 'EVI']
S2_BANDS = BANDS + INDICES

# NDWI here is the NIR-SWIR1 moisture index (B8, B11) = NDMI (Wilson & Sader 2002). It is NOT Gao
# (1996) NDWI, which pairs 0.86 um with 1.24 um rather than B11 (1.61 um). The band is labeled
# NDWI in the pipeline and published as NDMI.

# Season windows. filterDate's end date is exclusive, so each end is the first day after the season:
# spring is April and May, summer June to August, fall September and October.
SEASONS = {'spring': ('-04-01', '-06-01'),
           'summer': ('-06-01', '-09-01'),
           'fall':   ('-09-01', '-11-01')}
YEARS = [2021, 2022, 2023, 2024, 2025]      # the NFI-8 measurement period

# Cloud Score+ threshold. The monsoon leaves few clear scenes, and a stricter cut empties the summer
# composite over parts of the country.
CS_THRESHOLD = 0.60

# EVI is computed on the integer scale of S2_SR_HARMONIZED (reflectance x 10,000), so the soil term
# L = 1 of the reflectance-scale formula is 10,000 here.
EVI_L = 10000

BAND_NAMES = [f's2_{season}_{band}' for season in SEASONS for band in S2_BANDS]   # 42, season-major


def add_indices(img):
    """Add the four indices, computed once on the seasonal median, so the number of scenes a pixel
    kept after cloud masking cannot change what an index means between seasons."""
    return img.addBands([
        img.normalizedDifference(['B8', 'B4']).rename('NDVI'),
        img.normalizedDifference(['B8', 'B11']).rename('NDWI'),
        img.normalizedDifference(['B8', 'B5']).rename('NDRE'),
        img.expression('2.5*((N-R)/(N+6*R-7.5*B+L))',
                       {'N': img.select('B8'), 'R': img.select('B4'),
                        'B': img.select('B2'), 'L': EVI_L}).rename('EVI')])


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
