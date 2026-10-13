"""08e_publish_earthengine.py — publish the dataset as a public Earth Engine collection.

The archive is the citable copy. This is the one people can open: it ingests the maps as one public
ImageCollection, so a reader can look at them in a browser and use them in Earth Engine without
downloading several gigabytes.

Everything published here comes from the deposit, and only from the deposit. The ten deposited
rasters are staged as they are. The four attributes of the plot-index family are DERIVED here, by
substituting each cell's plot index with that plot's value from the deposited lookup table — the
same operation 05c performs, done from the archived inputs so that every layer in the viewer is
reproducible from the archive alone.

Deriving rather than uploading outputs/nfi8/attributes_ftm/ is not a convenience. Those rasters carry
a dominant-species legend of their own, built from whatever species the reference library happened
to contain (139 classes), so their codes mean different species from the delivered 37-class map. The
deposited lookup carries dominant_species_code in the published coding, so deriving from it puts both
families in one code space and lets one legend and one palette serve both.

The manifest is built from `data_records.csv`, which 08c writes by opening the rasters themselves.
Nothing about the grid, the data type or the nodata value is restated here; if the product changes,
this follows it.

Overview resampling matches what 06d bakes into the delivered files, so a zoomed-out view in Earth
Engine and a zoomed-out view in a GIS show the same thing:

    dominant species, plot index (discrete)   MODE
    everything continuous                     MEAN

MODE, not a sample, for the plot index: sampling keeps whichever single child cell it lands on and
leaves a hole wherever that one is empty, which thins the index map badly at low zoom. The mode of a
block is a plot that really was matched somewhere in it.

CONFIGURE  Both are yours, and neither is in this repository:
    EE_PROJECT    the Earth-Engine-enabled Cloud project the assets live under
    GCS_BUCKET    a Cloud Storage bucket in the same project, for staging

IN   outputs/nfi8/_deposit_v2/              the verified deposit, and its data_records.csv
OUT  outputs/nfi8/_ee_staging_v2/*.tif      the derived plot-index attribute rasters
     outputs/nfi8/_ee_manifests_v2/*.json   one ingestion manifest per image

Version 2 of the dataset goes into a collection of its own, korea-forest-attributes-v2, with its own
plot-id table, plot_ids_v2.

RUN  python pipeline/08e_publish_earthengine.py            # derive and write the manifests only
     python pipeline/08e_publish_earthengine.py --run      # create the bucket, stage, ingest
     python pipeline/08e_publish_earthengine.py --status   # how far ingestion has got
     python pipeline/08e_publish_earthengine.py --public   # let anyone read the collection

Everything runs through the Python client libraries, so neither the gcloud SDK nor the earthengine
command line has to be installed. `earthengine authenticate` once is enough: that credential carries
cloud-platform and devstorage.full_control, which covers the bucket too.

! Run 08a_verify_raster_integrity.py first. Publishing a corrupted raster is worse than not
  publishing it, because the damage is then downstream of a citable identifier.
! The deposit is read only. Derived rasters are written to a separate staging directory.
"""
import csv
import io
import json
import os
import sys

import numpy as np
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
DEPOSIT = os.path.join(ROOT, 'outputs', 'nfi8', '_deposit_v2')
STAGING = os.path.join(ROOT, 'outputs', 'nfi8', '_ee_staging_v2')
MANIFESTS = os.path.join(ROOT, 'outputs', 'nfi8', '_ee_manifests_v2')

COLLECTION = 'korea-forest-attributes-v2'
TIERS = ('ftm',)                    # the predictor set the deposit carries; version 2 file names carry no token
TABLE = 'plot_ids_v2'               # the plot-id table, read by the viewer

# The four attributes, as they are named in the delivered files and as they are carried in the
# deposited lookup table. dtype and nodata follow the delivered direct maps, so the two families
# are read the same way.
DERIVED = [
    # var name           lookup column               dtype      nodata  unit            range
    ('conifer_fraction', 'coniferous_fraction',      'float32', np.nan, '0-1 fraction', '0-1'),
    ('biomass',          'aboveground_biomass_Mgha', 'float32', np.nan, 'Mg ha-1',      '>=0'),
    ('shannon',          'shannon_diversity',        'float32', np.nan, 'unitless',     '>=0'),
    ('dominant_species', 'dominant_species_code',    'uint8',   255,    'class code',   '0-36'),
]
# One band id per quantity, shared by both families, so switching family in a viewer does not
# change what the band is called. Keyed on the variable text data_records.csv carries and on the
# names used above for the derived layers.
BAND = {
    'coniferous fraction': 'conifer_fraction',      'conifer_fraction': 'conifer_fraction',
    'dominant species': 'dominant_species',         'dominant_species': 'dominant_species',
    'aboveground biomass': 'biomass',               'biomass': 'biomass',
    'tree-species diversity (Shannon H′)': 'shannon',                'shannon': 'shannon',
    'coniferous fraction SD': 'conifer_fraction_sd',
    'biomass SD': 'biomass_sd',
    'Shannon SD': 'shannon_sd',
    'dominant-species uncertainty': 'dominant_species_uncertainty',
    'matched NFI plot index': 'plot_index',
    'random-forest distance to the matched reference plot': 'nearest_reference_distance',
}
# Discrete values, where averaging a block would invent a code that means something else. The
# uncertainty layers carry a class word in their name and are continuous, so they are excluded
# first — the same guard 06d applies before it treats anything as a class raster.
DISCRETE = ('dominant_species', 'imputation', 'plot_index')

PROJECT = os.environ.get('EE_PROJECT', '')
BUCKET = os.environ.get('GCS_BUCKET', '')
if not PROJECT or not BUCKET:
    raise SystemExit('set EE_PROJECT and GCS_BUCKET first (see the header)')

RUN = '--run' in sys.argv
BUCKET_LOCATION = os.environ.get('GCS_LOCATION', 'asia-northeast3')


def pyramiding(band):
    if band.endswith(('_sd', '_uncertainty')) or band == 'nearest_reference_distance':
        return 'MEAN'
    return 'MODE' if band in DISCRETE else 'MEAN'


def band_name(variable):
    try:
        return BAND[variable]
    except KeyError:
        raise SystemExit('no band id for %r - add it to BAND' % variable)


def read_lookup(tier):
    path = os.path.join(DEPOSIT, 'plot_lookup.csv')
    with io.open(path, encoding='utf-8-sig') as fh:
        rows = list(csv.DictReader(fh))
    idx = [int(r['plot_index']) for r in rows]
    if sorted(idx) != list(range(1, len(rows) + 1)):
        raise SystemExit('%s: plot_index is not 1..N in order' % os.path.basename(path))
    return rows


def derive(tier):
    """One raster per attribute, by substituting the plot index with the matched plot's value."""
    src = os.path.join(DEPOSIT, 'plot_index_30m.tif')
    rows = read_lookup(tier)
    n = len(rows)
    out_paths = []
    for var, col, dtype, nodata, _u, _r in DERIVED:
        dst = os.path.join(STAGING, 'imputed_%s_30m.tif' % var)
        out_paths.append((var, dst))
        if os.path.exists(dst):
            print('  %-42s already staged' % os.path.basename(dst))
            continue
        # index 0 is unused, so slot 0 carries the nodata value and an empty cell resolves to it.
        # A plot with an empty value in the lookup keeps nodata rather than becoming a zero.
        lut = np.full(n + 1, nodata, dtype=dtype)
        blank = 0
        for r in rows:
            v = str(r[col]).strip()
            if not v:
                blank += 1
                continue
            lut[int(r['plot_index'])] = np.dtype(dtype).type(float(v))
        with rasterio.open(src) as s:
            prof = s.profile.copy()
            prof.update(dtype=dtype, nodata=nodata, compress='LZW', tiled=True,
                        blockxsize=256, blockysize=256, BIGTIFF='IF_SAFER')
            prof.pop('photometric', None)
            with rasterio.open(dst, 'w', **prof) as d:
                for _, win in s.block_windows(1):
                    a = s.read(1, window=win)
                    valid = np.isfinite(a)
                    if s.nodata is not None and np.isfinite(s.nodata):
                        valid &= a != s.nodata                  # the index is uint16 with 0 as nodata
                    i = np.where(valid, np.rint(a), 0).astype('int64')
                    np.clip(i, 0, n, out=i)
                    d.write(lut[i].astype(dtype), 1, window=win)
        print('  %-42s derived from plot_index_30m.tif%s'
              % (os.path.basename(dst),
                 '  (%d plots have no value; left nodata)' % blank if blank else ''))
    return out_paths


def main():
    records = os.path.join(DEPOSIT, 'data_records.csv')
    if not os.path.exists(records):
        raise SystemExit('missing %s - run 08c and 08d first' % records)
    with io.open(records, encoding='utf-8-sig') as fh:
        rows = list(csv.DictReader(fh))

    os.makedirs(STAGING, exist_ok=True)
    os.makedirs(MANIFESTS, exist_ok=True)

    print('== deriving the plot-index attribute rasters ==')
    images = []                                  # (source path, stem, band, pyramiding, properties)
    for tier in TIERS:
        derived = dict(derive(tier))
        for var, _col, _dt, _nd, unit, vrange in DERIVED:
            dst = derived[var]
            images.append((dst, os.path.splitext(os.path.basename(dst))[0], band_name(var),
                           {'family': 'plot-index', 'variable': var,
                            'unit': unit, 'valid_range': vrange,
                            'source': 'derived from the deposited plot index and lookup table'}))

    print('\n== deposited rasters ==')
    for r in rows:
        src = os.path.join(DEPOSIT, r['file'])
        if not os.path.exists(src):
            raise SystemExit('missing %s - stage the deposit first' % src)
        images.append((src, os.path.splitext(r['file'])[0], band_name(r['variable']),
                       {'family': r['family'], 'variable': r['variable'], 'unit': r['unit'],
                        'valid_range': r['valid_range']}))

    coll = 'projects/%s/assets/%s' % (PROJECT, COLLECTION)
    plan = []
    for src, stem, band, props in images:
        blob = '%s/%s' % (COLLECTION, os.path.basename(src))
        with rasterio.open(src) as s:
            nodata = s.nodata
        manifest = {
            'name': '%s/%s' % (coll, stem),
            'tilesets': [{'id': stem, 'sources': [{'uris': ['gs://%s/%s' % (BUCKET, blob)]}]}],
            'bands': [{'id': band, 'tilesetId': stem, 'pyramidingPolicy': pyramiding(band)}],
            'properties': props,
        }
        if nodata is not None and np.isfinite(nodata):
            manifest['missingData'] = {'values': [float(nodata)]}
        with io.open(os.path.join(MANIFESTS, stem + '.json'), 'w', encoding='utf-8') as fh:
            fh.write(json.dumps(manifest, indent=2, ensure_ascii=False))
        plan.append((src, blob, manifest))
        print('  %-46s %-5s %s' % (os.path.basename(src), pyramiding(band), band))

    total = sum(os.path.getsize(s) for s, _, _ in plan)
    print('')
    print('%d images, %.1f GB. manifests -> %s' % (len(plan), total / 1e9, MANIFESTS))
    print('collection  %s' % coll)
    print('staging     gs://%s/%s/' % (BUCKET, COLLECTION))
    if not RUN:
        print('')
        print('Nothing has been uploaded. --run creates the bucket, stages and ingests.')
        return
    creds = credentials()
    upload(plan, creds)
    ingest(coll, plan, creds)
    ingest_tables(creds)


def credentials():
    """Application-default credentials when they are configured, otherwise the ones the Earth
    Engine command line already stored. That credential carries cloud-platform and
    devstorage.full_control, so it covers the bucket as well as the ingestion, and a reader who has
    run `earthengine authenticate` needs nothing further."""
    import google.auth
    try:
        creds, _ = google.auth.default(scopes=['https://www.googleapis.com/auth/cloud-platform'])
        return creds
    except Exception:
        pass
    import ee
    from google.oauth2.credentials import Credentials
    path = os.path.join(os.path.expanduser('~'), '.config', 'earthengine', 'credentials')
    if not os.path.exists(path):
        raise SystemExit('no credentials: run `earthengine authenticate` first')
    d = json.load(io.open(path, encoding='utf-8'))
    return Credentials(None, refresh_token=d['refresh_token'],
                       token_uri='https://oauth2.googleapis.com/token',
                       client_id=ee.oauth.CLIENT_ID, client_secret=ee.oauth.CLIENT_SECRET,
                       scopes=d.get('scopes'))


def upload(plan, creds):
    """Stage every image in Cloud Storage. A blob already there at the same size is left alone, so
    an interrupted run continues instead of starting over."""
    from google.cloud import storage
    client = storage.Client(project=PROJECT, credentials=creds)
    bucket = client.bucket(BUCKET)
    if not bucket.exists():
        bucket = client.create_bucket(bucket, location=BUCKET_LOCATION)
        print('created gs://%s in %s' % (BUCKET, BUCKET_LOCATION))
    print('')
    print('== staging to Cloud Storage ==')
    for i, (src, blob_name, _) in enumerate(plan, 1):
        size = os.path.getsize(src)
        blob = bucket.blob(blob_name)
        if blob.exists() and blob.size == size:
            print('  [%2d/%d] %-46s already staged' % (i, len(plan), os.path.basename(src)))
            continue
        blob.chunk_size = 32 * 1024 * 1024
        blob.upload_from_filename(src)
        print('  [%2d/%d] %-46s %6.0f MB' % (i, len(plan), os.path.basename(src), size / 1e6))


def ingest(coll, plan, creds):
    """Create the collection, then start one ingestion task per image."""
    import ee
    ee.Initialize(credentials=creds, project=PROJECT)
    try:
        ee.data.createAsset({'type': 'ImageCollection'}, coll)
        print('created collection %s' % coll)
    except Exception as e:
        if 'exist' not in str(e).lower():
            raise
        print('collection already exists: %s' % coll)
    print('')
    print('== starting ingestion ==')
    for i, (_src, _blob, manifest) in enumerate(plan, 1):
        name = manifest['name'].rsplit('/', 1)[-1]
        task = ee.data.newTaskId()[0]
        ee.data.startIngestion(task, manifest, allow_overwrite=True)
        print('  [%2d/%d] %-46s task %s' % (i, len(plan), name, task))
    print('')
    print('Ingestion is asynchronous. Poll it with --status.')


def ingest_tables(creds):
    """A table mapping the plot index to the NFI plot id, for the viewer to look up.

    Earth Engine ingests a CSV as a table only with an x and a y column, so every row carries 0, 0.
    Those are placeholders: no plot location enters the table."""
    import ee
    from google.cloud import storage
    ee.Initialize(credentials=creds, project=PROJECT)
    bucket = storage.Client(project=PROJECT, credentials=creds).bucket(BUCKET)
    print('')
    print('== plot-id tables ==')
    for tier in TIERS:
        rows = read_lookup(tier)
        name = TABLE
        path = os.path.join(STAGING, name + '.csv')
        with io.open(path, 'w', encoding='utf-8', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(['plot_index', 'nfi_plot_id', 'lon', 'lat'])
            for r in rows:
                w.writerow([r['plot_index'], r['nfi_plot_id'], 0, 0])
        blob = '%s/%s.csv' % (COLLECTION, name)
        bucket.blob(blob).upload_from_filename(path)
        task = ee.data.newTaskId()[0]
        ee.data.startTableIngestion(task, {
            'name': 'projects/%s/assets/%s' % (PROJECT, name),
            'sources': [{'uris': ['gs://%s/%s' % (BUCKET, blob)],
                         'xColumn': 'lon', 'yColumn': 'lat'}]}, allow_overwrite=True)
        print('  %-16s %6d rows  task %s' % (name, len(rows), task))


def status():
    """How far ingestion has got, and whether the collection is readable by anyone."""
    import ee
    ee.Initialize(credentials=credentials(), project=PROJECT)
    coll = 'projects/%s/assets/%s' % (PROJECT, COLLECTION)
    try:
        assets = ee.data.listAssets({'parent': coll}).get('assets', [])
    except Exception as e:
        raise SystemExit('cannot list %s: %s' % (coll, e))
    running = [op for op in ee.data.listOperations() if not op.get('done')]
    acl = ee.data.getAssetAcl(coll)
    print(coll)
    print('  %d image(s) ingested' % len(assets))
    print('  %d task(s) still running' % len(running))
    print('  readable by anyone: %s' % bool(acl.get('all_users_can_read')))
    return assets, running


def make_public():
    """Share the collection with everyone, which is what a public app reads through.

    Only the collection takes an ACL. An image inside a collection cannot hold one of its own --
    asking for its policy is rejected outright -- because it inherits the collection's. So one call
    covers every image, and the count below is what proves it did.
    """
    import ee
    ee.Initialize(credentials=credentials(), project=PROJECT)
    coll = 'projects/%s/assets/%s' % (PROJECT, COLLECTION)
    for asset in [coll] + ['projects/%s/assets/%s' % (PROJECT, TABLE % t) for t in TIERS]:
        try:
            acl = ee.data.getAssetAcl(asset)
        except Exception as e:
            print('  skip    %s (%s)' % (asset, str(e)[:50]))
            continue
        acl['all_users_can_read'] = True
        ee.data.setAssetAcl(asset, json.dumps(acl))
        print('  public  %s' % asset)
    n = len(ee.data.listAssets({'parent': coll}).get('assets', []))
    print('  %d image(s) inherit the collection' % n)


if __name__ == '__main__':
    if '--status' in sys.argv:
        status()
    elif '--public' in sys.argv:
        make_public()
    else:
        main()
