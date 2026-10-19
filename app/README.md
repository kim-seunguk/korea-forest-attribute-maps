# Viewer

This directory holds the browser viewer over the dataset's public Earth Engine image collection, so
the maps can be viewed without downloading them.

| | |
|---|---|
| viewer | <https://ee-korea-forest-maps.projects.earthengine.app/view/korea-forest-maps> |
| collection | `projects/ee-korea-forest-maps/assets/korea-forest-attributes-v2` |
| plot ids | `projects/ee-korea-forest-maps/assets/plot_ids_v2` |

| file | what it is |
|---|---|
| `layers.json` | the band each published raster carries, and the display range of each band |
| `viewer.js` | the Earth Engine App script |

`pipeline/08e_publish_earthengine.py` creates the collection.

`layers.json` holds the band of each published raster and the display range of each band. The
block at the top of `viewer.js` holds the same ranges and the species colors and names from
`data/species_colors.csv`.

## What is published

Fourteen images. Ten are the deposited rasters, unchanged. Four are the attributes of the plot-index
family, which `08e` derives from the deposited plot index and its lookup table by the same
substitution `05c` performs, so every layer here is reproducible from the archive alone.

Deriving from the deposited lookup puts both families on the published 37-class coding, so one
legend and one palette serve both.

## Display ranges

The ranges in `layers.json` are p5 to p95 of the values, **pooled across every layer that carries the
band**. A fraction or an index keeps its definitional 0 to 1, and a
class raster gets no range at all.

Because the range is pooled, a color means the same value on both sides of the wipe. The plot-index
family has the wider spread, so it reaches further into both ends of the scale.

`range_from` on each band records which rule produced it.

## Species colors

The species colors come from `data/species_colors.csv`, which also colors the map figures: a
color for each species named in their legends, and one shared color for the rest. The palette has
a row for each of the 37 class codes, 0 to 36.
Clicking a cell reports the species name. Per-class accuracy ships with the dataset in
`dominant_species_table.csv`, and the map is far more reliable for the common species than for the
rare ones.
