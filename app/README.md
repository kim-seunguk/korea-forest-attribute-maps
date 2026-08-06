# Viewer

The dataset published as a public Earth Engine image collection, with a browser viewer over it, so
that the maps can be looked at without downloading them.

| | |
|---|---|
| viewer | <https://ee-korea-forest-maps.projects.earthengine.app/view/korea-forest-maps> |
| collection | `projects/ee-korea-forest-maps/assets/korea-forest-attributes` |
| plot ids | `projects/ee-korea-forest-maps/assets/plot_ids_{base,ftm}` |

| file | what it is |
|---|---|
| `layers.json` | the band each published raster carries, and the display range of each band |
| `viewer.js` | the Earth Engine App script |

`pipeline/08e_publish_earthengine.py` is what creates the collection.

## What is published

Twenty of the images are the deposited rasters, unchanged. Eight are the four attributes of the
plot-index family, derived by `08e` from the deposited plot index and its lookup table — the same
substitution `05c` performs, done from the archived inputs so that every layer here is reproducible
from the archive alone.

They are derived rather than taken from `attributes_*/` for a reason. Those rasters carry a
dominant-species legend built per tier from whichever species the reference library happened to
contain, so the same code means a different species in each tier and a different species again in
the delivered maps. The deposited lookup carries the published 39-class coding, so deriving from it
puts both families in one code space, and one legend and one palette then serve both.

## Display ranges

The published figures use a histogram-equalised stretch, which has no fixed minimum and maximum to
copy, so the ranges in `layers.json` are measured instead: p5 to p95 of the values, **pooled across
every layer that carries the band**. A fraction or an index keeps its definitional 0 to 1, and a
class raster gets no range at all.

Pooling is what makes the side-by-side comparison mean anything. With a stretch per layer, the same
colour would stand for a different value on each side of the wipe. It also leaves the real
difference visible: the plot-index family holds a wider spread than the direct family, and on a
shared scale that shows as reaching further into both ends rather than being normalised away.

`range_from` on each band records which rule produced it.

## Species colours

From `data/species_colors.csv`, the same table the figures are drawn from: a colour for each species
named in the map figures' legends, and one shared colour for the rest. Clicking a cell reports the
full 39-class name from `dominant_species_table.csv`, which also carries the per-class accuracy —
the map is far more reliable for the common species than for the rare ones.
