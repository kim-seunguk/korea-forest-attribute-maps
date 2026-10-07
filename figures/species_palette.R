# species_palette.R — the display symbology for the dominant-species maps.
#
# Reads data/species_colors.csv, which gives a colour to each species named in the map figures'
# legends. The rest collapse into a single "Others" class, drawn in the CSV's last colour.
#
# This is figure symbology, and the same table is embedded in the delivered dominant-species raster
# as its internal colour table (06d), so a GIS colours that map on open with these colours.
#
#   source(here::here("figures", "species_palette.R"))
#
# Exposes CODE2LABEL (code to display label), DISP_COLORS (label to hex), DISP_LEVELS (legend order)
# and DOMSP_SCI (Korean species name to scientific name).
suppressPackageStartupMessages({ library(here); library(data.table) })
.sym <- fread(here("data","species_colors.csv"), encoding="UTF-8")
CODE2LABEL  <- setNames(.sym$label, as.character(.sym$code))          # map code -> "sci_name" or "Others"
.disp       <- rbind(unique(.sym[major == TRUE,  .(label, hex)]),     # majors in the CSV's genus order ...
                     unique(.sym[major == FALSE, .(label, hex)]))     # ... then the single "Other N species" row
DISP_LEVELS <- .disp$label
DISP_COLORS <- setNames(.disp$hex, .disp$label)

.leg      <- fread(here("data","species_codes.csv"), encoding="UTF-8")
DOMSP_SCI <- setNames(.leg$sci_name, .leg$species)                    # back-compat helper
