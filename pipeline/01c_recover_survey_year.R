# 01c_recover_survey_year.R — recover the year each plot was actually measured.
# ------------------------------------------------------------------------------------------------
# The extracted records carry the inventory cycle but not the year within it, and the cycle spans
# five years. Stage 07i reads the year to test whether pooling plots measured across the period
# against one image composite biases the result, so it is recovered here.
#
# The year column is auto-detected because its name differs between releases of the source data.
#
# In : data/raw/nfi_curated/NFI_Full.RData
# Out: data/interim/nfi8/plot_year.csv   (Plot_ID, Survey_year)
#
#   source(here::here("pipeline", "01c_recover_survey_year.R"))   # run in RStudio
# ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("R", "paths.R"))

CYCLE <- 8L
rdata <- require_file(here("data", "raw", "nfi_curated", "NFI_Full.RData"))
load(rdata); stopifnot(exists("tree")); setDT(tree)

# ---- 1. detect the survey-year / date column ----
cand <- c("Year","year","YEAR","Survey_year","SurveyYear","Measure_year","Inv_year","Inventory_year",
          "측정연도","조사년도","조사연도","조사년","조사일자","조사일","조사년월일","Date","date")
nm  <- names(tree)
hit <- intersect(cand, nm)
if (!length(hit)) hit <- nm[grepl("year|연도|일자|date|조사년|조사일", nm, ignore.case = TRUE)]
cat("candidate year/date columns:", if (length(hit)) paste(hit, collapse = ", ") else "(none)", "\n")
if (!length(hit)) { cat("\nnames(tree):\n"); print(nm)
  stop("No year or date column found; set `yc` by hand from the names printed above.") }
yc <- hit[1]; cat("using column:", yc, "  (sample:", paste(utils::head(tree[[yc]], 3), collapse = ", "), ")\n")

# ---- 2. derive 4-digit year (handles numeric year or date string) ----
to_year <- function(x) {
  if (is.numeric(x)) { y <- as.integer(x); y[y < 1990 | y > 2030] <- NA; return(y) }
  y <- suppressWarnings(as.integer(substr(gsub("[^0-9]", "", as.character(x)), 1, 4)))  # leading 4 digits
  y[y < 1990 | y > 2030] <- NA; y
}
tree[, Survey_year := to_year(get(yc))]
cat("\nSurvey_year by cycle (verify cycle 8 = 2021-2025):\n")
print(tree[, .(n = .N, yr_min = min(Survey_year, na.rm = TRUE),
               yr_max = max(Survey_year, na.rm = TRUE), na = sum(is.na(Survey_year))), by = Nth][order(Nth)])

# ---- 3. per-plot year (a plot's trees share one survey year → mode of non-NA) ----
mode1 <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_integer_)
                       as.integer(names(sort(table(x), decreasing = TRUE))[1]) }
py <- tree[Nth == CYCLE, .(Survey_year = mode1(Survey_year)), by = .(Plot_ID = Subplot_new)]
o  <- here("data", "interim", VERSION); dir.create(o, showWarnings = FALSE, recursive = TRUE)
fwrite(py, file.path(o, "plot_year.csv"))
cat(sprintf("%s: %d plots | year range %s | missing = %d\n", VERSION, nrow(py),
            paste(range(py$Survey_year, na.rm = TRUE), collapse = "-"), sum(is.na(py$Survey_year))))
