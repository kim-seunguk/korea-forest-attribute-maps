# R/run_utils.R — shared run instrumentation (timing · progress · adaptive sizing).
#
# Instrumentation for the heavy runs. Every long stage should print its per-stage wall time and a
# total breakdown, show brief progress, and let a slow pilot shrink the next run automatically.
# This file is the one place those helpers live; source it AFTER R/paths.R. Base R only, plus here
# for the timings path.
#
# Typical use:
#   source(here::here("R", "paths.R")); source(here::here("R", "run_utils.R"))
#   cfg <- run_sizing("nfi8", "pilot_ladder", pilot = PILOT, budget_min = 30)   # auto-shrinks if last pilot was slow
#   run_banner("pilot_ladder", PILOT, cfg, n = n)
#   tm <- run_timer("pilot_ladder", version = "nfi8", budget_min = 30)
#   tm$stage("load");  ... ; tm$stage("folds"); ... ; tm$stage("cv")
#   for (i in seq_len(K)) { progress(i, K, sprintf("fold %d", i)); ... }
#   tm$done(n = n, params = sprintf("K=%d ntree=%d nboot=%d", cfg$K, cfg$ntree, cfg$nboot))

# ---- 1. progress + timestamped log --------------------------------------------------------
# log_stage(...): a timestamped line, for stages long enough that a bare message says too little.
log_stage <- function(...) cat(sprintf("[%s]", format(Sys.time(), "%H:%M:%S")), ..., "\n")
# progress(i, n, label): one-line "[i/n] label" for fold / target / predictor-set loops.
progress  <- function(i, n, label = "") cat(sprintf("[%d/%d] %s\n", i, n, label))

# ---- 2. run timer (per-stage + cumulative + total breakdown) ------------------------------
# Returns an object with $stage(name) and $done(n, params). Each $stage() closes the previous
# stage (printing its minutes + cumulative); $done() prints the total + per-stage breakdown and,
# if version given, appends to outputs/<version>/run_timings.csv. budget_min > 0 → over-budget warn.
run_timer <- function(run = "run", version = NULL, budget_min = Inf, log = TRUE) {
  e <- new.env(parent = emptyenv())
  e$run <- run; e$version <- version; e$budget_min <- budget_min; e$log <- log
  e$t_run0 <- Sys.time(); e$t_stage0 <- e$t_run0; e$cur <- NA_character_; e$stages <- list()
  flush_stage <- function() {
    if (!is.na(e$cur)) {
      sec <- as.numeric(difftime(Sys.time(), e$t_stage0, units = "secs"))
      prev <- if (is.null(e$stages[[e$cur]])) 0 else e$stages[[e$cur]]
      e$stages[[e$cur]] <- prev + sec
      if (e$log) cat(sprintf("  [%-16s] %6.1f min  (cum %.1f min)\n", e$cur, sec / 60,
                             as.numeric(difftime(Sys.time(), e$t_run0, units = "mins"))))
    }
  }
  e$stage <- function(name) {
    flush_stage(); e$cur <- name; e$t_stage0 <- Sys.time()
    if (e$log) log_stage(sprintf(">> %s", name)); invisible()
  }
  e$done <- function(n = NA, params = "") {
    flush_stage()
    total <- as.numeric(difftime(Sys.time(), e$t_run0, units = "mins"))
    brk <- if (length(e$stages)) paste(sprintf("%s=%.1f", names(e$stages),
                                               unlist(e$stages) / 60), collapse = "  ") else "-"
    cat(sprintf("===== [%s] total %.1f min | %s =====\n", e$run, total, brk))
    if (!is.null(e$version)) record_timing(e$version, e$run, e$stages, n, params)
    if (is.finite(e$budget_min) && total > e$budget_min)
      cat(sprintf("WARN over budget: %.1f min > %.0f min — design the NEXT run smaller.\n",
                  total, e$budget_min))
    invisible(total)
  }
  e
}

# ---- 3. persist timings → outputs/<version>/run_timings.csv (append) ----------------------
record_timing <- function(version, run, stages, n = NA, params = "") {
  outd <- here::here("outputs", version); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(outd, "run_timings.csv"); ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  rows <- if (length(stages))
    data.frame(timestamp = ts, run = run, stage = names(stages),
               sec = round(unlist(stages), 1), n = n, params = params,
               row.names = NULL, stringsAsFactors = FALSE)
  else data.frame(timestamp = character(0), run = character(0), stage = character(0),
                  sec = numeric(0), n = integer(0), params = character(0))
  rows <- rbind(rows, data.frame(timestamp = ts, run = run, stage = "TOTAL",
                                 sec = round(sum(unlist(stages)), 1), n = n, params = params))
  utils::write.table(rows, f, sep = ",", append = file.exists(f),
                     col.names = !file.exists(f), row.names = FALSE, qmethod = "double")
  invisible(f)
}

# ---- 4. adaptive run sizing ---------------------------------------------------------------
# estimate_total(): project total wall-time from one measured unit (e.g. fold 1), print it, and
# flag if it would blow the budget — call it after the first fold to abort/shrink early.
estimate_total <- function(unit_sec, n_units, budget_min = Inf, label = "run") {
  est_min <- unit_sec * n_units / 60
  cat(sprintf("[est] %s: unit %.1fs x %d = ~%.1f min%s\n", label, unit_sec, n_units, est_min,
              if (is.finite(budget_min)) sprintf(" (budget %.0f min)", budget_min) else ""))
  over <- is.finite(budget_min) && est_min > budget_min
  if (over) cat(sprintf("WARN projected over budget — shrink ~x%.2f.\n", budget_min / est_min))
  list(est_min = est_min, over = over, shrink = if (over) budget_min / est_min else 1)
}

# run_sizing(): pick K, ntree and nboot for this run. If this is a pilot and the last recorded pilot
# for the same (version, run) ran over budget_min, the pilot defaults shrink, floored at a third, so
# the next pilot is cheaper. run_timings.csv is what carries that memory between runs.
#
# K stays at the frozen fold count in both settings. The fold loop runs 1..K over folds that were
# frozen at 10, so a smaller K would quietly evaluate a fraction of the plots and report it as if it
# were the whole sample. A pilot buys its speed from ntree and nboot alone.
run_sizing <- function(version, run, pilot = TRUE, budget_min = Inf,
                       pilot_defaults = list(K = 10L, ntree = 200L, nboot = 300L),
                       full_defaults  = list(K = 10L, ntree = 500L, nboot = 1000L)) {
  cfg <- if (pilot) pilot_defaults else full_defaults
  if (pilot && is.finite(budget_min)) {
    f <- here::here("outputs", version, "run_timings.csv")
    if (file.exists(f)) {
      d <- utils::read.csv(f, stringsAsFactors = FALSE)
      d <- d[d$run == run & d$stage == "TOTAL", , drop = FALSE]
      if (nrow(d)) {
        last <- d$sec[nrow(d)] / 60
        if (last > budget_min) {
          s <- max(1 / 3, budget_min / last)
          cfg$ntree <- max(100L, as.integer(round(cfg$ntree * s)))   # trees and resamples only; the folds are frozen
          cfg$nboot <- max(100L, as.integer(round(cfg$nboot * s)))
          cat(sprintf("[adapt] last pilot %.1f min > budget %.0f min -> shrink x%.2f: ntree=%d nboot=%d\n",
                      last, budget_min, s, cfg$ntree, cfg$nboot))
        }
      }
    }
  }
  cfg
}

# run_banner(): one-line echo of the active config so it's obvious what setting is running.
run_banner <- function(run, pilot, cfg, n = NA, extra = "") {
  cat(sprintf("\n===== %s [%s]  K=%s ntree=%s nboot=%s  n=%s %s =====\n", run,
              if (isTRUE(pilot)) "PILOT" else "FULL", cfg$K, cfg$ntree, cfg$nboot, n, extra))
}

# register_cluster(): a PSOCK cluster whose workers inherit the master's library paths. Each worker
# is a fresh R session and starts with the default .libPaths(), so packages installed in a user
# library or under renv are invisible to it and foreach fails to load. Passing the master's paths in
# is what makes the parallel stages work; every one of them registers its cluster through here.
register_cluster <- function(n_max = 16L) {
  cl <- parallel::makeCluster(min(n_max, max(1L, parallel::detectCores() - 1L)))
  parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())
  doParallel::registerDoParallel(cl)
  cl
}
