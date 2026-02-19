install_witchtools <- function() {
  if (!"remotes" %in% rownames(installed.packages())) {
    install.packages("remotes", repos = "http://cloud.r-project.org")
  }
  remotes::install_github("witch-team/witchtools")
}

if (!"witchtools" %in% rownames(installed.packages())) {
  install_witchtools()
  if (!requireNamespace("witchtools")) stop("Package witchtools not found")
}
if (packageVersion("witchtools") < "0.4.0") {
  cat("Need a more recent version of witchtools. Trying update.\n")
  install_witchtools()
  cat(paste("Installed version:", packageVersion("witchtools"), "\n"))
  if (packageVersion("witchtools") < "0.4.0") {
    stop("Please install witchtools version >= 0.4.0.\n")
  }
}

library(witchtools)

pkgs <- c("data.table", "stringr", "docopt", "ggplot2")
invisible(lapply(pkgs, require_package))
require_gdxtools()

# ------------------------------------------------------------------
# Script summary (synthetic):
# 1) Load matching pulse/base and srmpulse/srm simulations.
# 2) In chunks, build compact impulse effects (Dtemp, Dforc, SRMminus)
#    on relative time t_rel=t-pulse_time and save to disk.
# 3) Reload compact tables and compute year-by-year quantiles by gas.
# 4) Compute coherent percentile trajectories by selecting real runs
#    that minimize RMSE to the percentile curve.
#    This includes Dtemp, Dforc, and SRMminus.
# 5) Overlay -SRM (from srmpulsemasked) on forcing panel and plot:
#    a) coherent-trajectory ribbons/median
#    b) year-by-year ribbons/median
# ------------------------------------------------------------------

extract_names_dt <- function(files_vec) {
  DT <- data.table(gdx = files_vec)
  DT[, file := str_remove_all(gdx, ".*(?=RCP)|.gdx")]
  DT[, gas := str_extract(file, "(?<=GAS).+?(?=_)")]
  DT[, rcp := str_extract(file, "(?<=RCP).+?(?=_)")]
  DT[, cool_rate := str_extract(str_remove(file, "RCP"), "(?<=RC).+?(?=_)")]
  DT[, pulse_time := str_extract(gdx, "(?<=PT).+?(?=_)")]
  DT[, geo_end := str_extract(str_remove(file, "ECS"), "(?<=EC).+?(?=_)")]
  DT[, geo_start := str_extract(file, "(?<=BC).+?(?=_)")]
  DT[, ecs := str_extract(file, "(?<=ECS).+?(?=_)")]
  DT[, tcr := str_extract(file, "(?<=TCR).+?(?=_)")]
  DT[, term := str_extract(file, "(?<=TER).+?(?=_)")]
  DT[, experiment := str_extract(file, "(?<=EXP).+?(?=_)")]
  DT
}

sanitize_dt <- function(DT) {
  setDT(DT)
  if ("t" %in% names(DT)) DT[, t := as.numeric(t)]
  if ("gdx" %in% names(DT)) DT <- DT[, !"gdx", with = FALSE]
  DT <- DT[t <= 480]
  unique(DT)
}

replace_num_na0 <- function(DT) {
  numcols <- names(DT)[vapply(DT, is.numeric, logical(1))]
  for (cname in numcols) set(DT, which(is.na(DT[[cname]])), cname, 0)
  invisible(NULL)
}

quantile_summary <- function(DT, value_col, group_cols) {
  out <- DT[, .(
    n = .N,
    median = median(get(value_col), na.rm = TRUE),
    p05 = as.numeric(quantile(get(value_col), 0.05, na.rm = TRUE, type = 8)),
    p25 = as.numeric(quantile(get(value_col), 0.25, na.rm = TRUE, type = 8)),
    p75 = as.numeric(quantile(get(value_col), 0.75, na.rm = TRUE, type = 8)),
    p95 = as.numeric(quantile(get(value_col), 0.95, na.rm = TRUE, type = 8))
  ), by = group_cols]
  out[, variable := value_col]
  setcolorder(out, c(group_cols, "variable", "n", "median", "p05", "p25", "p75", "p95"))
  out
}

select_percentile_runs <- function(DT, value_col, group_cols, run_cols, probs) {
  # Select coherent trajectories:
  # for each percentile p and group (e.g., gas), pick the simulation run
  # whose full path is closest to the p-quantile curve across t_rel.
  key_cols <- c(group_cols, run_cols)
  work <- DT[, c(key_cols, "t", "t_rel", value_col), with = FALSE]

  best_list <- list()
  traj_list <- list()

  for (p in probs) {
    # 1) Build pointwise percentile target curve q_p(t_rel)
    qcurve <- work[, .(
      q_value = as.numeric(quantile(get(value_col), p, na.rm = TRUE, type = 8))
    ), by = c(group_cols, "t_rel")]

    # 2) Compute squared deviation from target at each time
    dev <- merge(work, qcurve, by = c(group_cols, "t_rel"), all = FALSE)
    dev[, sqe := (get(value_col) - q_value)^2]

    # 3) Collapse to one distance per run: RMSE over time
    rmse <- dev[, .(rmse = sqrt(mean(sqe, na.rm = TRUE))), by = key_cols]
    setorder(rmse, rmse)
    # 4) Best run = minimum RMSE in each group
    best <- rmse[, .SD[1], by = group_cols]
    best[, `:=`(variable = value_col, percentile = p)]

    # 5) Return full trajectory of selected run
    traj <- merge(
      work[, .(t, t_rel, value = get(value_col)), by = key_cols],
      best[, c(key_cols, "variable", "percentile", "rmse"), with = FALSE],
      by = key_cols,
      all = FALSE
    )

    best_list[[length(best_list) + 1]] <- best
    traj_list[[length(traj_list) + 1]] <- traj
  }

  list(
    best = rbindlist(best_list, use.names = TRUE, fill = TRUE),
    traj = rbindlist(traj_list, use.names = TRUE, fill = TRUE)
  )
}

build_effects_pair <- function(TATM, tot_forcing, pulse_exp, base_exp, pair_label, temp_keys, forc_keys) {
  temp_p <- TATM[experiment == pulse_exp]
  temp_b <- TATM[experiment == base_exp]
  forc_p <- tot_forcing[experiment == pulse_exp]
  forc_b <- tot_forcing[experiment == base_exp]

  if (nrow(temp_p) == 0 || nrow(temp_b) == 0 || nrow(forc_p) == 0 || nrow(forc_b) == 0) {
    return(data.table())
  }

  pcols_t <- setdiff(names(temp_p), c("file", "experiment"))
  bcols_t <- setdiff(names(temp_b), c("file", "experiment"))
  keys_t <- intersect(temp_keys, intersect(setdiff(pcols_t, "value"), setdiff(bcols_t, "value")))

  pcols_f <- setdiff(names(forc_p), c("file", "experiment"))
  bcols_f <- setdiff(names(forc_b), c("file", "experiment"))
  keys_f <- intersect(forc_keys, intersect(setdiff(pcols_f, "value"), setdiff(bcols_f, "value")))

  if (length(keys_t) == 0 || length(keys_f) == 0) {
    stop(paste0("No valid merge keys for pair ", pair_label))
  }

  restore_from_pulse <- function(DT, cols) {
    for (cc in cols) {
      cp <- paste0(cc, "_pulse")
      if (!(cc %in% names(DT)) && cp %in% names(DT)) {
        DT[, (cc) := get(cp)]
      } else if (cc %in% names(DT) && cp %in% names(DT)) {
        DT[is.na(get(cc)), (cc) := get(cp)]
      }
    }
    DT
  }

  temp_eff <- merge(
    temp_p[, ..pcols_t],
    temp_b[, ..bcols_t],
    by = keys_t,
    all = FALSE,
    suffixes = c("_pulse", "_base")
  )
  setnames(temp_eff, c("value_pulse", "value_base"), c("temp_pulse", "temp_base"), skip_absent = TRUE)
  temp_eff <- restore_from_pulse(temp_eff, c("gas", "pulse_time", "cool_rate", "geo_start", "geo_end", "term"))

  forc_eff <- merge(
    forc_p[, ..pcols_f],
    forc_b[, ..bcols_f],
    by = keys_f,
    all = FALSE,
    suffixes = c("_pulse", "_base")
  )
  setnames(forc_eff, c("value_pulse", "value_base"), c("forc_pulse", "forc_base"), skip_absent = TRUE)
  forc_eff <- restore_from_pulse(forc_eff, c("gas", "pulse_time", "cool_rate", "geo_start", "geo_end", "term"))

  merge_keys <- intersect(
    setdiff(names(temp_eff), c("temp_pulse", "temp_base")),
    setdiff(names(forc_eff), c("forc_pulse", "forc_base"))
  )

  out <- merge(temp_eff, forc_eff, by = merge_keys, all = FALSE)
  out <- restore_from_pulse(out, c("gas", "pulse_time", "cool_rate", "geo_start", "geo_end", "term"))
  if (!"pulse_time" %in% names(out)) {
    pulse_time_cols <- grep("^pulse_time", names(out), value = TRUE)
    if (length(pulse_time_cols) > 0) out[, pulse_time := get(pulse_time_cols[1])]
  }
  if (!"gas" %in% names(out)) {
    gas_cols <- grep("^gas", names(out), value = TRUE)
    if (length(gas_cols) > 0) out[, gas := get(gas_cols[1])]
  }
  replace_num_na0(out)

  if (!"pulse_time" %in% names(out)) {
    stop(paste0("Missing pulse_time after merging for pair ", pair_label))
  }

  out <- out[t >= as.numeric(pulse_time)]
  out[, t_rel := as.numeric(t) - as.numeric(pulse_time)]
  out[, Dtemp := temp_pulse - temp_base]
  out[, Dforc := forc_pulse - forc_base]
  out[, pair := pair_label]
  out
}

'Analyze pulse effects from Monte Carlo outputs and compute time-aligned statistics for two pairings: pulse-base and srmpulse-srm.

Usage:
  Analyze_montecarlo_pulse_effects.R [-o <results>] [--res <output_folder>] [--hpc <run_hpc>] [--traj_q <traj_q>] [--chunk <chunk>] [--plot <plot_results>]

Options:
--res <results>        Results folder(s) with .gdx outputs. Separate multiple folders with -
-o <output_folder>     Output folder 
--hpc <run_hpc>        T/F if running on HPC 
--traj_q <traj_q>      Comma-separated percentile(s) for coherent-run selection 
--chunk <chunk>        Number of matched scenario-groups processed per batch ì
--plot <plot_results>  T/F to generate plots 
' -> doc

opts <- docopt(doc, version = "Analyze_montecarlo_pulse_effects")

res <- ifelse(is.null(opts[["res"]]), "Results_montecarlo", as.character(opts[["res"]]))
res <- str_split(res, "-")[[1]]
output_folder <- ifelse(is.null(opts[["o"]]), "Montecarlo", as.character(opts[["o"]]))
run_hpc <- ifelse(is.null(opts[["hpc"]]), F, as.logical(opts[["hpc"]]))
traj_q <- ifelse(is.null(opts[["traj_q"]]), "0.05,0.25,0.5,0.75,0.95", as.character(opts[["traj_q"]]))
traj_probs <- as.numeric(trimws(str_split(traj_q, ",")[[1]]))
chunk_n <- ifelse(is.null(opts[["chunk"]]), 400L, as.integer(opts[["chunk"]]))
plot_results <- ifelse(is.null(opts[["plot"]]), T, as.logical(opts[["plot"]]))
if (any(is.na(traj_probs)) || any(traj_probs <= 0 | traj_probs >= 1)) {
  stop("--traj_q must contain numbers strictly between 0 and 1, e.g. 0.75 or 0.25,0.5,0.75")
}
if (is.na(chunk_n) || chunk_n <= 0) stop("--chunk must be a positive integer")
if (is.na(plot_results)) stop("--plot must be T or F")

if (!dir.exists(output_folder)) dir.create(output_folder)
if (any(!dir.exists(res))) stop("Some result folders do not exist")

if (run_hpc == FALSE) {
  igdx()
} else {
  igdx("/work/cmcc/pa12520/gams40.4_linux_x64_64_sfx")
  logfile <- file(file.path(output_folder, "r_console_pulse_effects.log"), open = "wt")
  sink(logfile)
  sink(logfile, type = "message")
}

filelist <- unlist(lapply(res, function(folder) list.files(folder, pattern = ".gdx$", full.names = TRUE)))
filelist <- filelist[str_detect(filelist, "EXP")]
if (length(filelist) == 0) stop("No .gdx files found in the selected results folders")

cat("Sanitizing filenames...\n")
sanitized_names <- extract_names_dt(filelist)
scenario_names <- setdiff(names(sanitized_names), c("gdx", "file", "experiment"))

k_pulse <- setdiff(scenario_names, c("term", "cool_rate", "geo_start", "geo_end"))
k_base <- setdiff(scenario_names, c("term", "cool_rate", "geo_start", "geo_end", "gas", "pulse_time"))
k_no_term <- setdiff(scenario_names, "term")

all_pulse <- unique(sanitized_names[experiment == "pulse", ..k_pulse])
all_srmpulse <- unique(sanitized_names[experiment == "srmpulse", ..k_no_term])
if (nrow(all_pulse) == 0 && nrow(all_srmpulse) == 0) stop("No pulse or srmpulse scenarios found")

dt_pulse <- sanitized_names[experiment == "pulse"]
dt_base <- sanitized_names[experiment == "base"]
dt_srmpulse <- sanitized_names[experiment == "srmpulse"]
dt_srm <- sanitized_names[experiment == "srm"]
dt_srmpmasked <- sanitized_names[experiment == "srmpulsemasked"]

setkeyv(dt_pulse, k_pulse)
setkeyv(dt_base, k_base)
setkeyv(dt_srmpulse, k_no_term)
setkeyv(dt_srm, k_no_term)
setkeyv(dt_srmpmasked, k_no_term)

# Build matched scenario-groups once, then process in chunks to limit RAM.
job_list <- list()

for (i in seq_len(nrow(all_pulse))) {
  exp_row <- all_pulse[i]
  f_pulse <- dt_pulse[exp_row, on = k_pulse, mult = "first"]$gdx
  exp_base <- exp_row[, ..k_base]
  f_base <- dt_base[exp_base, on = k_base, mult = "first"]$gdx
  if (!is.na(f_pulse) && !is.na(f_base)) {
    job_list[[length(job_list) + 1]] <- data.table(
      pair = "pulse-base",
      f_pulse = f_pulse,
      f_base = f_base,
      f_srm_masked = NA_character_
    )
  }
}

for (i in seq_len(nrow(all_srmpulse))) {
  exp_row <- all_srmpulse[i]
  f_srmpulse <- dt_srmpulse[exp_row, on = k_no_term, mult = "first"]$gdx
  f_srm <- dt_srm[exp_row, on = k_no_term, mult = "first"]$gdx
  f_srmpmasked <- dt_srmpmasked[exp_row, on = k_no_term, mult = "first"]$gdx
  if (!is.na(f_srmpulse) && !is.na(f_srm)) {
    job_list[[length(job_list) + 1]] <- data.table(
      pair = "srmpulse-srm",
      f_pulse = f_srmpulse,
      f_base = f_srm,
      f_srm_masked = f_srmpmasked
    )
  }
}

if (length(job_list) == 0) stop("No matched file pairs found for pulse-base or srmpulse-srm")
pair_jobs <- rbindlist(job_list, use.names = TRUE, fill = TRUE)

effects_file <- file.path(output_folder, "pulse_effects_timeseries.csv")
srm_timeseries_file <- file.path(output_folder, "pulse_effects_srmminus_timeseries.csv")
output_files <- c(
  effects_file,
  srm_timeseries_file,
  file.path(output_folder, "pulse_effects_stats_by_gas.csv"),
  file.path(output_folder, "pulse_effects_percentile_runs.csv"),
  file.path(output_folder, "pulse_effects_percentile_trajectories.csv"),
  file.path(output_folder, "pulse_effects_srmminus_percentile_runs.csv"),
  file.path(output_folder, "pulse_effects_srmminus_percentile_trajectories.csv"),
  file.path(output_folder, "pulse_effects_srmminus_p50_trajectories.csv")
)
for (ff in output_files[file.exists(output_files)]) file.remove(ff)

effects_written <- FALSE
srm_written <- FALSE
n_chunks <- ceiling(nrow(pair_jobs) / chunk_n)

for (chunk_id in seq_len(n_chunks)) {
  i_start <- (chunk_id - 1L) * chunk_n + 1L
  i_end <- min(chunk_id * chunk_n, nrow(pair_jobs))
  jobs_chunk <- pair_jobs[i_start:i_end]

  files_chunk <- unique(c(jobs_chunk$f_pulse, jobs_chunk$f_base, jobs_chunk$f_srm_masked))
  files_chunk <- files_chunk[!is.na(files_chunk)]
  if (length(files_chunk) == 0) next

  cat("Processing chunk", chunk_id, "of", n_chunks, "with", length(files_chunk), "files...\n")

  sanitized_chunk <- sanitized_names[gdx %in% files_chunk]
  if (nrow(sanitized_chunk) == 0) next

  TATM_raw <- gdxtools::batch_extract("TATM", files = files_chunk)$TATM
  FORC_raw <- gdxtools::batch_extract("FORCING", files = files_chunk)$FORCING
  SRM_raw <- gdxtools::batch_extract("SRM", files = files_chunk)$SRM
  if (is.null(TATM_raw) || is.null(FORC_raw)) next

  TATM <- setDT(TATM_raw)
  TATM <- merge(TATM, sanitized_chunk, by = "gdx", all = FALSE)
  TATM <- sanitize_dt(TATM)

  FORC <- setDT(FORC_raw)
  FORC <- merge(FORC, sanitized_chunk, by = "gdx", all = FALSE)
  FORC <- sanitize_dt(FORC)

  SRM_var <- data.table()
  if (!is.null(SRM_raw)) {
    SRM_var <- setDT(SRM_raw)
    SRM_var <- merge(SRM_var, sanitized_chunk, by = "gdx", all = FALSE)
    SRM_var <- sanitize_dt(SRM_var)
  }

  tot_forcing <- FORC[, .(value = sum(value)), by = c("t", "file", scenario_names, "experiment")]

  effects_pb <- build_effects_pair(
    TATM, tot_forcing, "pulse", "base", "pulse-base",
    temp_keys = c("t", "rcp", "ecs", "tcr"),
    forc_keys = c("t", "rcp", "ecs", "tcr")
  )
  effects_ss <- build_effects_pair(
    TATM, tot_forcing, "srmpulse", "srm", "srmpulse-srm",
    temp_keys = c("t", "gas", "rcp", "ecs", "tcr", "cool_rate", "pulse_time", "geo_start", "geo_end"),
    forc_keys = c("t", "gas", "rcp", "ecs", "tcr", "cool_rate", "pulse_time", "geo_start", "geo_end")
  )

  effects_chunk <- rbindlist(list(effects_pb, effects_ss), use.names = TRUE, fill = TRUE)
  if (nrow(effects_chunk) > 0) {
    keep_eff <- intersect(
      c("pair", "gas", "rcp", "ecs", "tcr", "pulse_time", "cool_rate", "geo_start", "geo_end", "term", "t", "t_rel", "Dtemp", "Dforc"),
      names(effects_chunk)
    )
    fwrite(effects_chunk[, ..keep_eff], file = effects_file, append = effects_written)
    effects_written <- TRUE
  }

  if (nrow(SRM_var) > 0) {
    srm_chunk <- SRM_var[experiment == "srmpulsemasked"]
    if (nrow(srm_chunk) > 0) {
      srm_chunk <- srm_chunk[t >= as.numeric(pulse_time)]
      srm_chunk[, `:=`(
        pair = "srmpulse-srm",
        t_rel = as.numeric(t) - as.numeric(pulse_time),
        SRMminus = -value
      )]
      keep_srm <- intersect(
        c("pair", "gas", "rcp", "ecs", "tcr", "pulse_time", "cool_rate", "geo_start", "geo_end", "term", "t", "t_rel", "SRMminus"),
        names(srm_chunk)
      )
      fwrite(srm_chunk[, ..keep_srm], file = srm_timeseries_file, append = srm_written)
      srm_written <- TRUE
    }
  }

  rm(TATM_raw, FORC_raw, SRM_raw, TATM, FORC, SRM_var, tot_forcing, effects_pb, effects_ss, effects_chunk)
  invisible(gc())
}

if (!effects_written) stop("No effects computed after chunk processing")

cat("Reloading compact effects and computing statistics...\n")
effects <- fread(effects_file)
srm_masked <- if (file.exists(srm_timeseries_file)) fread(srm_timeseries_file) else data.table()

stats_by_gas <- rbindlist(list(
  quantile_summary(effects, "Dtemp", c("gas", "t_rel")),
  quantile_summary(effects, "Dforc", c("gas", "t_rel"))
), use.names = TRUE, fill = TRUE)

id_candidates <- c("pair", "gas", "rcp", "ecs", "tcr", "pulse_time", "cool_rate", "geo_start", "geo_end", "term")
id_cols <- intersect(id_candidates, names(effects))
group_cols <- intersect(c("gas"), id_cols)
run_cols <- setdiff(id_cols, group_cols)

sel_temp <- select_percentile_runs(effects, "Dtemp", group_cols, run_cols, traj_probs)
sel_forc <- select_percentile_runs(effects, "Dforc", group_cols, run_cols, traj_probs)
selected_runs <- rbindlist(list(sel_temp$best, sel_forc$best), use.names = TRUE, fill = TRUE)
selected_traj <- rbindlist(list(sel_temp$traj, sel_forc$traj), use.names = TRUE, fill = TRUE)

srm_plot_line <- data.table()
srm_stats_by_gas <- data.table()
srm_selected_runs <- data.table()
srm_selected_traj <- data.table()

if (nrow(srm_masked) > 0) {
  id_cols_srm <- intersect(id_candidates, names(srm_masked))
  group_cols_srm <- intersect(c("gas"), id_cols_srm)
  run_cols_srm <- setdiff(id_cols_srm, group_cols_srm)
  sel_srm <- select_percentile_runs(srm_masked, "SRMminus", group_cols_srm, run_cols_srm, traj_probs)
  srm_selected_runs <- sel_srm$best
  srm_selected_traj <- sel_srm$traj
  srm_plot_line <- srm_selected_traj[t_rel < 200 & percentile == 0.5, .(gas, variable = "Dforc", t_rel, value)]
  srm_stats_by_gas <- quantile_summary(srm_masked, "SRMminus", c("gas", "t_rel"))
  srm_stats_by_gas[, variable := "Dforc"]
}

cat("Writing outputs...\n")
fwrite(stats_by_gas, file = file.path(output_folder, "pulse_effects_stats_by_gas.csv"))
fwrite(selected_traj, file = file.path(output_folder, "pulse_effects_percentile_trajectories.csv"))
if (nrow(srm_selected_runs) > 0) {
  fwrite(srm_stats_by_gas, file = file.path(output_folder, "pulse_effects_srmminus_stats_by_gas.csv"))
  fwrite(srm_selected_traj, file = file.path(output_folder, "pulse_effects_srmminus_percentile_trajectories.csv")) }

if (plot_results) {
  plot_traj <- selected_traj[
    t_rel < 200 & percentile %in% c(0.05, 0.5, 0.95),
    .(gas, variable, t_rel, percentile, value)
  ]

  plot_line <- plot_traj[percentile == 0.5]
  plot_ribbon <- dcast(
    plot_traj[percentile %in% c(0.05, 0.95)],
    gas + variable + t_rel ~ percentile,
    value.var = "value"
  )

  plot_stats <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = plot_ribbon,
      ggplot2::aes(x = t_rel, ymin = `0.05`, ymax = `0.95`, fill = gas),
      alpha = 0.25,
      color = NA
    ) +
    ggplot2::geom_line(
      data = plot_line,
      ggplot2::aes(x = t_rel, y = value, color = gas),
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      data = srm_plot_line,
      ggplot2::aes(x = t_rel, y = -value, color = gas),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(variable ~ ., scales = "free_y")
  print(plot_stats)

  # Same figure using year-by-year percentiles (not coherent trajectories).
  plot_stats_yby <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = stats_by_gas[t_rel < 200],
      ggplot2::aes(x = t_rel, ymin = p05, ymax = p95, fill = gas),
      alpha = 0.25,
      color = NA
    ) +
    ggplot2::geom_line(
      data = stats_by_gas[t_rel < 200],
      ggplot2::aes(x = t_rel, y = median, color = gas),
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      data = srm_stats_by_gas[t_rel < 200],
      ggplot2::aes(x = t_rel, y = -median, color = gas),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(variable ~ ., scales = "free_y")
  print(plot_stats_yby)
} else {
  cat("Skipping plots because --plot is set to F.\n")
}

cat("Saved outputs to", output_folder, "\n")

if (run_hpc == TRUE) {
  sink(type = "message")
  sink()
  close(logfile)
}
