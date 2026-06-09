#### Plot differences between two GDX files.
####
#### Default variables:
####   TATM,W_EMI:co2,W_EMI:ch4
####
#### Supported variable specifications:
####   TATM                  whole one-dimensional symbol indexed by t
####   W_EMI:co2             shorthand for the first non-time dimension
####   W_EMI('co2')          same shorthand, using GAMS-style member syntax
####   W_EMI[ghg=co2]        explicit dimension selector
####   FORCING[ghg=ch4]
####
#### Each comma-separated entry may combine several symbols with the operators
#### +,-,*,/ (operators are assumed never to appear inside a symbol name). The
#### arithmetic is evaluated on the experiment and base levels separately and then
#### differenced, so products and ratios are handled correctly:
####   W_EMI('co2')+CO2fromCH4    plots the difference of (W_EMI('co2') + CO2fromCH4)
####
#### Any dimension that a spec does not pin down is summed over (time is never
#### summed). Passing W_EMI alone therefore plots W_EMI summed across all ghgs.
#### In that aggregated case each gas is first converted to GtCO2-equivalent using
#### AR6 GWP-100 (and the model's Mt/Gt and N2-equivalent unit conventions), so the
#### sum is in GtCO2e/yr. A pinned gas (e.g. W_EMI:ch4) keeps its native unit.
####
#### CLI usage from the repository root:
####   Rscript Plots/gdx_difference_profiles.R --experiment=project --base=base
####
#### Experiment and base names are given WITHOUT the .gdx extension and WITHOUT a
#### directory: both the extension and the --folder location are added for you.
#### --folder (default "Results") is where inputs are read from and outputs written.
####
#### Multiple experiments can be compared against the same base scenario by
#### separating their names with a dash. Each experiment is differenced against
#### --base and drawn as its own colored scenario in every facet:
####   Rscript Plots/gdx_difference_profiles.R --experiment="a-b" --base=base --folder=Results
#### (Note: experiment names therefore cannot themselves contain a dash.)
####
#### A single experiment can override the global --base with an "experiment:base"
#### suffix; the override applies only to that experiment. For example
####   --experiment="scen1-scen2:scen3" --base=base
#### differences scen1 against base but scen2 against scen3.
#### (Experiment/base names therefore cannot contain a colon either.)
####
#### Optional CLI arguments:
####   --folder=Results              folder for inputs and outputs
####   --vars=TATM,W_EMI:co2,W_EMI:ch4
####   --tmax=100
####   --title=Experiment relative to base
####     (the variable plot is written to <folder>/gdx_difference_profiles.png)
####
#### Social value plot (in addition to the physical-variable plot):
####   --plot_social_value=yes        also draw the discounted social value plot
####   --normalization=auto           normalization factor; single value applied to
####                                  all experiments, or dash-separated one-per-experiment.
####                                  When omitted, each experiment is normalized by its own
####                                  default: the sum over time of the W_EMI difference
####                                  (experiment minus base) across all GHGs, in tonnes
####                                  CO2e (AR6 GWP-100, honouring the Mt/Gt and
####                                  N2-equivalent unit conventions).
####   --social_value_tmax=480
####   --discount_rate=0.02
####   --cumulative=yes               also show the cumulative profile (faceted)
####   --social_value_title=...
####     (the social value plot is written to <folder>/social_value_profiles.png)
#### Each experiment is coloured as its own scenario, mirroring the variable plot.

suppressPackageStartupMessages({
  library(tidyverse)
  library(gdxtools)
})

DEFAULT_VARIABLES <- c("TATM", "W_EMI:co2", "W_EMI:ch4")
HOWARD_STERNER_DAMAGE_SHARE <- 0.595 * 1.25 / 100
DEFAULT_GLOBAL_OUTPUT_USD <- 105e12

# GWP-100 from IPCC AR6 (WG1, Ch.7, Table 7.15). CH4 has separate fossil (29.8)
# and non-fossil (27.2) values; the model's W_EMI('ch4') is total methane, so a
# single value is used here (default: non-fossil, since biogenic CH4 dominates).
AR6_GWP100 <- c(co2 = 1, ch4 = 27.2, n2o = 273)

# Molar masses needed to undo the model's emission unit conventions, from
# Model/parameters.gms.
GHG_MM <- c(co2 = 44.01, ch4 = 16.04, n2o = 44.013, n2 = 28.013)

# Per-gas factor that turns a W_EMI(ghg) value into GtCO2-equivalent per year.
# Accounts for (i) AR6 GWP-100, (ii) the Mt->Gt step for the non-CO2 gases
# (W_EMI is GtCO2 but MtCH4 / MtN2O), and (iii) that W_EMI('n2o') is carried in
# N2-equivalent mass, so it is rescaled to actual N2O mass. Gases without a GWP
# (h2o, o3trop) get 0 - they also carry no emissions in this model.
wemi_co2e_factor <- function(ghg) {
  unit_to_gt <- ifelse(ghg == "co2", 1, 1e-3)
  mass_rescale <- ifelse(ghg == "n2o", GHG_MM[["n2o"]] / GHG_MM[["n2"]], 1)
  gwp <- unname(AR6_GWP100[ghg])
  gwp[is.na(gwp)] <- 0
  gwp * unit_to_gt * mass_rescale
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

parse_args <- function(args) {
  opts <- list()
  for (arg in args) {
    if (!str_starts(arg, "--") || !str_detect(arg, "=")) next
    key <- str_remove(str_extract(arg, "^--[^=]+"), "^--")
    value <- str_remove(arg, "^--[^=]+=")
    opts[[key]] <- value
  }
  opts
}

split_experiment_specs <- function(x) {
  if (is.null(x) || length(x) == 0) return(character())
  raw <- paste(x, collapse = "-")
  parts <- str_trim(strsplit(raw, "-", fixed = TRUE)[[1]])
  parts[nzchar(parts)]
}

experiment_label <- function(path) {
  str_remove(basename(path), "\\.gdx$")
}

#### Resolve a GDX name into a usable path: append the .gdx extension when it is
#### missing, and look the file up inside `folder` unless the name already carries
#### its own directory component. Idempotent, so it is safe to call more than once.
resolve_gdx <- function(name, folder = ".") {
  name <- if (str_detect(name, "\\.gdx$")) name else paste0(name, ".gdx")
  if (str_detect(name, "[/\\\\]")) name else file.path(folder, name)
}

#### Split the dash-separated experiment list and honour an optional per-experiment
#### base override written as "experiment:base" (e.g. "scen1-scen2:scen3" differences
#### scen1 against the global --base but scen2 against scen3). Returns resolved
#### experiment and base paths as parallel vectors.
resolve_experiment_bases <- function(experiment_gdx, base_gdx, folder = ".") {
  experiments <- if (length(experiment_gdx) > 1) experiment_gdx else split_experiment_specs(experiment_gdx)
  if (length(experiments) == 0) {
    stop("Provide at least one experiment GDX file.", call. = FALSE)
  }

  exp_names <- str_remove(experiments, ":.*$")
  overrides <- ifelse(str_detect(experiments, ":"), str_remove(experiments, "^[^:]*:"), NA_character_)
  default_base <- resolve_gdx(base_gdx, folder = folder)

  list(
    experiment = map_chr(exp_names, resolve_gdx, folder = folder),
    base = vapply(
      overrides,
      function(ov) if (is.na(ov)) default_base else resolve_gdx(ov, folder = folder),
      character(1)
    )
  )
}

split_variable_specs <- function(x) {
  if (is.null(x) || length(x) == 0) return(character())
  raw <- paste(x, collapse = ",")
  chars <- strsplit(raw, "", fixed = TRUE)[[1]]
  depth <- 0
  token <- ""
  out <- character()

  for (ch in chars) {
    if (ch %in% c("[", "(")) depth <- depth + 1
    if (ch %in% c("]", ")")) depth <- depth - 1

    if (ch == "," && depth == 0) {
      out <- c(out, str_trim(token))
      token <- ""
    } else {
      token <- paste0(token, ch)
    }
  }

  out <- c(out, str_trim(token))
  out[nzchar(out)]
}

strip_quotes <- function(x) {
  str_remove_all(str_trim(x), "^['\"]|['\"]$")
}

parse_selector_text <- function(selector_text) {
  selector_text <- str_trim(selector_text)
  if (!nzchar(selector_text)) return(list(named = list(), positional = character()))

  parts <- split_variable_specs(selector_text)
  named <- list()
  positional <- character()

  for (part in parts) {
    if (str_detect(part, "=")) {
      key <- str_trim(str_extract(part, "^[^=]+"))
      value <- strip_quotes(str_remove(part, "^[^=]+="))
      named[[key]] <- value
    } else {
      positional <- c(positional, strip_quotes(part))
    }
  }

  list(named = named, positional = positional)
}

parse_variable_spec <- function(spec) {
  spec <- str_trim(spec)

  bracket_match <- str_match(spec, "^([A-Za-z_][A-Za-z0-9_]*)\\[(.*)\\]$")
  if (!is.na(bracket_match[1, 1])) {
    selectors <- parse_selector_text(bracket_match[1, 3])
    return(list(raw = spec, symbol = bracket_match[1, 2], selectors = selectors))
  }

  call_match <- str_match(spec, "^([A-Za-z_][A-Za-z0-9_]*)\\((.*)\\)$")
  if (!is.na(call_match[1, 1])) {
    selectors <- parse_selector_text(call_match[1, 3])
    return(list(raw = spec, symbol = call_match[1, 2], selectors = selectors))
  }

  colon_match <- str_match(spec, "^([A-Za-z_][A-Za-z0-9_]*):(.+)$")
  if (!is.na(colon_match[1, 1])) {
    selectors <- list(named = list(), positional = strip_quotes(colon_match[1, 3]))
    return(list(raw = spec, symbol = colon_match[1, 2], selectors = selectors))
  }

  list(raw = spec, symbol = spec, selectors = list(named = list(), positional = character()))
}

choose_variable_specs <- function(variables = NULL,
                                  default = DEFAULT_VARIABLES,
                                  ask = interactive()) {
  if (!is.null(variables) && length(variables) > 0) {
    return(split_variable_specs(variables))
  }

  if (ask) {
    prompt <- paste0(
      "Variables to plot [",
      paste(default, collapse = ","),
      "]: "
    )
    answer <- readline(prompt)
    if (nzchar(str_trim(answer))) return(split_variable_specs(answer))
  }

  default
}

as_source <- function(gdx, experiment_gdx, base_gdx) {
  case_when(
    normalizePath(gdx, winslash = "/", mustWork = FALSE) ==
      normalizePath(experiment_gdx, winslash = "/", mustWork = FALSE) ~ "experiment",
    normalizePath(gdx, winslash = "/", mustWork = FALSE) ==
      normalizePath(base_gdx, winslash = "/", mustWork = FALSE) ~ "base",
    TRUE ~ NA_character_
  )
}

apply_selectors <- function(data, spec) {
  selector_cols <- setdiff(names(data), c("gdx", "source", "value", "t"))

  for (key in names(spec$selectors$named)) {
    if (!key %in% names(data)) {
      stop(
        "Variable specification '", spec$raw, "' refers to dimension '", key,
        "', but available dimensions are: ", paste(selector_cols, collapse = ", "),
        call. = FALSE
      )
    }
    data <- data %>% filter(.data[[key]] == spec$selectors$named[[key]])
  }

  if (length(spec$selectors$positional) > 0) {
    if (length(selector_cols) < length(spec$selectors$positional)) {
      stop(
        "Variable specification '", spec$raw, "' has too many positional selectors.",
        call. = FALSE
      )
    }
    for (i in seq_along(spec$selectors$positional)) {
      key <- selector_cols[[i]]
      data <- data %>% filter(.data[[key]] == spec$selectors$positional[[i]])
    }
  }

  data
}

format_series_name <- function(spec) {
  named_values <- unlist(spec$selectors$named, use.names = FALSE)
  positional_values <- spec$selectors$positional
  selected_values <- c(named_values, positional_values)

  if (length(selected_values) == 0) return(spec$symbol)
  paste0(spec$symbol, "(", paste(selected_values, collapse = ","), ")")
}

selected_member <- function(spec, dimension_name = NULL) {
  if (!is.null(dimension_name) && dimension_name %in% names(spec$selectors$named)) {
    return(spec$selectors$named[[dimension_name]])
  }
  if (length(spec$selectors$positional) > 0) return(spec$selectors$positional[[1]])
  NA_character_
}

gdx_symbol_unit <- function(spec) {
  symbol <- spec$symbol
  member <- selected_member(spec)

  if (symbol %in% c("TATM", "TSLOW", "TFAST", "target_temp")) return("deg C")

  if (symbol == "W_EMI") {
    return(case_when(
      member == "co2" ~ "GtCO2/yr",
      member == "ch4" ~ "MtCH4/yr",
      member == "n2o" ~ "MtN2O/yr",
      TRUE ~ "GtCO2e/yr"
    ))
  }

  if (symbol %in% c("FORCING", "forcing_exogenous", "forcing_srm", "SRM")) {
    return("W/m2")
  }

  if (symbol == "CONC") {
    return(case_when(
      member == "co2" ~ "ppm",
      member %in% c("ch4", "n2o") ~ "ppb",
      TRUE ~ "concentration"
    ))
  }

  if (symbol %in% c("CO2fromCH4", "CO2toCH4")) return("GtC/yr")
  if (symbol %in% c("CUMEMI", "C_SINKS", "C_ATM")) return("GtC")
  if (symbol == "IRF") return("yr")
  if (symbol == "FF_CH4") return("fraction")

  ""
}

format_facet_label <- function(series, unit) {
  ifelse(nzchar(unit), paste0(series, " [", unit, "]"), series)
}

#### Split a single variable token into operands and the +,-,*,/ operators that
#### join them, ignoring operators inside [] or () (selectors). Operators are
#### assumed never to appear inside variable names.
split_expression <- function(token) {
  chars <- strsplit(token, "", fixed = TRUE)[[1]]
  depth <- 0
  operands <- character()
  operators <- character()
  current <- ""

  for (ch in chars) {
    if (ch %in% c("[", "(")) depth <- depth + 1
    if (ch %in% c("]", ")")) depth <- depth - 1

    if (depth == 0 && ch %in% c("+", "-", "*", "/")) {
      operands <- c(operands, str_trim(current))
      operators <- c(operators, ch)
      current <- ""
    } else {
      current <- paste0(current, ch)
    }
  }

  operands <- c(operands, str_trim(current))
  list(operands = operands[nzchar(operands)], operators = operators)
}

#### Combine aligned operand vectors with their operators, honouring the usual
#### precedence (* and / before + and -). `values` is a list of numeric vectors,
#### `operators` a character vector of length length(values) - 1.
apply_arithmetic <- function(values, operators) {
  if (length(values) == 1) return(values[[1]])
  if (length(operators) != length(values) - 1) {
    stop("Malformed arithmetic expression across variables.", call. = FALSE)
  }

  vals <- values
  ops <- operators
  i <- 1
  while (i <= length(ops)) {
    if (ops[i] %in% c("*", "/")) {
      vals[[i]] <- if (ops[i] == "*") vals[[i]] * vals[[i + 1]] else vals[[i]] / vals[[i + 1]]
      vals[[i + 1]] <- NULL
      ops <- ops[-i]
    } else {
      i <- i + 1
    }
  }

  result <- vals[[1]]
  for (j in seq_along(ops)) {
    result <- if (ops[j] == "+") result + vals[[j + 1]] else result - vals[[j + 1]]
  }
  result
}

#### Extract a single symbol and return its experiment/base levels per t, summing
#### over any dimension that the spec does not pin down (never over time).
extract_symbol_levels <- function(spec, experiment_gdx, base_gdx) {
  files <- c(experiment_gdx, base_gdx)
  extracted <- gdxtools::batch_extract(spec$symbol, files = files)[[spec$symbol]]

  if (is.null(extracted)) {
    stop("Could not extract symbol '", spec$symbol, "' from the GDX files.", call. = FALSE)
  }

  data <- extracted %>%
    as_tibble() %>%
    mutate(
      t = as.numeric(t),
      source = as_source(gdx, experiment_gdx, base_gdx)
    )

  data <- apply_selectors(data, spec)

  if (!"t" %in% names(data)) {
    stop("Symbol '", spec$symbol, "' is not indexed by t.", call. = FALSE)
  }

  # When W_EMI is aggregated across gases (no member pinned), convert each gas to
  # GtCO2-equivalent with its AR6 GWP before summing, so the total is meaningful.
  summed_over_ghg <- length(spec$selectors$named) == 0 &&
    length(spec$selectors$positional) == 0
  if (spec$symbol == "W_EMI" && summed_over_ghg && "ghg" %in% names(data)) {
    data <- data %>% mutate(value = value * wemi_co2e_factor(ghg))
  }

  data %>%
    group_by(source, t) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    pivot_wider(names_from = source, values_from = value)
}

#### Extract one comma-separated token, which may combine several symbols with
#### +,-,*,/. The arithmetic is evaluated on experiment and base levels
#### separately, then differenced, so that products and ratios are handled
#### correctly (not just sums).
extract_expression_delta <- function(token, experiment_gdx, base_gdx,
                                     experiment = experiment_label(experiment_gdx),
                                     year_offset = 2019) {
  parsed <- split_expression(token)
  specs <- map(parsed$operands, parse_variable_spec)

  levels_list <- map(specs, extract_symbol_levels, experiment_gdx, base_gdx)

  joined <- reduce(seq_along(levels_list), function(acc, i) {
    lv <- levels_list[[i]] %>%
      rename(!!paste0("exp_", i) := experiment, !!paste0("base_", i) := base)
    if (is.null(acc)) lv else inner_join(acc, lv, by = "t")
  }, .init = NULL)

  exp_cols <- map(seq_along(specs), ~ joined[[paste0("exp_", .x)]])
  base_cols <- map(seq_along(specs), ~ joined[[paste0("base_", .x)]])

  units <- map_chr(specs, gdx_symbol_unit)
  unit <- if (length(unique(units)) == 1) units[[1]] else ""
  series <- if (length(specs) == 1) format_series_name(specs[[1]]) else token

  tibble(
    t = joined$t,
    year = year_offset + joined$t,
    scenario = experiment,
    base = experiment_label(base_gdx),
    variable = token,
    series = series,
    unit = unit,
    facet_label = format_facet_label(series, unit),
    value = apply_arithmetic(exp_cols, parsed$operators) -
      apply_arithmetic(base_cols, parsed$operators)
  ) %>%
    filter(!is.na(value))
}

build_gdx_difference_data <- function(experiment_gdx,
                                      base_gdx,
                                      folder = ".",
                                      variables = NULL,
                                      tmax = 100,
                                      year_offset = 2019,
                                      ask_variables = interactive()) {
  eb <- resolve_experiment_bases(experiment_gdx, base_gdx, folder = folder)
  scenario_levels <- experiment_label(eb$experiment)

  tokens <- choose_variable_specs(variables, ask = ask_variables)

  plot_data <- map2_dfr(eb$experiment, eb$base, function(exp_path, base_path) {
    tokens %>%
      map_dfr(
        extract_expression_delta,
        experiment_gdx = exp_path,
        base_gdx = base_path,
        year_offset = year_offset
      )
  })

  plot_data %>%
    filter(t <= tmax) %>%
    mutate(
      scenario = factor(scenario, levels = unique(scenario_levels)),
      series = factor(series, levels = unique(series)),
      facet_label = factor(facet_label, levels = unique(facet_label))
    )
}

plot_gdx_difference <- function(experiment_gdx,
                                base_gdx,
                                folder = ".",
                                variables = NULL,
                                output_file = NULL,
                                tmax = 100,
                                year_offset = 2019,
                                title = "Experiment relative to base",
                                ask_variables = interactive()) {
  plot_data <- build_gdx_difference_data(
    experiment_gdx = experiment_gdx,
    base_gdx = base_gdx,
    folder = folder,
    variables = variables,
    tmax = tmax,
    year_offset = year_offset,
    ask_variables = ask_variables
  )

  if (nrow(plot_data) == 0) {
    stop("No data available after applying variable and tmax filters.", call. = FALSE)
  }

  base_map <- plot_data %>% distinct(scenario, base)
  caption_base <- if (n_distinct(base_map$base) == 1) {
    paste0("Base: ", base_map$base[[1]])
  } else {
    paste0(
      "Base per scenario: ",
      paste0(base_map$scenario, " vs ", base_map$base, collapse = "; ")
    )
  }

  p <- ggplot(plot_data, aes(x = year, y = value, color = scenario, group = scenario)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
    geom_line(linewidth = 0.9) +
    facet_wrap(~facet_label, ncol = 1, scales = "free_y") +
    labs(
      title = title,
      subtitle = "Difference is experiment minus base",
      x = NULL,
      y = NULL,
      color = NULL,
      caption = paste0(
        "Experiments: ", paste(levels(plot_data$scenario), collapse = ", "),
        "\n", caption_base
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),
      plot.caption = element_text(hjust = 0, size = 8)
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(output_file, p, width = 8, height = max(4, 2.2 * length(unique(plot_data$facet_label))), dpi = 300)
    message("Wrote ", output_file)
  }

  invisible(list(plot = p, data = plot_data))
}

extract_temperature_pair <- function(experiment_gdx,
                                     base_gdx,
                                     year_offset = 2019) {
  gdxtools::batch_extract("TATM", files = c(experiment_gdx, base_gdx))$TATM %>%
    as_tibble() %>%
    mutate(
      t = as.numeric(t),
      year = year_offset + t,
      source = as_source(gdx, experiment_gdx, base_gdx)
    ) %>%
    select(source, t, year, value) %>%
    pivot_wider(names_from = source, values_from = value) %>%
    filter(!is.na(experiment), !is.na(base)) %>%
    rename(Texp = experiment, Tbase = base) %>%
    arrange(t)
}

evaluate_strategy_social_value <- function(experiment_gdx,
                                           base_gdx,
                                           normalization_factor,
                                           discount_rate = 0.02,
                                           tmax = 480,
                                           discount_start_t = 1,
                                           year_offset = 2019,
                                           damage_share = HOWARD_STERNER_DAMAGE_SHARE,
                                           global_output = DEFAULT_GLOBAL_OUTPUT_USD,
                                           output_growth = 0.015,
                                           a = NULL) {
  if (missing(normalization_factor) || is.null(normalization_factor) ||
      is.na(normalization_factor) || normalization_factor == 0) {
    stop("Provide a non-zero normalization_factor.", call. = FALSE)
  }

  temp <- extract_temperature_pair(
    experiment_gdx = experiment_gdx,
    base_gdx = base_gdx,
    year_offset = year_offset
  ) %>%
    filter(t <= tmax)

  if (nrow(temp) == 0) {
    stop("No TATM data available after applying tmax.", call. = FALSE)
  }

  temp <- temp %>%
    mutate(
      periods_from_discount_start = t - discount_start_t,
      output = global_output * (1 + output_growth)^periods_from_discount_start,
      damage_scale = if (is.null(a)) output * damage_share else a,
      damage_experiment = damage_scale * Texp^2,
      damage_base = damage_scale * Tbase^2,
      damage_difference = damage_experiment - damage_base,
      social_value = damage_base - damage_experiment,
      discount_factor = 1 / (1 + discount_rate)^periods_from_discount_start,
      discounted_damage_difference = damage_difference * discount_factor,
      discounted_social_value = social_value * discount_factor,
      discounted_social_value_per_unit = discounted_social_value / normalization_factor,
      cumulative_social_value_per_unit = cumsum(discounted_social_value_per_unit)
    )

  total_damage_difference <- sum(temp$discounted_damage_difference, na.rm = TRUE)
  total_social_value <- sum(temp$discounted_social_value, na.rm = TRUE)

  list(
    summary = tibble(
      social_value_per_unit = total_social_value / normalization_factor,
      npv_damage_difference_per_unit = total_damage_difference / normalization_factor,
      social_value_total = total_social_value,
      npv_damage_difference_total = total_damage_difference,
      normalization_factor = normalization_factor,
      discount_rate = discount_rate,
      discount_start_t = discount_start_t,
      tmax = tmax,
      damage_share = if (is.null(a)) damage_share else NA_real_,
      global_output = if (is.null(a)) global_output else NA_real_,
      output_growth = if (is.null(a)) output_growth else NA_real_,
      a = if (is.null(a)) global_output * damage_share else a
    ),
    time_profile = temp %>%
      select(
        t,
        year,
        Texp,
        Tbase,
        damage_experiment,
        damage_base,
        damage_difference,
        social_value,
        discount_factor,
        discounted_damage_difference,
        discounted_social_value,
        discounted_social_value_per_unit,
        cumulative_social_value_per_unit
      )
  )
}

plot_strategy_social_value <- function(social_value_result,
                                       output_file = NULL,
                                       cumulative = TRUE,
                                       title = "Discounted social value profile") {
  profile <- social_value_result$time_profile

  if (is.null(profile) || nrow(profile) == 0) {
    stop("social_value_result must contain a non-empty time_profile.", call. = FALSE)
  }

  plot_data <- profile %>%
    select(
      year,
      discounted_social_value_per_unit,
      cumulative_social_value_per_unit
    ) %>%
    pivot_longer(
      cols = c(discounted_social_value_per_unit, cumulative_social_value_per_unit),
      names_to = "series",
      values_to = "value"
    ) %>%
    filter(cumulative | series == "discounted_social_value_per_unit") %>%
    mutate(
      series = recode(
        series,
        discounted_social_value_per_unit = "Annual discounted social value",
        cumulative_social_value_per_unit = "Cumulative discounted social value"
      )
    )

  p <- ggplot(plot_data, aes(x = year, y = value, color = series)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
    geom_line(linewidth = 0.9) +
    labs(
      title = title,
      subtitle = "Values are normalized by the supplied project unit",
      x = NULL,
      y = "Net-present dollars per unit",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(output_file, p, width = 8, height = 4.8, dpi = 300)
    message("Wrote ", output_file)
  }

  invisible(p)
}

plot_strategy_social_value_from_gdx <- function(experiment_gdx,
                                                base_gdx,
                                                normalization_factor,
                                                output_file = NULL,
                                                discount_rate = 0.02,
                                                tmax = 480,
                                                discount_start_t = 1,
                                                year_offset = 2025,
                                                damage_share = HOWARD_STERNER_DAMAGE_SHARE,
                                                global_output = DEFAULT_GLOBAL_OUTPUT_USD,
                                                output_growth = 0,
                                                a = NULL,
                                                cumulative = FALSE,
                                                title = "Discounted social value profile") {
  result <- evaluate_strategy_social_value(
    experiment_gdx = experiment_gdx,
    base_gdx = base_gdx,
    normalization_factor = normalization_factor,
    discount_rate = discount_rate,
    tmax = tmax,
    discount_start_t = discount_start_t,
    year_offset = year_offset,
    damage_share = damage_share,
    global_output = global_output,
    output_growth = output_growth,
    a = a
  )

  plot <- plot_strategy_social_value(
    social_value_result = result,
    output_file = output_file,
    cumulative = cumulative,
    title = title
  )

  invisible(list(plot = plot, summary = result$summary, time_profile = result$time_profile))
}

#### Default normalization for one experiment/base pair: the sum over time of the
#### W_EMI difference (experiment minus base) aggregated across all GHGs and
#### expressed in tonnes of CO2-equivalent (AR6 GWP-100, honouring the model's
#### Mt/Gt and N2-equivalent unit conventions). This mirrors the aggregated W_EMI
#### series of the variable plot (which is in GtCO2e) but is converted to tCO2e, so
#### the social value is reported per tonne of cumulative emission change. The sign
#### follows the "experiment minus base" convention.
GTCO2E_PER_TCO2E <- 1e9

default_emissions_normalization <- function(experiment_gdx, base_gdx) {
  delta <- extract_expression_delta("W_EMI", experiment_gdx, base_gdx)
  total <- sum(delta$value, na.rm = TRUE) * GTCO2E_PER_TCO2E
  if (!is.finite(total) || total == 0) {
    stop(
      "Default emissions normalization is zero or undefined for experiment '",
      experiment_label(experiment_gdx), "' against base '",
      experiment_label(base_gdx), "'; pass an explicit --normalization.",
      call. = FALSE
    )
  }
  total
}

build_social_value_data <- function(experiment_gdx,
                                    base_gdx,
                                    folder = ".",
                                    normalization_factor = NULL,
                                    discount_rate = 0.02,
                                    tmax = 480,
                                    discount_start_t = 1,
                                    year_offset = 2025,
                                    damage_share = HOWARD_STERNER_DAMAGE_SHARE,
                                    global_output = DEFAULT_GLOBAL_OUTPUT_USD,
                                    output_growth = 0,
                                    a = NULL) {
  eb <- resolve_experiment_bases(experiment_gdx, base_gdx, folder = folder)
  experiments <- eb$experiment

  norm <- if (is.null(normalization_factor)) {
    map2_dbl(experiments, eb$base, default_emissions_normalization)
  } else if (length(normalization_factor) == 1) {
    rep(normalization_factor, length(experiments))
  } else if (length(normalization_factor) == length(experiments)) {
    normalization_factor
  } else {
    stop(
      "normalization_factor must be a single value or one value per experiment (",
      length(experiments), " given).",
      call. = FALSE
    )
  }

  scenario_levels <- experiment_label(experiments)

  results <- pmap(list(experiments, eb$base, norm), function(experiment_gdx, base_path, nf) {
    evaluate_strategy_social_value(
      experiment_gdx = experiment_gdx,
      base_gdx = base_path,
      normalization_factor = nf,
      discount_rate = discount_rate,
      tmax = tmax,
      discount_start_t = discount_start_t,
      year_offset = year_offset,
      damage_share = damage_share,
      global_output = global_output,
      output_growth = output_growth,
      a = a
    )
  })

  time_profile <- map2_dfr(results, scenario_levels, function(res, scenario) {
    res$time_profile %>% mutate(scenario = scenario, .before = 1)
  }) %>%
    mutate(scenario = factor(scenario, levels = unique(scenario_levels)))

  summary <- map2_dfr(results, scenario_levels, function(res, scenario) {
    res$summary %>% mutate(scenario = scenario, .before = 1)
  })

  list(time_profile = time_profile, summary = summary)
}

plot_social_value_scenarios <- function(social_value_data,
                                        output_file = NULL,
                                        cumulative = FALSE,
                                        title = "Discounted social value profile") {
  profile <- social_value_data$time_profile

  if (is.null(profile) || nrow(profile) == 0) {
    stop("social_value_data must contain a non-empty time_profile.", call. = FALSE)
  }

  plot_data <- profile %>%
    select(
      scenario,
      year,
      discounted_social_value_per_unit,
      cumulative_social_value_per_unit
    ) %>%
    pivot_longer(
      cols = c(discounted_social_value_per_unit, cumulative_social_value_per_unit),
      names_to = "metric",
      values_to = "value"
    ) %>%
    filter(cumulative | metric == "discounted_social_value_per_unit") %>%
    mutate(
      metric = recode(
        metric,
        discounted_social_value_per_unit = "Annual discounted social value",
        cumulative_social_value_per_unit = "Cumulative discounted social value"
      ),
      metric = factor(
        metric,
        levels = c("Annual discounted social value", "Cumulative discounted social value")
      )
    )

  p <- ggplot(plot_data, aes(x = year, y = value, color = scenario, group = scenario)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
    geom_line(linewidth = 0.9) +
    labs(
      title = title,
      subtitle = "Social value is base damages minus experiment damages, discounted to present value",
      x = NULL,
      y = "Net-present dollars per unit",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0)
    )

  if (cumulative) {
    p <- p + facet_wrap(~metric, ncol = 1, scales = "free_y")
  }

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    n_panels <- if (cumulative) 2 else 1
    ggsave(output_file, p, width = 8, height = max(4, 2.4 * n_panels), dpi = 300)
    message("Wrote ", output_file)
  }

  invisible(list(plot = p, data = plot_data))
}

#### -------------------------------------------------------------------------
#### Emission metrics: GWP, GTP and relative social cost
#### -------------------------------------------------------------------------
#### Compute pulse emission metrics for one or more `experiment` scenarios
#### relative to a `co2` reference pulse, both taken as the marginal difference
#### from a common `base` scenario. `experiment_gdx` may be a character vector or
#### a single dash-separated string ("scen1-scen2"); every experiment is compared
#### against the same base and the same co2 reference, and the result carries a
#### `scenario` column identifying each. For each time horizon (default 20/50/100
#### yr after the pulse) the function returns:
####
####   GWP(H) = (integral of the experiment forcing difference over [pulse,
####            pulse+H]) / (same integral for the co2 pulse).
####            Forcing is the total radiative forcing summed across all gases
####            (so indirect effects - CO2 from CH4 oxidation, ozone,
####            stratospheric water vapour - are included), in W/m2.
####
####   GTP(H) = (experiment temperature difference at year pulse+H) /
####            (co2 temperature difference at the same year). End-point
####            temperature ratio, the conventional GTP definition.
####
####   social_cost_fraction(H) = (net present value of the experiment's climate
####            damage difference accumulated to pulse+H) / (same for the co2
####            pulse). The quadratic-damage / discounting machinery is shared
####            with evaluate_strategy_social_value(). An additional
####            social_cost_fraction_infinite column gives the same ratio over
####            the full model horizon (`infinite_tmax`); it is constant across
####            horizons.
####
#### Forcing/temperature deltas are reused from extract_expression_delta(), which
#### evaluates "scenario minus base" on the summed FORCING and on TATM. The
#### integral over time uses the model time step (tstep, 1 yr) as the quadrature
#### weight (rectangle rule); pre-pulse differences are ~0 so the window start is
#### immaterial. The social-cost ratios are invariant to discount_start_t, since
#### it scales numerator and denominator identically.
compute_emission_metrics <- function(experiment_gdx,
                                      co2_gdx,
                                      base_gdx,
                                      folder = ".",
                                      horizons = c(20, 50, 100),
                                      pulse_time = 5,
                                      tstep = 1,
                                      year_offset = 2019,
                                      infinite_tmax = 500,
                                      discount_rate = 0.02,
                                      discount_start_t = NULL,
                                      damage_share = HOWARD_STERNER_DAMAGE_SHARE,
                                      global_output = DEFAULT_GLOBAL_OUTPUT_USD,
                                      output_growth = 0,
                                      a = NULL) {
  # `experiment_gdx` may be a character vector or a single dash-separated string
  # (e.g. "scen1-scen2"); every experiment is compared against the same base and
  # the same co2 reference pulse.
  experiments <- if (length(experiment_gdx) > 1) experiment_gdx else split_experiment_specs(experiment_gdx)
  if (length(experiments) == 0) {
    stop("Provide at least one experiment GDX file.", call. = FALSE)
  }
  exp_paths <- map_chr(experiments, resolve_gdx, folder = folder)
  co2_path <- resolve_gdx(co2_gdx, folder = folder)
  base_path <- resolve_gdx(base_gdx, folder = folder)
  discount_start_t <- discount_start_t %||% pulse_time

  # Integral of a delta profile over the closed window [t_start, t_end], using
  # tstep as the quadrature weight.
  integrate_window <- function(df, t_start, t_end) {
    d <- df %>% filter(t >= t_start, t <= t_end)
    sum(d$value, na.rm = TRUE) * tstep
  }

  # Value of a delta profile at a single time point (linear interpolation when
  # the exact t is absent; with tstep = 1 the point is normally present).
  value_at <- function(df, t_point) {
    hit <- df$value[df$t == t_point]
    if (length(hit) > 0) return(hit[[1]])
    stats::approx(df$t, df$value, xout = t_point, rule = 2)$y
  }

  # Discounted, cumulative climate-damage difference (scenario minus base) using
  # the shared quadratic-damage convention. cumulative[t] is the NPV of damages
  # accumulated up to and including t.
  damage_profile <- function(scen_gdx) {
    extract_temperature_pair(scen_gdx, base_path, year_offset = year_offset) %>%
      filter(t <= infinite_tmax) %>%
      arrange(t) %>%
      mutate(
        periods_from_discount_start = t - discount_start_t,
        output = global_output * (1 + output_growth)^periods_from_discount_start,
        damage_scale = if (is.null(a)) output * damage_share else a,
        damage_difference = damage_scale * (Texp^2 - Tbase^2),
        discount_factor = 1 / (1 + discount_rate)^periods_from_discount_start,
        discounted_damage_difference = damage_difference * discount_factor,
        cumulative_npv = cumsum(discounted_damage_difference)
      )
  }

  cumulative_npv_to <- function(df, t_end) {
    d <- df %>% filter(t <= t_end)
    if (nrow(d) == 0) return(NA_real_)
    tail(d$cumulative_npv, 1)
  }

  # The co2 reference pulse and base are shared across experiments, so the co2
  # forcing/temperature/damage quantities are computed only once.
  forcing_co2 <- extract_expression_delta("FORCING", co2_path, base_path) %>% arrange(t)
  temp_co2 <- extract_expression_delta("TATM", co2_path, base_path) %>% arrange(t)
  damage_co2 <- damage_profile(co2_path)
  sc_co2_infinite <- cumulative_npv_to(damage_co2, infinite_tmax)

  map2_dfr(exp_paths, experiments, function(exp_path, experiment_name) {
    # Marginal differences from the common base for this experiment.
    forcing_exp <- extract_expression_delta("FORCING", exp_path, base_path) %>% arrange(t)
    temp_exp <- extract_expression_delta("TATM", exp_path, base_path) %>% arrange(t)
    damage_exp <- damage_profile(exp_path)
    sc_exp_infinite <- cumulative_npv_to(damage_exp, infinite_tmax)

    map_dfr(horizons, function(H) {
      t_end <- pulse_time + H

      agwp_exp <- integrate_window(forcing_exp, pulse_time, t_end)
      agwp_co2 <- integrate_window(forcing_co2, pulse_time, t_end)

      dT_exp <- value_at(temp_exp, t_end)
      dT_co2 <- value_at(temp_co2, t_end)

      sc_exp <- cumulative_npv_to(damage_exp, t_end)
      sc_co2 <- cumulative_npv_to(damage_co2, t_end)

      tibble(
        scenario = experiment_label(exp_path),
        horizon = H,
        pulse_time = pulse_time,
        end_t = t_end,
        end_year = year_offset + t_end,
        agwp_experiment = agwp_exp,
        agwp_co2 = agwp_co2,
        GWP = agwp_exp / agwp_co2,
        agtp_experiment = dT_exp,
        agtp_co2 = dT_co2,
        GTP = dT_exp / dT_co2,
        social_cost_experiment = sc_exp,
        social_cost_co2 = sc_co2,
        social_cost_fraction = sc_exp / sc_co2,
        social_cost_fraction_infinite = sc_exp_infinite / sc_co2_infinite
      )
    })
  })
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  if (is.null(args$experiment) || is.null(args$base)) {
    stop(
      "Provide both --experiment=experiment_name and --base=base_name ",
      "(the .gdx extension and --folder are added for you).",
      call. = FALSE
    )
  }

  folder <- args$folder %||% "Results"

  plot_gdx_difference(
    experiment_gdx = args$experiment,
    base_gdx = args$base,
    folder = folder,
    variables = args$vars,
    output_file = file.path(folder, "gdx_difference_profiles.png"),
    tmax = as.numeric(args$tmax %||% "100"),
    title = args$title %||% "Experiment relative to base",
    ask_variables = is.null(args$vars) && interactive()
  )

  if (tolower(args$plot_social_value %||% "no") %in% c("yes", "true", "1")) {
    # When --normalization is omitted (or "auto"), each experiment is normalized by
    # its own cumulative GtCO2e emission difference against its base.
    normalization <- if (is.null(args$normalization) ||
                         tolower(args$normalization) == "auto") {
      NULL
    } else {
      vals <- as.numeric(split_experiment_specs(args$normalization))
      if (length(vals) == 0 || any(is.na(vals))) {
        stop(
          "--normalization must be numeric (a single value, or dash-separated values, ",
          "one per experiment), or \"auto\".",
          call. = FALSE
        )
      }
      vals
    }

    social_value_data <- build_social_value_data(
      experiment_gdx = args$experiment,
      base_gdx = args$base,
      folder = folder,
      normalization_factor = normalization,
      discount_rate = as.numeric(args$discount_rate %||% "0.02"),
      tmax = as.numeric(args$social_value_tmax %||% "480")
    )

    print(social_value_data$summary)

    message("Normalization factor (tCO2e, experiment minus base):")
    norm_table <- social_value_data$summary %>%
      select(scenario, normalization_factor)
    print(norm_table)

    plot_social_value_scenarios(
      social_value_data = social_value_data,
      output_file = file.path(folder, "social_value_profiles.png"),
      cumulative = tolower(args$cumulative %||% "no") %in% c("yes", "true", "1"),
      title = args$social_value_title %||% "Discounted social value profile"
    )
  }
}
