
library(tidyverse)
require(patchwork)

# Single control file for plotting settings / result folders. Resolve relative to
# this script if launched via Rscript, else assume the project working directory.
.sp <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- if (length(.sp) == 1) dirname(.sp) else getwd()
while (!file.exists(file.path(.root, "all_parameters.R")) && dirname(.root) != .root)
  .root <- dirname(.root)
if (!file.exists(file.path(.root, "all_parameters.R")))
  stop("Cannot locate all_parameters.R (Paper_SAI root).")
source(file.path(.root, "all_parameters.R"))

pow10_labels <- scales::label_math(10^.x)
output_folder <- RESULTS_FOLDER_MAIN
damnpv <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "npc_output")), read.csv))
scc <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "sccnosrm_output")), read.csv)) %>% rename(scc=scc_nosrm)
scc_srm <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "scc_output")), read.csv)) %>% rename(scc_srm=scc)
all_cols <- names(damnpv)[1:20]
remove_outliers <- FIG_OUTLIER_COLS
check_densities <- damnpv %>%
  filter(gas=="co2") %>%
  filter_at(remove_outliers, ~ (.x <= quantile(.x, FIG_OUTLIER_QHI, na.rm = TRUE) & .x >= quantile(.x, FIG_OUTLIER_QLO, na.rm = TRUE) ) ) %>%
  select_at(setdiff(all_cols,c("gas") ))
check_densities %>% 
  select_at(remove_outliers) %>% 
  pivot_longer(remove_outliers) %>% 
  ggplot() +
  geom_density(aes(x=value) ) +
  facet_wrap(name~.,scales="free")
  
  
# filter scc and npc>0 and scenarios with both CH4 and CO2
damnpv <- damnpv %>% 
  inner_join(check_densities) %>% 
#  inner_join(scc %>% select(-damnpv,-ozpnpv,-pulse_time)) %>% 
#  inner_join(scc_srm %>% select(-damnpv,-ozpnpv,-pulse_time)) %>% 
  filter(npc_srm>0 ) %>% 
  group_by_at(setdiff(all_cols,c("gas")) ) %>% 
  filter(n()==2 & !any(duplicated(gas))) %>% 
  ungroup() %>% unique()

damnorm <- damnpv %>% unique() %>% 
  mutate(geo_start=ifelse(cool_rate==0,2500,geo_start),
         geo_end =ifelse(cool_rate==0,max(geo_end)+100,geo_end )) %>% 
  mutate_if(is.integer, as.numeric) %>%
  group_by(gas) %>% 
  mutate(npc_std=(npc_srm-median(npc_srm,na.rm=TRUE))/mad(npc_srm,na.rm=TRUE),
         npc_norm=npc_srm/median(npc_srm,na.rm=TRUE))


tail_inputs <- tibble::tribble(
  ~input, ~label, ~column,
  "discount", "Discount", "discount",
  "sai_angle", "SAI angle", "sai_angle",
  "ecs_k", "ECS", "ecs_k",
  "damage_alpha", "Damages", "damage_alpha",
  "termination_year", "Term. year", "termination_year"
)

density_tail_inputs <- tibble::tribble(
  ~gas, ~input, ~label, ~column,
  "co2", "damage_alpha", "Damages", "damage_alpha",
  "co2", "sai_angle", "SAI angle", "sai_angle",
  "co2", "discount", "Discount", "discount",
  "co2", "ecs_k", "ECS", "ecs_k",
  "ch4", "ozone_mortality", "Ozone mortality", "ozone_mortality",
  "ch4", "pulse_year", "Pulse year", "pulse_year",
  "ch4", "sai_angle", "SAI angle", "sai_angle",
  "ch4", "vsl_m", "VSL", "vsl_m"
)

prefix_bins <- function(values) {
  sprintf("%02d: %s", seq_along(values), values)
}

strip_bin_prefix <- function(values) {
  sub("^\\d+:\\s*", "", values)
}

bin_tail_input <- function(x, max_bins = 5) {
  ux <- sort(unique(x[!is.na(x)]))
  if (length(ux) <= 8) {
    raw_labels <- if (max(abs(ux), na.rm = TRUE) >= 100) {
      sprintf("%.0f", ux)
    } else {
      sprintf("%.1f", ux)
    }
    ordered_labels <- prefix_bins(raw_labels)
    label_map <- setNames(ordered_labels, raw_labels)
    raw_values <- if (max(abs(ux), na.rm = TRUE) >= 100) {
      sprintf("%.0f", x)
    } else {
      sprintf("%.1f", x)
    }
    return(unname(label_map[raw_values]))
  }

  cut_values <- as.character(ggplot2::cut_number(x, n = min(max_bins, length(ux)), dig.lab = 4))
  raw_labels <- unique(cut_values[order(x)])
  label_map <- setNames(prefix_bins(raw_labels), raw_labels)
  unname(label_map[cut_values])
}

tail_data <- damnorm %>%
  ungroup() %>%
  mutate(discount = delta * FIG2_DELTA_PCT,
         sai_angle = theta,
         ecs_k = ecs / FIG2_ECS_TENTHS,
         damage_alpha = alpha * FIG2_ALPHA_PCT,
         termination_year = FIG2_TERM_YEAR_BASE + term,
         ozone_mortality = mortality_ozone,
         pulse_year = FIG2_TERM_YEAR_BASE + pulse_time,
         vsl_m = vsl / 1e6) %>%
  group_by(gas) %>%
  mutate(bad_tail = npc_norm >= quantile(npc_norm, FIG2_TAIL_QUANTILE, na.rm = TRUE)) %>%
  ungroup()

tail_input_order <- tail_inputs$input
tail_input_labels <- tail_inputs$label
names(tail_input_labels) <- tail_inputs$input
tail_input_columns <- tail_inputs$column
names(tail_input_columns) <- tail_inputs$input
tail_input_labels[density_tail_inputs$input] <- density_tail_inputs$label
tail_input_columns[density_tail_inputs$input] <- density_tail_inputs$column
density_tail_input_order <- unique(density_tail_inputs$input)
tail_definition_order <- unique(c(tail_input_order, density_tail_input_order))

tail_condition_quantile <- 0.25

tail_subset_key <- function(inputs) {
  if (length(inputs) == 0) return("")
  paste(sort(inputs), collapse = "|")
}

all_tail_subsets <- function(inputs) {
  purrr::map(0:length(inputs), ~ utils::combn(inputs, .x, simplify = FALSE)) %>%
    purrr::flatten()
}

filter_bad_tail <- function(data, definitions, inputs) {
  out <- data
  for (input_name in inputs) {
    input_def <- definitions %>% filter(input == input_name)
    if (nrow(input_def) != 1) stop("Missing bad-tail definition for ", input_name)
    input_col <- tail_input_columns[[input_name]]
    if (input_def$direction == "high") {
      out <- out %>% filter(.data[[input_col]] >= input_def$cutoff)
    } else {
      out <- out %>% filter(.data[[input_col]] <= input_def$cutoff)
    }
  }
  out
}

tail_bad_definitions <- purrr::map_dfr(unique(tail_data$gas), function(gas_value) {
  gas_data <- tail_data %>%
    filter(gas == gas_value) %>%
    mutate(log_npc_norm = log10(npc_norm))

  purrr::map_dfr(tail_definition_order, function(input_name) {
    input_col <- tail_input_columns[[input_name]]
    input_values <- gas_data[[input_col]]
    lower_cutoff <- quantile(input_values, tail_condition_quantile,
                             na.rm = TRUE, type = QUANTILE_TYPE)
    upper_cutoff <- quantile(input_values, 1 - tail_condition_quantile,
                             na.rm = TRUE, type = QUANTILE_TYPE)
    lower_median <- median(gas_data$log_npc_norm[input_values <= lower_cutoff],
                           na.rm = TRUE)
    upper_median <- median(gas_data$log_npc_norm[input_values >= upper_cutoff],
                           na.rm = TRUE)
    bad_direction <- ifelse(upper_median >= lower_median, "high", "low")

    tibble(gas = gas_value,
           input = input_name,
           label = tail_input_labels[[input_name]],
           direction = bad_direction,
           cutoff = ifelse(bad_direction == "high", upper_cutoff, lower_cutoff))
  })
})

tail_subsets <- all_tail_subsets(tail_input_order)

tail_subset_stats <- purrr::map_dfr(unique(tail_data$gas), function(gas_value) {
  gas_data <- tail_data %>%
    filter(gas == gas_value) %>%
    mutate(log_npc_norm = log10(npc_norm))
  gas_definitions <- tail_bad_definitions %>%
    filter(gas == gas_value)

  purrr::map_dfr(tail_subsets, function(input_subset) {
    subset_data <- filter_bad_tail(gas_data, gas_definitions, input_subset)

    tibble(gas = gas_value,
           subset = tail_subset_key(input_subset),
           subset_size = length(input_subset),
           n = nrow(subset_data),
           median_log = ifelse(nrow(subset_data) == 0, NA_real_,
                               median(subset_data$log_npc_norm, na.rm = TRUE)),
           p95_log = ifelse(nrow(subset_data) == 0, NA_real_,
                            quantile(subset_data$log_npc_norm, 0.95,
                                     na.rm = TRUE, type = QUANTILE_TYPE)))
  })
})

tail_shapley_scores <- purrr::map_dfr(unique(tail_data$gas), function(gas_value) {
  gas_stats <- tail_subset_stats %>%
    filter(gas == gas_value)
  n_inputs <- length(tail_input_order)

  purrr::map_dfr(tail_input_order, function(input_name) {
    other_inputs <- setdiff(tail_input_order, input_name)

    all_tail_subsets(other_inputs) %>%
      purrr::map_dfr(function(input_subset) {
        subset_size <- length(input_subset)
        subset_weight <- factorial(subset_size) *
          factorial(n_inputs - subset_size - 1) / factorial(n_inputs)
        baseline <- gas_stats %>%
          filter(subset == tail_subset_key(input_subset))
        augmented <- gas_stats %>%
          filter(subset == tail_subset_key(c(input_subset, input_name)))

        tibble(weight = subset_weight,
               marginal_median_log = augmented$median_log - baseline$median_log,
               marginal_p95_log = augmented$p95_log - baseline$p95_log)
      }) %>%
      summarise(gas = gas_value,
                input = input_name,
                median_score = sum(weight * marginal_median_log, na.rm = TRUE),
                p95_score = sum(weight * marginal_p95_log, na.rm = TRUE),
                .groups = "drop")
  })
})

tail_origin_colors <- c(A = "#D81B60",
                        B = "#1E88E5",
                        C = "#FFC107",
                        D = "#004D40")

tail_top_inputs <- density_tail_inputs %>%
  select(gas, input) %>%
  left_join(tail_bad_definitions, by = c("gas", "input")) %>%
  group_by(gas) %>%
  mutate(origin = LETTERS[row_number()],
         origin_color = unname(tail_origin_colors[origin])) %>%
  ungroup()

blend_tail_colors <- function(origins) {
  if (length(origins) == 1) return(unname(tail_origin_colors[origins]))
  rgb_values <- grDevices::col2rgb(tail_origin_colors[origins]) / 255
  mixed <- rowSums(rgb_values)
  mixed <- mixed / max(mixed) * 0.82
  grDevices::rgb(mixed[1], mixed[2], mixed[3])
}

contrast_text_color <- function(colors) {
  rgb_values <- grDevices::col2rgb(colors) / 255
  luminance <- 0.299 * rgb_values[1, ] + 0.587 * rgb_values[2, ] + 0.114 * rgb_values[3, ]
  ifelse(luminance < 0.48, "white", "black")
}

tail_selected_subsets <- purrr::map_dfr(unique(tail_data$gas), function(gas_value) {
  gas_data <- tail_data %>%
    filter(gas == gas_value) %>%
    mutate(log_npc_norm = log10(npc_norm))
  gas_definitions <- tail_bad_definitions %>%
    filter(gas == gas_value)
  gas_top_definitions <- tail_top_inputs %>%
    filter(gas == gas_value) %>%
    arrange(origin)

  candidate_stats <- all_tail_subsets(gas_top_definitions$input) %>%
    purrr::keep(~ length(.x) %in% c(2, 3)) %>%
    purrr::map_dfr(function(input_subset) {
      subset_data <- filter_bad_tail(gas_data, gas_definitions, input_subset)
      origin_subset <- gas_top_definitions$origin[match(input_subset, gas_top_definitions$input)]

      tibble(gas = gas_value,
             input = tail_subset_key(input_subset),
             origin = tail_subset_key(origin_subset),
             subset_size = length(input_subset),
             median_log = median(subset_data$log_npc_norm, na.rm = TRUE),
             p95_log = quantile(subset_data$log_npc_norm, 0.95,
                                na.rm = TRUE, type = QUANTILE_TYPE),
             n = nrow(subset_data))
    })

  purrr::map_dfr(c(2, 3), function(subset_size_value) {
    size_stats <- candidate_stats %>%
      filter(subset_size == subset_size_value)
    median_pick <- size_stats %>%
      arrange(desc(median_log), desc(p95_log), desc(n), input) %>%
      slice_head(n = 1) %>%
      mutate(selection_reason = "median")
    p95_pick <- size_stats %>%
      arrange(desc(p95_log), desc(median_log), desc(n), input) %>%
      slice_head(n = 1) %>%
      mutate(selection_reason = "p95")

    bind_rows(median_pick, p95_pick) %>%
      group_by(gas, input, origin, subset_size, median_log, p95_log, n) %>%
      summarise(selection_reason = paste(selection_reason, collapse = "+"),
                .groups = "drop")
  })
})

tail_condition_label <- function(definitions, inputs) {
  definitions %>%
    filter(input %in% inputs) %>%
    mutate(input = factor(input, levels = inputs)) %>%
    arrange(input) %>%
    transmute(label = paste(label, direction, sep = " ")) %>%
    pull(label) %>%
    paste(collapse = " + ")
}

tail_density_overlay_data <- purrr::map_dfr(unique(tail_data$gas), function(gas_value) {
  gas_data <- tail_data %>%
    filter(gas == gas_value) %>%
    mutate(log_npc_norm = log10(npc_norm),
           gas_label = ifelse(gas_value == "ch4", "CH4", "CO2"))
  gas_definitions <- tail_bad_definitions %>%
    filter(gas == gas_value)
  gas_top_definitions <- tail_top_inputs %>%
    filter(gas == gas_value) %>%
    arrange(origin)
  gas_top_inputs <- gas_top_definitions$input
  gas_selected_subsets <- tail_selected_subsets %>%
    filter(gas == gas_value) %>%
    pull(input)

  full_density_data <- gas_data %>%
    transmute(gas, gas_label, log_npc_norm,
              condition = "Full distribution",
              condition_type = "full",
              input = NA_character_,
              origin = NA_character_,
              subset_size = 0L,
              tail_count = "Full")

  subset_density_data <- all_tail_subsets(gas_top_inputs) %>%
    purrr::keep(~ length(.x) > 0) %>%
    purrr::keep(~ length(.x) %in% c(1, 4) ||
                  tail_subset_key(.x) %in% gas_selected_subsets) %>%
    purrr::map_dfr(function(input_subset) {
      subset_size <- length(input_subset)
      origin_subset <- gas_top_definitions$origin[match(input_subset, gas_top_definitions$input)]
      origin_subset <- sort(origin_subset)
      filter_bad_tail(gas_data, gas_definitions, input_subset) %>%
        transmute(gas, gas_label, log_npc_norm,
                  condition = paste(origin_subset, collapse = "+"),
                  condition_type = "tail_subset",
                  input = tail_subset_key(input_subset),
                  origin = tail_subset_key(origin_subset),
                  subset_size = subset_size,
                  tail_count = paste(subset_size, ifelse(subset_size == 1, "tail", "tails")),
                  definition = tail_condition_label(gas_definitions, input_subset))
    })

  bind_rows(full_density_data, subset_density_data)
})

tail_density_conditions <- tail_density_overlay_data %>%
  filter(condition_type != "full") %>%
  distinct(gas, gas_label, condition, origin, subset_size, tail_count) %>%
  mutate(plot_color = ifelse(subset_size == 4,
                             "black",
                             purrr::map_chr(strsplit(origin, "\\|"), blend_tail_colors))) %>%
  group_by(gas, condition) %>%
  slice_head(n = 1) %>%
  ungroup()
tail_density_linetypes <- c("1 tail" = "solid",
                            "2 tails" = "dashed",
                            "3 tails" = "dotdash",
                            "4 tails" = "longdash")

tail_density_counts <- tail_density_overlay_data %>%
  filter(is.finite(log_npc_norm)) %>%
  count(gas, gas_label, condition, condition_type, subset_size,
        tail_count, origin, name = "sample_n")
tail_density_full_counts <- tail_density_counts %>%
  filter(condition_type == "full") %>%
  select(gas, full_n = sample_n)

compute_tail_density_curve <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < 2) return(tibble(x = numeric(), density = numeric()))
  density_values <- stats::density(values, adjust = 1, n = 512)
  tibble(x = density_values$x, density = density_values$y)
}

tail_density_curves <- tail_density_overlay_data %>%
  filter(is.finite(log_npc_norm)) %>%
  group_by(gas, gas_label, condition, condition_type,
           subset_size, tail_count, origin) %>%
  group_modify(~ compute_tail_density_curve(.x$log_npc_norm)) %>%
  ungroup() %>%
  left_join(tail_density_counts,
            by = c("gas", "gas_label", "condition", "condition_type",
                   "subset_size", "tail_count", "origin")) %>%
  left_join(tail_density_full_counts, by = "gas") %>%
  group_by(gas, condition) %>%
  mutate(sample_share = sample_n / full_n,
         height_weight = sample_share,
         weighted_density = density / max(density, na.rm = TRUE) * height_weight) %>%
  ungroup()

tail_combo_key <- function(origins) {
  paste(sort(origins), collapse = "+")
}

make_tail_density_gas_plot <- function(gas_value, gas_title) {
  gas_conditions <- tail_density_conditions %>%
    filter(gas == gas_value)
  gas_palette <- setNames(gas_conditions$plot_color, gas_conditions$condition)
  gas_curves <- tail_density_curves %>%
    filter(gas == gas_value)

  gas_combo_nodes <- gas_conditions %>%
    filter(subset_size %in% 1:3) %>%
    distinct(condition, origin, subset_size, tail_count, plot_color) %>%
    mutate(origin_list = strsplit(origin, "\\|"),
           x = subset_size,
           text_color = contrast_text_color(plot_color),
           node_label = gsub("\\+", "", condition)) %>%
    arrange(subset_size, condition)
  gas_combo_column_n <- gas_combo_nodes %>%
    count(subset_size, name = "column_n")
  gas_combo_max_column_n <- max(gas_combo_column_n$column_n)
  gas_combo_nodes <- gas_combo_nodes %>%
    left_join(gas_combo_column_n, by = "subset_size") %>%
    group_by(subset_size) %>%
    mutate(y = rev(seq_len(n())) + (gas_combo_max_column_n - column_n) / 2) %>%
    ungroup() %>%
    select(-column_n)
  gas_combo_lookup <- gas_combo_nodes %>%
    select(condition, target_x = x, target_y = y)

  gas_pair_edges <- gas_combo_nodes %>%
    filter(subset_size == 2) %>%
    select(x, y, origin_list) %>%
    tidyr::unnest_longer(origin_list, values_to = "target") %>%
    left_join(gas_combo_lookup, by = c("target" = "condition")) %>%
    transmute(x = x - 0.22, y, xend = target_x + 0.22, yend = target_y)

  gas_selected_pairs <- gas_combo_nodes %>%
    filter(subset_size == 2) %>%
    pull(condition)
  gas_triple_edges <- gas_combo_nodes %>%
    filter(subset_size == 3) %>%
    select(x, y, origin_list) %>%
    mutate(target = purrr::map(origin_list, function(origins) {
      selected_pairs <- purrr::map_chr(utils::combn(origins, 2, simplify = FALSE),
                                       tail_combo_key)
      selected_pairs <- intersect(selected_pairs, gas_selected_pairs)
      if (length(selected_pairs) == 0) origins else selected_pairs
    })) %>%
    select(-origin_list) %>%
    tidyr::unnest_longer(target) %>%
    left_join(gas_combo_lookup, by = c("target" = "condition")) %>%
    transmute(x = x - 0.22, y, xend = target_x + 0.22, yend = target_y)

  gas_combo_edges <- bind_rows(gas_pair_edges, gas_triple_edges) %>%
    filter(!is.na(xend), !is.na(yend))
  gas_combo_headers <- tibble(x = 1:3,
                              y = gas_combo_max_column_n + 1,
                              label = c("1 tail", "2 tails", "3 tails"),
                              tail_count = label)
  gas_combo_ymax <- gas_combo_max_column_n + 1.3

  gas_combo_legend_plot <- ggplot() +
    geom_segment(data = gas_combo_edges,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey72",
                 linewidth = 0.25,
                 arrow = grid::arrow(length = grid::unit(0.045, "inches"),
                                     type = "closed")) +
    geom_segment(data = gas_combo_headers,
                 aes(x = x - 0.22, y = y - 0.3, xend = x + 0.22, yend = y - 0.3,
                     linetype = tail_count),
                 color = "grey20",
                 linewidth = 0.6) +
    geom_tile(data = gas_combo_nodes,
              aes(x = x, y = y, fill = plot_color, linetype = tail_count),
              width = 0.36,
              height = 0.36,
              color = "grey15",
              linewidth = 0.35) +
    geom_text(data = gas_combo_nodes,
              aes(x = x, y = y, label = node_label, color = text_color),
              size = 2.3,
              fontface = "bold") +
    geom_text(data = gas_combo_headers,
              aes(x = x, y = y, label = label),
              size = 2.8,
              fontface = "bold") +
    scale_fill_identity() +
    scale_color_identity() +
    scale_linetype_manual(values = tail_density_linetypes, guide = "none") +
    coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0.55, gas_combo_ymax), clip = "off") +
    labs(title = paste(gas_title, "combination key")) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5,
                                    margin = margin(b = 2)))

  gas_origin_definition_lines <- tail_top_inputs %>%
    filter(gas == gas_value) %>%
    mutate(definition = paste0(origin, ": ", label, " ", direction)) %>%
    arrange(origin) %>%
    mutate(row_id = row_number(),
           text_color = contrast_text_color(origin_color))

  gas_origin_key_plot <- ggplot(gas_origin_definition_lines,
                                aes(y = -row_id)) +
    geom_segment(data = tibble(y = 0),
                 aes(x = 0, xend = 0.28, y = y, yend = y),
                 inherit.aes = FALSE,
                 color = "black",
                 linetype = tail_density_linetypes[["4 tails"]],
                 linewidth = 0.8) +
    geom_text(data = tibble(y = 0, label = "A+B+C+D: all four tails"),
              aes(x = 0.38, y = y, label = label),
              inherit.aes = FALSE,
              hjust = 0,
              size = 2.35) +
    geom_tile(aes(x = 0, fill = origin_color),
              width = 0.22,
              height = 0.5,
              color = "grey15",
              linewidth = 0.25) +
    geom_text(aes(x = 0, label = origin, color = text_color),
              size = 2.15,
              fontface = "bold") +
    geom_text(aes(x = 0.34, label = definition),
              hjust = 0,
              size = 2.2) +
    scale_fill_identity() +
    scale_color_identity() +
    coord_cartesian(xlim = c(-0.12, 3.9), ylim = c(-4.8, 0.7), clip = "off") +
    labs(title = "Origin definitions") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 9,
                                    margin = margin(b = 2)))

  gas_custom_density_legend <- gas_combo_legend_plot / gas_origin_key_plot +
    patchwork::plot_layout(heights = c(2.1, 1.15))

  gas_density_panel <- ggplot() +
    geom_area(data = gas_curves %>% filter(condition_type == "full"),
              aes(x = x, y = weighted_density),
              color = "grey55",
              fill = "grey82",
              alpha = 0.55,
              linewidth = 0.45) +
    geom_line(data = gas_curves %>%
                filter(condition_type != "full", subset_size < 4),
              aes(x = x, y = weighted_density, color = condition, linetype = tail_count),
              linewidth = 0.8) +
    geom_line(data = gas_curves %>%
                filter(condition_type != "full", subset_size == 4),
              aes(x = x, y = weighted_density, color = condition, linetype = tail_count),
              linewidth = 1.25) +
    scale_x_continuous(labels = pow10_labels) +
    scale_color_manual(values = gas_palette, guide = "none") +
    scale_linetype_manual(values = tail_density_linetypes, guide = "none") +
    xlab("Normalized present cost") +
    ylab("Sample-size weighted density height") +
    labs(title = gas_title) +
    ggpubr::theme_pubr() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  gas_density_panel + gas_custom_density_legend +
    patchwork::plot_layout(widths = c(4.6, 1.35))
}

tail_bad_density_ch4 <- make_tail_density_gas_plot("ch4", "CH4")
tail_bad_density_co2 <- make_tail_density_gas_plot("co2", "CO2")
tail_bad_density_plot <- tail_bad_density_ch4 / tail_bad_density_co2 +
  patchwork::plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(
    title = "Cost distributions under gas-specific bad input tails",
    subtitle = "Grey: full distribution. Colored curves: singles plus worst 2- and 3-tail combinations by median and/or p95; heights are weighted by sample share."
  )

tail_pair_index <- tidyr::crossing(row_input = tail_input_order,
                                   col_input = tail_input_order) %>%
  mutate(row_id = match(row_input, tail_input_order),
         col_id = match(col_input, tail_input_order)) %>%
  filter(row_id != col_id)

tail_offdiag <- purrr::pmap_dfr(tail_pair_index,
                                function(row_input, col_input, row_id, col_id) {
  gas_value <- ifelse(row_id > col_id, "ch4", "co2")

  tail_data %>%
    filter(gas == gas_value) %>%
    transmute(bad_tail,
              x_bin = bin_tail_input(.data[[tail_input_columns[[col_input]]]]),
              y_bin = bin_tail_input(.data[[tail_input_columns[[row_input]]]])) %>%
    group_by(x_bin, y_bin) %>%
    summarise(tail_share = mean(bad_tail, na.rm = TRUE),
              n = n(),
              .groups = "drop") %>%
    mutate(row_var = factor(tail_input_labels[[row_input]], levels = tail_input_labels),
           col_var = factor(tail_input_labels[[col_input]], levels = tail_input_labels),
           gas_label = ifelse(gas_value == "ch4", "CH4", "CO2"))
})

tail_share_limit <- max(0.6, ceiling(max(tail_offdiag$tail_share, na.rm = TRUE) * 10) / 10)
tail_density_data <- tail_data %>%
  filter(gas == "ch4")

tail_cell_theme <- theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, face = "bold", angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2))

make_tail_heatmap_cell <- function(row_input, col_input, row_id, col_id) {
  gas_label <- ifelse(row_id > col_id, "CH4", "CO2")
  gas_fill <- ifelse(gas_label == "CH4", "#F8766D", "#00BFC4")
  cell_data <- tail_offdiag %>%
    filter(row_var == tail_input_labels[[row_input]],
           col_var == tail_input_labels[[col_input]])

  ggplot(cell_data, aes(x = x_bin, y = y_bin, alpha = tail_share)) +
    geom_tile(fill = gas_fill, color = "white", linewidth = 0.15) +
    geom_text(aes(label = ifelse(tail_share >= 0.2,
                                 scales::percent(tail_share, accuracy = 1), "")),
              color = "white", size = 1.7, show.legend = FALSE) +
    scale_x_discrete(labels = strip_bin_prefix) +
    scale_y_discrete(labels = strip_bin_prefix) +
    scale_alpha_continuous(range = c(0.08, 1),
                           limits = c(0, tail_share_limit),
                           labels = scales::percent,
                           name = "Bad-tail share") +
    guides(alpha = guide_legend(override.aes = list(fill = "grey35",
                                                    color = "white",
                                                    linewidth = 0.15))) +
    labs(title = ifelse(row_id == 1, tail_input_labels[[col_input]], ""),
         y = ifelse(col_id == 1, tail_input_labels[[row_input]], "")) +
    tail_cell_theme +
    theme(axis.text.x = if (row_id == length(tail_input_order)) {
            element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
          } else {
            element_blank()
          },
          axis.ticks.x = if (row_id == length(tail_input_order)) element_line() else element_blank(),
          axis.text.y = if (col_id == 1) element_text(size = 5) else element_blank(),
          axis.ticks.y = if (col_id == 1) element_line() else element_blank())
}

make_tail_density_cell <- function(input_name, cell_id) {
  input_values <- tail_density_data[[tail_input_columns[[input_name]]]]

  ggplot(tibble(value = input_values), aes(x = value)) +
    geom_density(fill = "grey75", color = "grey25", linewidth = 0.35, alpha = 0.9,
                 adjust = 1) +
    labs(title = ifelse(cell_id == 1, tail_input_labels[[input_name]], ""),
         y = ifelse(cell_id == 1, tail_input_labels[[input_name]], "")) +
    tail_cell_theme +
    theme(axis.text.x = if (cell_id == length(tail_input_order)) {
            element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
          } else {
            element_blank()
          },
          axis.ticks.x = if (cell_id == length(tail_input_order)) element_line() else element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

tail_matrix_cells <- purrr::map(seq_along(tail_input_order), function(row_id) {
  purrr::map(seq_along(tail_input_order), function(col_id) {
    row_input <- tail_input_order[[row_id]]
    col_input <- tail_input_order[[col_id]]

    if (row_id == col_id) {
      make_tail_density_cell(row_input, row_id)
    } else {
      make_tail_heatmap_cell(row_input, col_input, row_id, col_id)
    }
  })
}) %>%
  purrr::flatten()

tail_contrast_matrix <- patchwork::wrap_plots(tail_matrix_cells,
                                              ncol = length(tail_input_order),
                                              guides = "collect") +
  plot_annotation(title = "Bad-tail concentration by input combination",
                  subtitle = "Lower triangle: CH4 bad-tail share. Upper triangle: CO2 bad-tail share. Diagonal: marginal input distributions.") &
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9))

save_figure("fig_2_tail_contrast_matrix.png",tail_contrast_matrix,width=12,height=11,dpi=300)
save_figure("fig_2_tail_density_overlay.png",tail_bad_density_plot,width=16,height=9,dpi=300)
