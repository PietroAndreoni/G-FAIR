library(tidyverse)
library(patchwork)

# Locate Paper_SAI whether run from the project root, Plots/, or RStudio.
find_paper_root <- function() {
  starts <- c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
              unlist(lapply(sys.frames(), function(frame) frame$ofile)),
              getwd())
  for (start in starts[nzchar(starts)]) {
    directory <- if (dir.exists(start)) start else dirname(start)
    repeat {
      if (file.exists(file.path(directory, "all_parameters.R"))) {
        return(normalizePath(directory, "/", FALSE))
      }
      if (file.exists(file.path(directory, "Paper_SAI", "all_parameters.R"))) {
        return(normalizePath(file.path(directory, "Paper_SAI"), "/", FALSE))
      }
      if (identical(dirname(directory), directory)) break
      directory <- dirname(directory)
    }
  }
  stop("Cannot locate Paper_SAI/all_parameters.R.", call. = FALSE)
}

source(file.path(find_paper_root(), "all_parameters.R"))

# Read and pair the CH4/CO2 Monte Carlo observations used by the matrix.
npc_files <- list.files(RESULTS_FOLDER_MAIN, pattern = "npc_output", full.names = TRUE)
if (length(npc_files) == 0) {
  stop("No npc_output files found in ", RESULTS_FOLDER_MAIN, call. = FALSE)
}

damnpv <- bind_rows(lapply(npc_files, read.csv))
pairing_columns <- names(damnpv)[1:21]

tail_data <- damnpv %>%
  filter(npc_srm > 0) %>%
  group_by(across(all_of(setdiff(pairing_columns, "gas")))) %>%
  filter(n() == 2, !anyDuplicated(gas)) %>%
  ungroup() %>%
  distinct() %>%
  group_by(gas) %>%
  mutate(npc_norm = npc_srm / median(npc_srm, na.rm = TRUE)) %>%
  ungroup() %>%
  transmute(gas,
            npc_norm,
            discount = delta * 100,
            damage_alpha = alpha * 100,
            sai_angle = theta,
            term_prob = prob)

# The four inputs selected by FIG2_MATRIX_INPUTS in all_parameters.R.
input_catalog <- tibble::tribble(
  ~param,  ~input,          ~label,       ~distribution,
  "delta", "discount",     "Discount",  "uniform",
  "alpha", "damage_alpha", "Damages",   "other",
  "theta", "sai_angle",    "SAI angle", "other",
  "prob",  "term_prob",    "Term. prob.", "loguniform"
)

if (!setequal(FIG2_MATRIX_INPUTS, input_catalog$param)) {
  stop("Additional_figure_4.R expects FIG2_MATRIX_INPUTS to contain delta, alpha, theta, and prob.",
       call. = FALSE)
}

inputs <- input_catalog %>%
  arrange(match(param, FIG2_MATRIX_INPUTS))
input_order <- inputs$input
input_labels <- setNames(inputs$label, inputs$input)
input_distributions <- setNames(inputs$distribution, inputs$input)

prefix_bins <- function(values) {
  sprintf("%02d: %s", seq_along(values), values)
}

strip_bin_prefix <- function(values) {
  sub("^\\d+:\\s*", "", values)
}

bin_input <- function(values, max_bins = 5) {
  unique_values <- sort(unique(values[!is.na(values)]))
  if (length(unique_values) <= 8) {
    formatter <- if (max(abs(unique_values), na.rm = TRUE) >= 100) "%.0f" else "%.1f"
    raw_labels <- sprintf(formatter, unique_values)
    label_map <- setNames(prefix_bins(raw_labels), raw_labels)
    return(unname(label_map[sprintf(formatter, values)]))
  }

  cut_values <- as.character(ggplot2::cut_number(
    values,
    n = min(max_bins, length(unique_values)),
    dig.lab = 4
  ))
  raw_labels <- unique(cut_values[order(values)])
  unname(setNames(prefix_bins(raw_labels), raw_labels)[cut_values])
}

# Apply the same trimming used by the Figure 2 sensitivity panels.
input_trim <- setNames(lapply(input_order, function(input_name) {
  if (identical(input_distributions[[input_name]], "other")) {
    unname(quantile(tail_data[[input_name]],
                    c(FIG2_PANEL_QLO, FIG2_PANEL_QHI),
                    na.rm = TRUE))
  } else {
    NULL
  }
}), input_order)

trim_input <- function(input_name, values) {
  limits <- input_trim[[input_name]]
  if (!is.null(limits)) {
    values[values < limits[[1]] | values > limits[[2]]] <- NA
  }
  values
}

pair_index <- tidyr::crossing(row_input = input_order,
                              col_input = input_order) %>%
  mutate(row_id = match(row_input, input_order),
         col_id = match(col_input, input_order)) %>%
  filter(row_id != col_id)

mean_offdiag <- purrr::pmap_dfr(
  pair_index,
  function(row_input, col_input, row_id, col_id) {
    gas_value <- ifelse(row_id > col_id, "ch4", "co2")
    tail_data %>%
      filter(gas == gas_value) %>%
      transmute(
        metric = npc_norm,
        x_bin = bin_input(trim_input(col_input, .data[[col_input]])),
        y_bin = bin_input(trim_input(row_input, .data[[row_input]]))
      ) %>%
      filter(!is.na(x_bin), !is.na(y_bin)) %>%
      group_by(x_bin, y_bin) %>%
      summarise(mean_damnorm = mean(metric, na.rm = TRUE), .groups = "drop") %>%
      mutate(row_var = factor(input_labels[[row_input]], levels = input_labels),
             col_var = factor(input_labels[[col_input]], levels = input_labels))
  }
)

cell_theme <- theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, face = "bold", angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2))

make_heatmap_cell <- function(row_input, col_input, row_id, col_id) {
  gas_fill <- ifelse(row_id > col_id, "#F8766D", "#00BFC4")
  cell_data <- mean_offdiag %>%
    filter(row_var == input_labels[[row_input]],
           col_var == input_labels[[col_input]])

  ggplot(cell_data, aes(x = x_bin, y = y_bin)) +
    geom_tile(aes(alpha = mean_damnorm),
              fill = gas_fill, color = "white", linewidth = 0.15) +
    geom_text(aes(label = ifelse(mean_damnorm > 1.5,
                                 scales::number(mean_damnorm, accuracy = 0.1),
                                 "")),
              color = "white", size = 1.7, show.legend = FALSE) +
    scale_x_discrete(labels = strip_bin_prefix) +
    scale_y_discrete(labels = strip_bin_prefix) +
    scale_alpha_continuous(
      range = c(0.08, 1),
      limits = mean_alpha_limits,
      oob = scales::squish,
      labels = function(value) ifelse(
        value >= mean_alpha_limits[[2]],
        paste0(">", scales::number(value, accuracy = 1)),
        scales::number(value, accuracy = 0.1)
      ),
      name = "Mean of normalized present cost"
    ) +
    guides(alpha = guide_legend(
      nrow = 1,
      override.aes = list(fill = "grey35", color = "white", linewidth = 0.15)
    )) +
    labs(title = ifelse(row_id == 1, input_labels[[col_input]], ""),
         y = ifelse(col_id == 1, input_labels[[row_input]], "")) +
    cell_theme +
    theme(
      axis.text.x = if (row_id == length(input_order)) {
        element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
      } else {
        element_blank()
      },
      axis.ticks.x = if (row_id == length(input_order)) element_line() else element_blank(),
      axis.text.y = if (col_id == 1) element_text(size = 5) else element_blank(),
      axis.ticks.y = if (col_id == 1) element_line() else element_blank()
    )
}

ch4_data <- tail_data %>% filter(gas == "ch4")
pow10_labels <- scales::label_math(10^.x)

make_density_cell <- function(input_name, cell_id) {
  input_values <- ch4_data[[input_name]]
  distribution <- input_distributions[[input_name]]

  marginal_layers <- if (identical(distribution, "uniform")) {
    list(geom_bar(
      data = tibble(value = factor(round(input_values, 6))),
      aes(x = value),
      fill = "grey75", color = "grey25", linewidth = 0.35, alpha = 0.9
    ))
  } else if (identical(distribution, "loguniform")) {
    list(
      geom_density(
        data = tibble(value = log10(input_values)), aes(x = value),
        fill = "grey75", color = "grey25", linewidth = 0.35, alpha = 0.9
      ),
      scale_x_continuous(labels = pow10_labels)
    )
  } else {
    trimmed_values <- trim_input(input_name, input_values)
    list(geom_density(
      data = tibble(value = trimmed_values[!is.na(trimmed_values)]), aes(x = value),
      fill = "grey75", color = "grey25", linewidth = 0.35, alpha = 0.9
    ))
  }

  ggplot() +
    marginal_layers +
    labs(title = ifelse(cell_id == 1, input_labels[[input_name]], ""),
         y = ifelse(cell_id == 1, input_labels[[input_name]], "")) +
    cell_theme +
    theme(
      axis.text.x = if (cell_id == length(input_order)) {
        element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
      } else {
        element_blank()
      },
      axis.ticks.x = if (cell_id == length(input_order)) element_line() else element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

matrix_cells <- purrr::map(seq_along(input_order), function(row_id) {
  purrr::map(seq_along(input_order), function(col_id) {
    row_input <- input_order[[row_id]]
    col_input <- input_order[[col_id]]
    if (row_id == col_id) {
      make_density_cell(row_input, row_id)
    } else {
      make_heatmap_cell(row_input, col_input, row_id, col_id)
    }
  })
}) %>%
  purrr::flatten()

additional_figure_4 <- patchwork::wrap_plots(
  matrix_cells,
  ncol = length(input_order),
  guides = "collect"
) &
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

save_figure(
  "fig_2_tail_mean_damnorm_matrix.png",
  additional_figure_4,
  width = 2 + 2 * length(input_order),
  height = 1 + 2 * length(input_order),
  dpi = 300
)
