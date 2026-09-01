# =============================================================================
# Figure_3_lifetime.R
#
# Figure 3's abatement-cost densities, drawn twice: once for the classic
# single-period pulse and once for the emitting asset with an uncertain lifetime.
# Everything else -- the Harmsen MACC curve, the two cost definitions, the year
# facets, the axes and the cost unit -- is taken from Figure_3.R unchanged, so the
# only thing that differs between the two sets of curves is infra_life.
#
#   Rscript Paper_SAI/Plots/Figure_3_lifetime.R [pulse_stream] [var_stream]
#   defaults: Results_base_16072026  MC_base_var2
#
# PAIRING. The two streams are matched 1-for-1 on the full sampled design
# (everything except infra_life and the model outputs) and only matched draws are
# plotted, so the two densities are built from the same underlying draws and any
# difference between them is the lifetime and nothing else. An unmatched draw on
# either side is dropped rather than left to shift one density on its own.
#
# The pulse arm is READ FROM THE ARCHIVE, not re-run: Results_base_16072026 already
# holds those runs, and its design is a superset of the variable-lifetime one.
#
# Both streams divide by pulse_size = LIFETIME tonnes, which is what makes the x
# axis comparable to a MACC in $/ton. The per-capacity figure (npc_srm_cap) is a
# different question and is deliberately not plotted here.
# =============================================================================

require(witchtools)
pkgs <- c('tidyverse', 'stringr', 'arrow', 'data.table')
res <- lapply(pkgs, require_package)

.find_paper_root <- function() {
  starts <- c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
              unlist(lapply(sys.frames(), function(f) f$ofile)), getwd())
  for (s in starts[nzchar(starts)]) {
    d <- if (dir.exists(s)) s else dirname(s)
    repeat {
      if (file.exists(file.path(d, "all_parameters.R"))) return(normalizePath(d, "/", FALSE))
      if (file.exists(file.path(d, "Paper_SAI", "all_parameters.R")))
        return(normalizePath(file.path(d, "Paper_SAI"), "/", FALSE))
      if (identical(dirname(d), d)) break
      d <- dirname(d)
    }
  }
  stop("Cannot locate all_parameters.R (Paper_SAI control file).", call. = FALSE)
}
source(file.path(.find_paper_root(), "all_parameters.R"))

.args <- commandArgs(TRUE)

# Same resolution rule as compare_lifetime_streams.R / Figure_lifetime.R: a stream
# is named, not pathed, and may sit under Sampling/ or under Results/.
resolve_stream <- function(name) {
  cand <- c(file.path(MC_WORK_PARENT, name), file.path(RUN_RESULTS_PARENT, name), name)
  for (p in cand) if (file.exists(file.path(p, "npc_output.csv"))) return(p)
  stop("no npc_output.csv found for stream '", name, "'. Looked in:\n  ",
       paste(cand, collapse = "\n  "), call. = FALSE)
}
f_pulse <- resolve_stream(if (length(.args) >= 1) .args[1] else "Results_base_16072026")
f_var   <- resolve_stream(if (length(.args) >= 2) .args[2] else "MC_base_var2")
cat("pulse stream: ", f_pulse, "\nvar stream:   ", f_var, "\n", sep = "")

read_npc <- function(folder) {
  f <- list.files(path = folder, pattern = "npc_output")
  bind_rows(lapply(file.path(folder, f), read.csv)) %>% filter(gas == "ch4") %>% unique()
}
pulse <- read_npc(f_pulse)
var   <- read_npc(f_var)

# The archive predates the infrastructure axis and carries no infra_life column;
# every run in it is a single-period pulse, so L = 1 is exact, not an assumption.
if (!"infra_life" %in% names(pulse)) pulse$infra_life <- 1L
if (!"infra_life" %in% names(var))
  stop("the variable stream has no infra_life column -- wrong folder?", call. = FALSE)

# --- 1-for-1 pairing ---------------------------------------------------------
# The design key is every sampled input the two streams share, minus infra_life
# (the one axis that is meant to differ) and minus the outputs.
out_cols <- c("ID", "infra_life", "dirnpv", "srmpnpv", "ozpnpv", "masknpv", "damnpv",
              "pulse_size", "pulse_size_annual", "npc_srm", "npc_srm_cap")
key <- setdiff(intersect(names(pulse), names(var)), out_cols)
cat("pairing on", length(key), "design columns\n")

stopifnot(!any(duplicated(pulse[key])), !any(duplicated(var[key])))
paired_ids <- inner_join(pulse[key], var[key], by = key)
cat(sprintf("pulse draws %d | var draws %d | MATCHED %d\n",
            nrow(pulse), nrow(var), nrow(paired_ids)))
if (nrow(paired_ids) == 0) stop("no draws matched across the two streams.", call. = FALSE)

keep <- function(x, lab) x %>% inner_join(paired_ids, by = key) %>% mutate(stream = lab)
data <- bind_rows(
  keep(pulse, "Pulse (L = 1)"),
  keep(var,   "Emitting asset (L ~ logU 1-100)")) %>%
  filter(npc_srm > 0) %>%
  mutate(year = pulse_time + 2020) %>%
  filter(year %in% FIG3_MACC_YEARS) %>%
  mutate(stream = factor(stream, levels = c("Pulse (L = 1)",
                                            "Emitting asset (L ~ logU 1-100)")))
cat("plotted rows by stream x year:\n")
print(as.data.frame(data %>% count(stream, year)))

# --- the MACC curve, exactly as in Figure_3.R --------------------------------
baseline <- read_parquet(HARMSEN_BASELINE_FIG3)
macc     <- read_parquet(HARMSEN_MACC_FIG3)
macc_by_gas_w <- macc %>%
  filter(unit == "emissions") %>%
  full_join(baseline %>% rename(base = value)) %>%
  group_by(year, e, cost) %>%
  summarize(value = sum(value, na.rm = TRUE), base = sum(base, na.rm = TRUE)) %>%
  ungroup() %>% mutate(miu = (base - value) / base, e = tolower(e)) %>%
  mutate(cost = cost * AR4_GWP100[e] * C_PER_CO2 * USD_DEFLATOR_2010_2020) %>%
  select(year, e, cost, miu)

density_cost_unit <- DENSITY_COST_UNIT

x_axis <- scale_x_continuous(
  labels = ~paste(., . / CH4_GWP100, sep = "\n"),
  name = expression(atop("Abatement cost ($/ton" * CH[4] * ")",
                         "Abatement cost ($/ton" * CO[2] * "eq)")))
macc_line <- geom_line(
  data = macc_by_gas_w %>% filter(e == "ch4" & year %in% FIG3_MACC_YEARS),
  aes(y = miu * 100, x = cost), color = "black", linewidth = 1, linetype = 2)

# --- (1) Figure 3 itself, drawn once per stream ------------------------------
# Encoding is left exactly as published -- light blue = cost excluding ozone,
# dark blue = total cost -- and the stream becomes a facet row. Overlaying the two
# streams inside one panel was tried first and rejected: it needs a second visual
# channel for four curves, and solid-vs-dashed at these linewidths is not readable.
fig <- ggplot() +
  macc_line +
  geom_density(data = data,
               aes(x = (masknpv + damnpv + dirnpv + srmpnpv) / pulse_size,
                   y = after_stat(density * 100 * density_cost_unit)),
               adjust = 2, linetype = 1, linewidth = 1, color = "#6BAED6") +
  geom_density(data = data,
               aes(x = (masknpv + damnpv + dirnpv + srmpnpv + ozpnpv) / pulse_size,
                   y = after_stat(density * 100 * density_cost_unit)),
               adjust = 2, linetype = 1, linewidth = 1.5, color = "#08306B") +
  geom_hline(yintercept = 0) +
  geom_point(data = data, aes(x = npc_srm, y = 0), shape = 108, color = "#08306B") +
  facet_grid(stream ~ year) +
  theme_classic() +
  ylab("Emission reductions (% of baseline)\nDensity (% per $100/ton)") +
  theme(legend.position = "none", strip.text.y = element_text(size = 8)) +
  coord_cartesian(xlim = FIG3_XLIM, ylim = FIG3_YLIM) + x_axis

save_figure("fig_3_lifetime.png", fig, width = 12, height = 8, dpi = 300)
cat("saved fig_3_lifetime.png  (Figure 3 replicated per stream)\n")

# --- (2) the direct contrast -------------------------------------------------
# One curve per stream, total cost only, so the two are on the same axes with a
# single visual channel doing the work. Blue/orange is the pair that survives the
# common colour-vision deficiencies, and the rug carries the same colour.
fig_overlay <- ggplot() +
  macc_line +
  geom_density(data = data,
               aes(x = (masknpv + damnpv + dirnpv + srmpnpv + ozpnpv) / pulse_size,
                   y = after_stat(density * 100 * density_cost_unit),
                   color = stream),
               adjust = 2, linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 0) +
  geom_point(data = data, aes(x = npc_srm, y = 0, color = stream),
             shape = 108, show.legend = FALSE) +
  scale_color_manual(name = "", values = c("#08306B", "#D95F02")) +
  facet_wrap(year ~ .) +
  theme_classic() +
  ylab("Emission reductions (% of baseline)\nDensity (% per $100/ton)") +
  theme(legend.position = "top") +
  coord_cartesian(xlim = FIG3_XLIM, ylim = FIG3_YLIM) + x_axis

save_figure("fig_3_lifetime_overlay.png", fig_overlay, width = 12, height = 6, dpi = 300)
cat("saved fig_3_lifetime_overlay.png  (total cost, both streams)\n")

# --- the numbers behind the curves -------------------------------------------
cat("\nnpc_srm (USD/tCH4), matched draws only:\n")
print(as.data.frame(data %>% group_by(year, stream) %>%
  summarise(n = n(), med = median(npc_srm), p25 = quantile(npc_srm, .25),
            p75 = quantile(npc_srm, .75), p95 = quantile(npc_srm, .95),
            .groups = "drop")))

# The share of the CH4 MACC that SAI undercuts is what the figure is read for:
# the fraction of draws whose cost falls below a given abatement cost.
cat("\nshare of draws below a given abatement cost ($/tCH4):\n")
for (thr in c(1000, 2500, 5000)) {
  s <- data %>% group_by(year, stream) %>%
    summarise(pct = round(100 * mean(npc_srm < thr)), .groups = "drop") %>%
    mutate(threshold = thr)
  print(as.data.frame(s))
}
