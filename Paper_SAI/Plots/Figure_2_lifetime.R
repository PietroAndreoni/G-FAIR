# =============================================================================
# Figure_2_lifetime.R
#
# The top panel of Figure 2 -- the density of the normalized net present cost --
# drawn for the single-period pulse and for the emitting asset with an uncertain
# lifetime, one facet per gas. Construction follows the `density_plot` block of
# Figure_2.R: same log10 x axis with pow10 tick labels, same Silverman-style
# bandwidth, same rug, same coord_cartesian window, same ggpubr::theme_pubr.
#
#   Rscript Paper_SAI/Plots/Figure_2_lifetime.R [pulse_stream] [var_stream]
#   defaults: Results_16072026  MC_infra
#
# TWO DEVIATIONS FROM Figure_2.R, both required by the question being asked:
#
#  1. NORMALIZATION. Figure_2.R divides each gas by ITS OWN median. Doing that per
#     stream here would re-centre both curves on 1 and erase the very shift we are
#     trying to see. Instead both streams are divided by the PULSE median of that
#     gas, so the pulse curve is centred at 1 by construction and the asset curve
#     shows its true displacement as well as its shape.
#
#  2. COLOUR. Figure_2.R uses colour for the gas; here the gas is the facet and
#     colour carries the stream, so a single visual channel does one job. Blue and
#     orange are the pair that survives the common colour-vision deficiencies.
#
# PAIRING. Streams are matched 1-for-1 on the full sampled design (everything
# except infra_life and the model outputs), and, as in Figure_2.R, only draws that
# solved for BOTH gases are kept. Every draw behind the blue curve is therefore
# also behind the orange one, so any difference is the lifetime and nothing else.
# =============================================================================

library(tidyverse)

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
resolve_stream <- function(name) {
  cand <- c(file.path(MC_WORK_PARENT, name), file.path(RUN_RESULTS_PARENT, name), name)
  for (p in cand) if (file.exists(file.path(p, "npc_output.csv"))) return(p)
  stop("no npc_output.csv found for stream '", name, "'. Looked in:\n  ",
       paste(cand, collapse = "\n  "), call. = FALSE)
}
f_pulse <- resolve_stream(if (length(.args) >= 1) .args[1] else "Results_16072026")
f_var   <- resolve_stream(if (length(.args) >= 2) .args[2] else "MC_infra")
cat("pulse stream: ", f_pulse, "\nvar stream:   ", f_var, "\n", sep = "")

read_npc <- function(folder) bind_rows(lapply(
  file.path(folder, list.files(path = folder, pattern = "npc_output")), read.csv)) %>% unique()

pulse <- read_npc(f_pulse)
var   <- read_npc(f_var)
# The archive predates the infrastructure axis: every run in it is a pulse, so
# L = 1 is exact rather than an assumption.
if (!"infra_life" %in% names(pulse)) pulse$infra_life <- 1L
if (!"infra_life" %in% names(var))
  stop("the variable stream has no infra_life column -- wrong folder?", call. = FALSE)

out_cols <- c("ID", "infra_life", "dirnpv", "srmpnpv", "ozpnpv", "masknpv", "damnpv",
              "pulse_size", "pulse_size_annual", "npc_srm", "npc_srm_cap")
key <- setdiff(intersect(names(pulse), names(var)), out_cols)
cat("pairing on", length(key), "design columns\n")
stopifnot(!any(duplicated(pulse[key])), !any(duplicated(var[key])))

paired <- inner_join(pulse[key], var[key], by = key)
keep <- function(x, lab) x %>% inner_join(paired, by = key) %>% mutate(stream = lab)
both <- bind_rows(keep(pulse, "Pulse (L = 1)"),
                  keep(var,   "Emitting asset (L ~ logU 1-100)")) %>%
  filter(npc_srm > 0)

# Same both-gases requirement as Figure_2.R, applied to the PAIR: a draw is kept
# only if all four cells (2 gases x 2 streams) survived, so the facets are built
# from one common set of draws.
draw_key <- setdiff(key, "gas")
complete <- both %>% count(across(all_of(draw_key))) %>% filter(n == 4) %>% select(-n)
both <- both %>% inner_join(complete, by = draw_key)
cat("complete draws (both gases, both streams):", nrow(complete), "\n")
if (nrow(complete) == 0) stop("no draws matched across the two streams.", call. = FALSE)

# Deviation 1: normalize BOTH streams by the pulse median of that gas.
pulse_med <- both %>% filter(stream == "Pulse (L = 1)") %>%
  group_by(gas) %>% summarise(med = median(npc_srm, na.rm = TRUE), .groups = "drop")
damnorm <- both %>% left_join(pulse_med, by = "gas") %>%
  mutate(npc_norm = npc_srm / med,
         stream = factor(stream, levels = c("Pulse (L = 1)",
                                            "Emitting asset (L ~ logU 1-100)")),
         gas = factor(gas, levels = c("co2", "ch4"),
                      labels = c("CO2", "CH4"))) %>%
  filter(is.finite(npc_norm), npc_norm > 0)

# Bandwidth exactly as Figure_2.R derives it, from the pooled standardized costs.
npc_std <- damnorm %>% group_by(gas) %>%
  mutate(s = (npc_srm - median(npc_srm, na.rm = TRUE)) / mad(npc_srm, na.rm = TRUE)) %>%
  pull(s)
adjust_opt <- log(FIG_KDE_BW_CONST * sd(npc_std) * length(npc_std)^(FIG_KDE_BW_POW))
pow10_labels <- scales::label_math(10^.x)

density_plot <- ggplot(damnorm) +
  geom_density(aes(x = log10(npc_norm), color = stream),
               adjust = adjust_opt, linewidth = 1.5) +
  geom_point(aes(x = log10(npc_norm), color = stream,
                 y = -(as.numeric(stream) - 1) / 100),
             shape = 108) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3, color = "grey40") +
  scale_color_manual(name = "", values = c("#08306B", "#D95F02")) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  scale_x_continuous(labels = pow10_labels) +
  facet_wrap(gas ~ .) +
  xlab("Present cost, relative to the median pulse") + ylab("Density") +
  ggpubr::theme_pubr(legend = "top")

save_figure("fig_2_lifetime.png", density_plot, width = 10, height = 5, dpi = 300)
cat("saved fig_2_lifetime.png\n")

# --- the numbers behind the curves -------------------------------------------
cat("\nnpc_srm (USD/t), matched draws only:\n")
print(as.data.frame(damnorm %>% group_by(gas, stream) %>%
  summarise(n = n(), med = median(npc_srm), p25 = quantile(npc_srm, .25),
            p75 = quantile(npc_srm, .75), p95 = quantile(npc_srm, .95),
            mad = mad(npc_srm), .groups = "drop")))

cat("\npaired ratio asset/pulse, per draw:\n")
wide <- damnorm %>% select(all_of(key), stream, npc_srm) %>%
  pivot_wider(names_from = stream, values_from = npc_srm) %>%
  rename(np = `Pulse (L = 1)`, nv = `Emitting asset (L ~ logU 1-100)`) %>%
  filter(is.finite(np), is.finite(nv), np > 0, nv > 0) %>%
  mutate(ratio = nv / np)
print(as.data.frame(wide %>% group_by(gas) %>%
  summarise(n = n(), med = median(ratio), p25 = quantile(ratio, .25),
            p75 = quantile(ratio, .75), pct_above_1 = round(100 * mean(ratio > 1)),
            .groups = "drop")))

cat("\nKS test on log10(npc_norm), pulse vs asset:\n")
for (g in levels(damnorm$gas)) {
  a <- damnorm %>% filter(gas == g, stream == "Pulse (L = 1)") %>% pull(npc_norm)
  b <- damnorm %>% filter(gas == g, stream != "Pulse (L = 1)") %>% pull(npc_norm)
  k <- suppressWarnings(stats::ks.test(log10(a), log10(b)))
  cat(sprintf("  %s: D = %.4f, p = %.3g\n", g, k$statistic, k$p.value))
}
