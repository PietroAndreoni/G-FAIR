# =============================================================================
# Figure_lifetime.R
#
# Normalized net present cost of SAI substitution against the lifetime of the
# emitting infrastructure, one GAM smooth per gas. Built on the model of the
# sensitivity panels in Figure_2.R (`dr`, `alpha`, `theta`, `term`): same
# npc_norm definition, same log10 axes with pow10 tick labels, same GAM
# smoother, same coord_cartesian y window, same ggpubr::theme_pubr.
#
#   Rscript Paper_SAI/Plots/Figure_lifetime.R [var_folder] [pulse_folder]
#   defaults: MC_var  MC_pulse
#
# Two figures are written:
#
#   fig_lifetime_npc.png    the requested Figure_2-style sensitivity panel. Like
#                           every panel in Figure_2.R it is UNPAIRED: each draw
#                           has its own ecs/alpha/delta/rcp/cooling, and that
#                           between-draw variance is left in. It answers "across
#                           the ensemble, how does cost move with lifetime".
#
#   fig_lifetime_paired.png the same axis, but y is each draw's cost divided by
#                           the cost of THE SAME draw run as a single-period
#                           pulse. The pulse campaign is identical in every other
#                           sampled dimension (same Sobol net, same seed), so this
#                           differences out all between-draw variance and isolates
#                           the lifetime effect. Only produced if the pulse
#                           campaign is present.
#
# The x axis only carries information when the campaign was drawn with a
# non-degenerate lifetime support (Generate_montecarlo.R --lifetime LO:HI); a
# pulse campaign pins infra_life = 1 and the panel collapses to a point.
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

# Same resolution rule as compare_lifetime_streams.R: a stream is named, not
# pathed, and may live either under Sampling/ (a working campaign) or under
# Results/ (an archived one, e.g. Results_16072026). An explicit path also works.
resolve_stream <- function(name) {
  cand <- c(file.path(MC_WORK_PARENT, name), file.path(RUN_RESULTS_PARENT, name), name)
  for (p in cand) if (file.exists(file.path(p, "npc_output.csv"))) return(p)
  stop("no npc_output.csv found for stream '", name, "'. Looked in:\n  ",
       paste(cand, collapse = "\n  "), call. = FALSE)
}

work_var   <- resolve_stream(if (length(.args) >= 1) .args[1] else "MC_var")
work_pulse <- if (length(.args) >= 2) resolve_stream(.args[2]) else
              file.path(MC_WORK_PARENT, "MC_pulse")
cat("var stream:   ", work_var, "\npulse stream: ", work_pulse, "\n", sep = "")

read_npc <- function(folder) {
  f <- list.files(path = folder, pattern = "npc_output")
  if (!length(f)) stop("no npc_output in ", folder, call. = FALSE)
  bind_rows(lapply(file.path(folder, f), read.csv))
}

damnpv <- read_npc(work_var)
if (!"infra_life" %in% names(damnpv))
  stop("npc_output in ", work_var, " has no infra_life column.", call. = FALSE)

pow10_labels <- scales::label_math(10^.x)

# Same normalization as Figure_2.R: each gas divided by its own median, so the two
# gases share one y axis despite differing by orders of magnitude in level.
damnorm <- damnpv %>% unique() %>% select(-ID) %>%
  mutate_if(is.integer, as.numeric) %>%
  group_by(gas) %>%
  mutate(npc_norm = npc_srm / median(npc_srm, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(is.finite(npc_norm), npc_norm > 0, infra_life > 0)

# k = 7 as in the `dr` panel. These campaigns run in the hundreds of draws per gas
# rather than the 8192 those panels were tuned on, so the basis is kept small to
# avoid reading noise as structure.
lifetime <- ggplot(damnorm) +
  geom_smooth(aes(x = log10(infra_life), y = log10(npc_norm), color = gas),
              method = "gam", formula = y ~ s(x, bs = "cs", k = 7)) +
  scale_x_continuous(labels = pow10_labels) +
  scale_y_continuous(labels = pow10_labels) +
  xlab("Emitting infrastructure lifetime [yr]") + ylab("Normalized cost") +
  coord_cartesian(ylim = c(-0.5, 1)) +
  ggpubr::theme_pubr(legend = "right")

save_figure("fig_lifetime_npc.png", lifetime, width = 7, height = 5, dpi = 300)
cat("saved fig_lifetime_npc.png  (unpaired, Figure_2 style)\n")

# --- paired view -------------------------------------------------------------
if (file.exists(file.path(work_pulse, "npc_output.csv"))) {
  pulse <- read_npc(work_pulse) %>% select(ID, gas, npc_pulse = npc_srm)
  paired <- damnpv %>% select(ID, gas, infra_life, npc_var = npc_srm) %>%
    inner_join(pulse, by = c("ID", "gas")) %>%
    filter(is.finite(npc_var), is.finite(npc_pulse), npc_pulse > 0, npc_var > 0) %>%
    mutate(ratio = npc_var / npc_pulse)

  cat("paired (ID,gas) cells:", nrow(paired), "\n")

  lifetime_paired <- ggplot(paired) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
    geom_point(aes(x = log10(infra_life), y = log10(ratio), color = gas),
               alpha = 0.25, size = 1, show.legend = FALSE) +
    geom_smooth(aes(x = log10(infra_life), y = log10(ratio), color = gas),
                method = "gam", formula = y ~ s(x, bs = "cs", k = 7)) +
    scale_x_continuous(labels = pow10_labels) +
    scale_y_continuous(labels = pow10_labels) +
    xlab("Emitting infrastructure lifetime [yr]") +
    ylab("Cost relative to the same draw as a pulse") +
    ggpubr::theme_pubr(legend = "right")

  save_figure("fig_lifetime_paired.png", lifetime_paired, width = 7, height = 5, dpi = 300)
  cat("saved fig_lifetime_paired.png  (paired against the pulse campaign)\n")
} else {
  cat("no pulse campaign at ", work_pulse, " -- paired panel skipped\n", sep = "")
}
