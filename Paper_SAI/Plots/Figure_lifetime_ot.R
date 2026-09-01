# =============================================================================
# Figure_lifetime_ot.R
#
# Optimal-transport global sensitivity indices on log10(npc_norm), with the
# emitting-infrastructure lifetime included as an input. Same construction as the
# gsaot block of Figure_2.R (ot_indices_1d at M = GSA_M, irrelevance_threshold as
# the noise floor, one panel per gas, viridis fill by input category), so the
# result is directly comparable to the published importance figure -- the
# question it answers is where infra_life ranks among the existing inputs.
#
#   Rscript Paper_SAI/Plots/Figure_lifetime_ot.R [work_folder]   (default: MC_var)
#
# Two deviations from Figure_2.R, both deliberate:
#   * the y limit is taken from the data instead of the hard-coded c(0, 0.175).
#     That constant was set for the 8192-draw campaign; clipping here would
#     silently drop bars rather than show them.
#   * the irrelevance threshold is printed as well as drawn, because with ~1e3
#     draws rather than 8192 the noise floor sits high enough to matter when
#     reading the ranking.
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
library(gsaot); require(patchwork)

.args <- commandArgs(TRUE)
work  <- file.path(MC_WORK_PARENT, if (length(.args) >= 1) .args[1] else "MC_var")

damnpv <- bind_rows(lapply(
  file.path(work, list.files(path = work, pattern = "npc_output")), read.csv))
if (!"infra_life" %in% names(damnpv))
  stop("npc_output in ", work, " has no infra_life column.", call. = FALSE)

damnorm <- damnpv %>% unique() %>%
  mutate_if(is.integer, as.numeric) %>%
  group_by(gas) %>%
  mutate(npc_norm = npc_srm / median(npc_srm, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(is.finite(npc_norm), npc_norm > 0)

# Everything that is an OUTPUT or a bookkeeping column is dropped; what remains is
# the sampled input set. infra_life stays -- that is the point of this figure.
# `term` is dropped as in Figure_2.R: term_delta is a deterministic function of
# pulse and prob, both of which are already inputs, so keeping it double-counts.
drop_cols <- c("gas","ID","term","pulse_size","pulse_size_annual",
               "npc_srm","npc_srm_cap","npc_norm","npc_std",
               "dirnpv","srmpnpv","ozpnpv","masknpv","damnpv")

ot_for <- function(g) {
  dat <- damnorm %>% filter(gas == g) %>% select(-any_of(drop_cols))
  y   <- log10(damnorm %>% filter(gas == g) %>% pull(npc_norm))
  dat <- dat %>% select(where(~ n_distinct(.x) > 1))    # OT needs a varying input
  list(idx = ot_indices_1d(dat, y, M = GSA_M),
       lb  = irrelevance_threshold(y, M = GSA_M, solver = "1d"),
       n   = length(y))
}
ch4 <- ot_for("ch4"); co2 <- ot_for("co2")

# --- the CO2 - CH4 contrast --------------------------------------------------
# Same construction, but y is the per-draw gap between the two gases rather than
# either gas's level. Because the panels above use y = log10(npc_norm), the gap is
#   log10(npc_norm_co2) - log10(npc_norm_ch4) = log10(npc_co2 / npc_ch4) + const,
# so "difference of the normalized costs" and "ratio of the raw costs" are the same
# response up to an additive constant -- and OT indices are shift-invariant, so the
# panel does not depend on which of the two readings is meant.
#
# gas is no longer an input here: every sampled parameter is drawn per DRAW, not per
# gas, so the two rows of a draw share them. That is asserted rather than assumed.
ot_diff <- function() {
  wide <- damnorm %>% select(ID, gas, npc_norm) %>%
    pivot_wider(names_from = gas, values_from = npc_norm, names_prefix = "npc_") %>%
    filter(is.finite(npc_co2), is.finite(npc_ch4), npc_co2 > 0, npc_ch4 > 0)

  keep <- function(g) damnorm %>% filter(gas == g) %>%
    select(!any_of(setdiff(drop_cols, "ID"))) %>% arrange(ID)
  a <- keep("co2"); b <- keep("ch4")
  shared <- intersect(a$ID, b$ID)
  stopifnot(identical(as.data.frame(a[a$ID %in% shared, ]),
                      as.data.frame(b[b$ID %in% shared, ])))

  d <- inner_join(a, wide, by = "ID")
  y <- log10(d$npc_co2) - log10(d$npc_ch4)
  dat <- d %>% select(!any_of(c("ID", "npc_co2", "npc_ch4"))) %>%
    select(where(~ n_distinct(.x) > 1))
  list(idx = ot_indices_1d(dat, y, M = GSA_M),
       lb  = irrelevance_threshold(y, M = GSA_M, solver = "1d"),
       n   = length(y))
}
dif <- ot_diff()
cat(sprintf("draws: ch4 %d, co2 %d, paired %d | inputs: %d\n",
            ch4$n, co2$n, dif$n, length(ch4$idx$indices)))
cat(sprintf("irrelevance threshold: ch4 %.4f, co2 %.4f, diff %.4f\n",
            ch4$lb$indices, co2$lb$indices, dif$lb$indices))

input_categories <- c(ecs="Climate", tcr="Climate",
                      rcp="Socio-economic", pulse_time="Socio-economic",
                      dg="Socio-economic", infra_life="Socio-economic",
                      geo_start="SAI", geo_end="SAI", cool_rate="SAI", prob="SAI",
                      forctoTg="SAI", TgtoUSD="SAI",
                      alpha="Impacts", theta="Impacts",
                      mortality_ozone="Impacts", mortality_srm="Impacts",
                      delta="Normative", vsl="Normative", vsl_eta="Normative")

input_axis_labels <- c(ecs="plain(ECS)", tcr="plain(TCR)", rcp="plain(RCP)",
                       pulse_time="plain(yr)[pulse]", infra_life="plain(L)[asset]",
                       geo_start="plain(SAI)[start]", geo_end="plain(SAI)[end]",
                       cool_rate="plain(SAI)[rate]", alpha="alpha", delta="delta",
                       theta="theta", prob="plain(p)", dg="plain(g)",
                       mortality_ozone="mu[O[3]]", vsl="plain(VSL)",
                       vsl_eta="eta[plain(VSL)]", mortality_srm="mu[plain(SAI)]",
                       forctoTg="frac(plain(Tg) ~ plain(yr)^-1, plain(W) ~ plain(m)^-2)",
                       TgtoUSD="frac(plain(USD), plain(Tg) ~ plain(yr)^-1)")

ymax <- 1.15 * max(c(ch4$idx$indices, co2$idx$indices, dif$idx$indices), na.rm = TRUE)

panel <- function(o, lab) {
  d <- tibble(input = names(o$idx$indices), original = o$idx$indices)
  ggplot() +
    geom_bar(data = d, aes(x = reorder(input, -original), y = original,
                           fill = input_categories[input]),
             stat = "identity", color = "black") +
    geom_hline(yintercept = o$lb$indices) +
    ylab(lab) + xlab("") +
    scale_fill_viridis_d(name = "") + ylim(c(0, ymax)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2),
                     labels = function(x) parse(text = unname(input_axis_labels[x]))) +
    ggpubr::theme_pubr()
}

importances <- (panel(co2, expression("importance [" * log[10](CO[2]) * "]")) +
                panel(ch4, expression("importance [" * log[10](CH[4]) * "]")) +
                panel(dif, expression("importance [" * log[10](CO[2]/CH[4]) * "]"))) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

save_figure("fig_lifetime_ot.png", importances, width = 18, height = 6, dpi = 300)
cat("saved fig_lifetime_ot.png\n")

rank_of <- function(o, g) {
  i <- sort(o$idx$indices, decreasing = TRUE)
  cat(sprintf("%s: infra_life index %.4f -> rank %d of %d (threshold %.4f, %s)\n",
              g, o$idx$indices[["infra_life"]], which(names(i) == "infra_life"),
              length(i), o$lb$indices,
              ifelse(o$idx$indices[["infra_life"]] > o$lb$indices,
                     "ABOVE noise floor", "below noise floor")))
  cat("   top 5:", paste(sprintf("%s=%.3f", names(i)[1:5], i[1:5]), collapse = "  "), "\n")
}
rank_of(co2, "co2"); rank_of(ch4, "ch4"); rank_of(dif, "co2-ch4")
