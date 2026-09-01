# =============================================================================
# compare_lifetime_streams.R
#
# Compares two Monte Carlo campaigns that are identical draw for draw EXCEPT for
# the emitting-infrastructure lifetime:
#   pulse : infra_life == 1 for every draw (the classic single-period pulse)
#   var   : infra_life ~ round(LogUniform(1, 100))
# Same seed, same Sobol net, same 26 other columns -- so any difference in the
# outputs is attributable to the lifetime axis alone, paired by ID and gas.
#
#   Rscript Paper_SAI/Sampling/compare_lifetime_streams.R [pulse_folder] [var_folder]
#
# Reads the npc/scc outputs each Analyze_montecarlo.R run wrote into the two
# working folders under Sampling/, and prints paired statistics per gas.
# =============================================================================

# Locate the Paper_SAI folder (holds all_parameters.R) the same way every other
# script in the pipeline does, so this works under Rscript (--file), RStudio
# "Source" (sys.frame $ofile), `R -f`, and an interactive console whose working
# directory is at, under, or above the project. The earlier one-liner read
# commandArgs("--file=") directly and produced "NA/../all_parameters.R" whenever
# that argument was absent.
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
  stop("Cannot locate all_parameters.R (Paper_SAI control file); set the working ",
       "directory to the project root or the Paper_SAI folder.", call. = FALSE)
}
source(file.path(.find_paper_root(), "all_parameters.R"))
suppressMessages(library(data.table))

args <- commandArgs(TRUE)

# A stream is named by a FOLDER holding npc_output.csv / scc_output.csv /
# sccnosrm_output.csv. Accept a working folder under Sampling/, an archived
# analysis folder under Paper_SAI/Results/, or an explicit path -- so an archived
# campaign (e.g. Results_16072026) can serve as the pulse reference directly.
resolve_stream <- function(name) {
  cand <- c(file.path(MC_WORK_PARENT, name), file.path(RUN_RESULTS_PARENT, name), name)
  for (p in cand) if (file.exists(file.path(p, "npc_output.csv"))) return(p)
  stop("no npc_output.csv found for stream '", name, "'. Looked in:\n  ",
       paste(cand, collapse = "\n  "), call. = FALSE)
}

f_pulse <- resolve_stream(if (length(args) >= 1) args[1] else "MC_pulse")
f_var   <- resolve_stream(if (length(args) >= 2) args[2] else "MC_var")

read_stream <- function(folder, label) {
  rd <- function(f) {
    p <- file.path(folder, f)
    if (!file.exists(p)) stop("missing ", p, call. = FALSE)
    as.data.table(read.csv(p, stringsAsFactors = FALSE))
  }
  npc <- rd("npc_output.csv")
  scc <- rd("scc_output.csv")
  nos <- rd("sccnosrm_output.csv")

  # A campaign analysed before the lifetime axis existed carries neither
  # infra_life nor the *_cap columns. Such a campaign IS a single-period pulse,
  # so all four are reconstructed EXACTLY rather than approximated: at L = 1 the
  # lifetime-tonne and annual-capacity denominators are the same number.
  legacy <- !"infra_life" %in% names(npc)
  if (legacy) {
    message(folder, ": legacy analysis (no infra_life / *_cap columns); ",
            "treating it as infra_life = 1 and setting the capacity metrics ",
            "equal to the per-tonne ones, which is exact for a pulse.")
    npc[, infra_life := 1L]
    npc[, pulse_size_annual := pulse_size]
    npc[, npc_srm_cap := npc_srm]
    scc[, scc_cap := scc]
    nos[, scc_nosrm_cap := scc_nosrm]
  }

  keep_npc <- intersect(c("ID","gas","infra_life","pulse_size","pulse_size_annual",
                          "npc_srm","npc_srm_cap"), names(npc))
  d <- npc[, ..keep_npc]
  d <- merge(d, scc[, .(ID, gas, scc, scc_cap)], by = c("ID","gas"), all.x = TRUE)
  d <- merge(d, nos[, .(ID, gas, scc_nosrm, scc_nosrm_cap)], by = c("ID","gas"), all.x = TRUE)
  d <- unique(d)
  d[, stream := label][]
}

P <- read_stream(f_pulse, "pulse")
V <- read_stream(f_var,   "var")

cat("=====================================================================\n")
cat(" Emitting-infrastructure lifetime: pulse (L=1) vs variable (L~logU[1,100])\n")
cat("=====================================================================\n\n")
cat(sprintf("draws analysed   pulse: %d rows (%d IDs)   var: %d rows (%d IDs)\n",
            nrow(P), uniqueN(P$ID), nrow(V), uniqueN(V$ID)))

# ---- paired set: only IDs present in BOTH streams for BOTH gases -------------
key <- merge(unique(P[, .(ID, gas)]), unique(V[, .(ID, gas)]), by = c("ID","gas"))
P <- P[key, on = c("ID","gas")]; V <- V[key, on = c("ID","gas")]
cat(sprintf("paired (ID,gas) cells common to both streams: %d\n\n", nrow(key)))

setnames(V, setdiff(names(V), c("ID","gas")), paste0(setdiff(names(V), c("ID","gas")), "_v"))
M <- merge(P, V, by = c("ID","gas"))
stopifnot(nrow(M) == nrow(key))

# ---- sanity: the two streams really are the same experiment but for L --------
cat("--- structural checks ---\n")
chk <- function(lbl, ok) cat(sprintf("  %-62s %s\n", lbl, if (isTRUE(ok)) "PASS" else "FAIL"))
chk("pulse stream has infra_life == 1 everywhere", all(M$infra_life == 1))
chk("annual emission rate identical across streams (per gas)",
    isTRUE(all.equal(M$pulse_size_annual, M$pulse_size_annual_v)))
chk("lifetime tonnes = annual rate x infra_life (var stream)",
    isTRUE(all.equal(M$pulse_size_v, M$pulse_size_annual_v * as.numeric(M$infra_life_v))))
chk("L == 1 rows agree exactly between streams (scc)",
    { s <- M[infra_life_v == 1]; nrow(s) == 0 || isTRUE(all.equal(s$scc, s$scc_v)) })

cat("\n--- lifetime drawn in the variable stream ---\n")
print(M[, .(n = .N, min = min(infra_life_v), q25 = quantile(infra_life_v, .25),
            median = as.numeric(median(infra_life_v)), q75 = quantile(infra_life_v, .75),
            max = max(infra_life_v), mean = round(mean(infra_life_v), 1)), by = gas])

# ---- distributional comparison ----------------------------------------------
qs <- c(.05, .25, .5, .75, .95)
stat_block <- function(x) {
  x <- x[is.finite(x)]
  as.list(c(n = length(x), mean = mean(x), sd = sd(x),
            setNames(quantile(x, qs, names = FALSE), paste0("p", qs * 100))))
}

metrics <- list(
  "scc  (USD/tCO2e, per lifetime tonne)"        = c("scc", "scc_v"),
  "scc_cap  (USD per tonne/yr capacity)"        = c("scc_cap", "scc_cap_v"),
  "scc_nosrm  (USD/tCO2e, no SRM)"              = c("scc_nosrm", "scc_nosrm_v"),
  "npc_srm  (USD/tCO2e, net present cost)"      = c("npc_srm", "npc_srm_v"),
  "npc_srm_cap  (USD per tonne/yr capacity)"    = c("npc_srm_cap", "npc_srm_cap_v")
)

for (g in sort(unique(M$gas))) {
  cat("\n\n=====================  GAS =", toupper(g), " =====================\n")
  Mg <- M[gas == g]
  for (nm in names(metrics)) {
    a <- metrics[[nm]][1]; b <- metrics[[nm]][2]
    if (!all(c(a, b) %in% names(Mg))) next
    tab <- rbind(
      data.table(stream = "pulse (L=1)", as.data.table(stat_block(Mg[[a]]))),
      data.table(stream = "var (L~logU)", as.data.table(stat_block(Mg[[b]])))
    )
    cat("\n", nm, "\n", sep = "")
    print(tab[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 4) else z)], row.names = FALSE)

    # paired difference and ratio, on rows where both are finite and pulse != 0
    d <- Mg[is.finite(get(a)) & is.finite(get(b))]
    if (nrow(d) > 2) {
      ratio <- d[[b]] / d[[a]]
      ratio <- ratio[is.finite(ratio)]
      w <- suppressWarnings(wilcox.test(d[[b]], d[[a]], paired = TRUE))
      cat(sprintf("   paired: median ratio var/pulse = %.3f   [p25 %.3f, p75 %.3f]   Wilcoxon p = %.3g\n",
                  median(ratio), quantile(ratio, .25), quantile(ratio, .75), w$p.value))
    }
  }

  # does the gap track the lifetime?
  Mg2 <- Mg[is.finite(scc) & is.finite(scc_v) & scc != 0]
  if (nrow(Mg2) > 5) {
    cat(sprintf("\n   corr(log L, log(scc_var/scc_pulse)) = %.3f   (Spearman %.3f)\n",
                cor(log(Mg2$infra_life_v), log(pmax(Mg2$scc_v / Mg2$scc, 1e-12))),
                cor(log(Mg2$infra_life_v), log(pmax(Mg2$scc_v / Mg2$scc, 1e-12)), method = "spearman")))
    br <- Mg2[, .(n = .N,
                  scc_pulse = signif(median(scc), 4),
                  scc_var = signif(median(scc_v), 4),
                  ratio = signif(median(scc_v / scc), 3)),
              by = .(L_bin = cut(infra_life_v, c(0, 1, 3, 10, 32, 100),
                                 labels = c("1", "2-3", "4-10", "11-32", "33-100")))]
    cat("\n   median scc by lifetime bin:\n")
    print(br[order(L_bin)], row.names = FALSE)
  }
}
cat("\n")
