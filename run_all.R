# run_all.R -- reproduce every result and figure in the manuscript.
# Usage (from this folder):  Rscript run_all.R
# Quick test with few bootstrap resamples:  HI_BOOT=5 Rscript run_all.R
# Full run time is roughly 1-2 hours on a laptop, mostly in the bootstraps.
steps <- c("R/01_agridat_analysis.R", "R/02_cimmyt_clean.R", "R/03_cimmyt_analysis.R",
           "R/03_cimmyt_analysis.R std", "R/04_cimmyt_bootstrap.R", "R/05_deflation_rule.R", "R/06_figures.R")
for (s in steps) {
  cat("\n==========", s, "==========\n")
  parts <- strsplit(s, " ")[[1]]
  st <- system2("Rscript", c(shQuote(parts[1]), parts[-1]))
  if (st != 0) stop("Step failed: ", s)
}
sink("results/sessionInfo.txt"); print(sessionInfo()); sink()
