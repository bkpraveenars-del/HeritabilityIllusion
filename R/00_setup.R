# 00_setup.R -- install required packages and define shared helpers
# Heritability Illusion: single-environment heritability overstates selection response
# Author: Praveen Kumar (ORCID 0000-0002-7736-5986)

pkgs <- c("lme4", "agridat")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(lme4); library(agridat) })

# Every mixed model uses the bobyqa optimizer: the default optimizer converged to a
# spurious optimum in one fit during the original analysis.
CTL <- lmerControl(optimizer = "bobyqa")

# Named vector of variance components from an lmer fit
vc <- function(m) { v <- as.data.frame(VarCorr(m)); setNames(v$vcov, v$grp) }

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Number of bootstrap resamples (override with env var HI_BOOT for a quick test run)
BOOT_AGRIDAT <- as.integer(Sys.getenv("HI_BOOT", "200"))
BOOT_CIMMYT  <- as.integer(Sys.getenv("HI_BOOT", "60"))
