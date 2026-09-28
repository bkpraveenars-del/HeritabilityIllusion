# 04_cimmyt_bootstrap.R -- cluster bootstrap over environments for the inflation
# (mean single-env H2 minus one-environment H2) in each CIMMYT nursery.
# Note: resampled duplicate environments are identical copies, which reduces apparent
# G x E; bootstrap intervals for rho are therefore biased low and are not reported.
source("R/00_setup.R")
set.seed(2026)
D <- readRDS("results/cimmyt_clean.rds"); P <- read.csv("results/cimmyt_per_env.csv")
D <- D[D$env %in% P$env, ]

fitmet <- function(d) {
  d$env <- factor(d$env); d$er <- factor(paste(d$env, d$rep)); d$erb <- factor(paste(d$env, d$rep, d$blk))
  m <- suppressMessages(lmer(y ~ (1|env) + (1|er) + (1|erb) + (1|gen) + (1|gen:env), d, control = CTL))
  list(v = vc(m), rb = nrow(d) / (nlevels(droplevels(d$gen)) * nlevels(d$env)))
}
h1 <- function(f) f$v[["gen"]] / (f$v[["gen"]] + f$v[["gen:env"]] + f$v[["Residual"]] / f$rb)

out <- list()
for (n in unique(D$nursery)) {
  d <- D[D$nursery == n, ]; Pn <- P[P$nursery == n, ]; f <- fitmet(d)
  envs <- unique(d$env)
  bt <- replicate(BOOT_CIMMYT, {
    s  <- sample(envs, replace = TRUE)
    db <- do.call(rbind, lapply(seq_along(s), function(k) { z <- d[d$env == s[k], ]; z$env <- paste0("b", k); z }))
    fb <- tryCatch(fitmet(db), error = function(e) NULL)
    if (is.null(fb)) NA else mean(Pn$H2[match(s, Pn$env)]) - h1(fb) })
  out[[n]] <- data.frame(nursery = n, nenv = length(envs), H2_single = mean(Pn$H2), H2_true1 = h1(f),
                         rho = f$v[["gen:env"]] / f$v[["gen"]], inflation = mean(Pn$H2) - h1(f),
                         infl_lo = quantile(bt, .025, na.rm = TRUE), infl_hi = quantile(bt, .975, na.rm = TRUE),
                         nboot = sum(is.finite(bt)))
  message(n, " done")
}
O <- do.call(rbind, out); rownames(O) <- NULL
write.csv(O, "results/cimmyt_bootstrap.csv", row.names = FALSE)
print(format(O, digits = 3), row.names = FALSE)
