# 03_cimmyt_analysis.R -- single-environment and multi-environment H2, leave-one-
# environment-out (LOEO) selection, and robustness subsets for the CIMMYT ESWYT data.
# Run with argument "std" to standardise yield within environment:
#   Rscript R/03_cimmyt_analysis.R std
source("R/00_setup.R")
args <- commandArgs(trailingOnly = TRUE); STD <- length(args) > 0 && args[1] == "std"
tag  <- if (STD) "cimmyt_std" else "cimmyt"
D <- readRDS("results/cimmyt_clean.rds")
SA <- c("INDIA", "PAKISTAN", "NEPAL", "BANGLADESH")

# ---- single-environment fits (alpha-lattice: replicate fixed, incomplete block random)
per <- list()
for (e in unique(D$env)) {
  d <- droplevels(D[D$env == e, ])
  if (length(unique(d$rep)) < 2 || nlevels(d$gen) < 30) next
  d$rep <- factor(d$rep); d$rb <- factor(paste(d$rep, d$blk))
  f <- if (nlevels(d$rb) > nlevels(d$rep)) y ~ rep + (1|rb) + (1|gen) else y ~ rep + (1|gen)
  m <- tryCatch(suppressMessages(lmer(f, d, control = CTL)), error = function(z) NULL)
  if (is.null(m)) next
  v <- vc(m); r <- nrow(d) / nlevels(d$gen)
  per[[e]] <- data.frame(nursery = d$nursery[1], env = e, country = d$country[1],
                         sigma2_g = v[["gen"]], sigma2_e = v[["Residual"]], r = r,
                         H2 = v[["gen"]] / (v[["gen"]] + v[["Residual"]] / r),
                         mean = mean(d$y), cv = 100 * sqrt(v[["Residual"]]) / mean(d$y))
}
P <- do.call(rbind, per); rownames(P) <- NULL
D <- D[D$env %in% P$env, ]
if (STD) D$y <- ave(D$y, D$env, FUN = function(z) (z - mean(z)) / sd(z))

# ---- multi-environment fit
met <- function(d) {
  d$env <- factor(d$env); d$er <- factor(paste(d$env, d$rep)); d$erb <- factor(paste(d$env, d$rep, d$blk))
  m <- suppressMessages(lmer(y ~ (1|env) + (1|er) + (1|erb) + (1|gen) + (1|gen:env), d, control = CTL))
  v <- vc(m); rb <- nrow(d) / (nlevels(droplevels(d$gen)) * nlevels(d$env)); ne <- nlevels(d$env)
  c(sigma2_G = v[["gen"]], sigma2_GE = v[["gen:env"]], sigma2_e = v[["Residual"]], r = rb,
    H2_true1 = v[["gen"]] / (v[["gen"]] + v[["gen:env"]] + v[["Residual"]] / rb),
    H2_MET   = v[["gen"]] / (v[["gen"]] + v[["gen:env"]] / ne + v[["Residual"]] / (rb * ne)))
}

# ---- LOEO: select top 5 of 49 in one environment, measure advantage elsewhere
loeo <- function(d, Pn) {
  M <- tapply(d$y, list(d$gen, d$env), mean); out <- list()
  for (e in colnames(M)) {
    xi <- M[, e]; xo <- rowMeans(M[, colnames(M) != e, drop = FALSE], na.rm = TRUE)
    ok <- is.finite(xi) & is.finite(xo); xi <- xi[ok]; xo <- xo[ok]
    sel <- order(xi, decreasing = TRUE)[1:5]
    S <- mean(xi[sel]) - mean(xi); H <- Pn$H2[Pn$env == e]
    out[[e]] <- data.frame(env = e, pred = H * S, real = mean(xo[sel]) - mean(xo), r_cross = cor(xi, xo))
  }
  do.call(rbind, out)
}

run <- function(set, sub) {
  res <- list(); L <- list()
  for (n in unique(sub$nursery)) {
    d <- sub[sub$nursery == n, ]; Pn <- P[P$env %in% unique(d$env), ]
    if (nrow(Pn) < 4) next
    mm <- met(d); lo <- loeo(d, Pn); L[[n]] <- cbind(nursery = n, lo)
    res[[n]] <- data.frame(set = set, nursery = n, nenv = nrow(Pn),
      H2_single_mean = mean(Pn$H2), H2_single_median = median(Pn$H2),
      H2_true1 = mm[["H2_true1"]], H2_MET = mm[["H2_MET"]], rho = mm[["sigma2_GE"]] / mm[["sigma2_G"]],
      realised_over_predicted = sum(lo$real) / sum(lo$pred))
  }
  list(R = do.call(rbind, res), L = do.call(rbind, L))
}

A <- run("All environments", D)
S <- run("South Asia only",  D[D$country %in% SA, ])
C <- run("CV <= 20% only",   D[D$env %in% P$env[P$cv <= 20], ])
R <- rbind(A$R, S$R, C$R); R$inflation <- R$H2_single_mean - R$H2_true1; rownames(R) <- NULL

write.csv(R,   paste0("results/", tag, "_summary.csv"), row.names = FALSE)
write.csv(P,   paste0("results/", tag, "_per_env.csv"), row.names = FALSE)
write.csv(A$L, paste0("results/", tag, "_loeo.csv"),    row.names = FALSE)

print(format(R, digits = 2), row.names = FALSE)
L <- A$L
cat("\nEnvironments analysed:", nrow(P), " countries:", length(unique(P$country)), "\n")
cat("Median single-env H2:", round(median(P$H2), 2),
    " | India:", round(median(P$H2[P$country == "INDIA"]), 2), "(n =", sum(P$country == "INDIA"), ")\n")
cat("LOEO realised/predicted gain:", round(sum(L$real) / sum(L$pred), 3),
    " | environments with no gain:", sum(L$real <= 0),
    " | median cross-env r:", round(median(L$r_cross, na.rm = TRUE), 2),
    " | median sqrt(H2):", round(median(sqrt(P$H2)), 2), "\n")
