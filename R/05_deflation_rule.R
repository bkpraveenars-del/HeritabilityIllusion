# 05_deflation_rule.R -- test H2_single / (1 + rho) against the direct REML estimate
# of the one-environment H2, by trial balance (share of filled G x E cells).
source("R/00_setup.R")
A <- read.csv("results/agridat_summary.csv")
C <- read.csv("results/cimmyt_summary.csv"); C <- C[C$set == "All environments", ]
D <- readRDS("results/cimmyt_clean.rds"); P <- read.csv("results/cimmyt_per_env.csv")
D <- D[D$env %in% P$env, ]
fillC <- sapply(C$nursery, function(n) { x <- D[D$nursery == n, ]; mean(table(droplevels(x$gen), x$env) > 0) })

X <- data.frame(source = c(rep("agridat", nrow(A)), rep("CIMMYT", nrow(C))),
                id = c(A$ds, C$nursery), H2_single = c(A$H2_single_mean, C$H2_single_mean),
                rho = c(A$rho, C$rho), H2_true1 = c(A$H2_true_1env, C$H2_true1), fill = c(A$fill, fillC))
X <- X[is.finite(X$rho), ]                     # linder.wheat has sigma2_G = 0
X$H2_deflated <- X$H2_single / (1 + X$rho); X$error <- X$H2_deflated - X$H2_true1
X$balanced <- X$fill >= 0.8
write.csv(X, "results/deflation_rule.csv", row.names = FALSE)
for (b in c(TRUE, FALSE)) {
  e <- X$error[X$balanced == b]
  cat(if (b) "Balanced  " else "Unbalanced", " n =", length(e), " RMSE =", round(sqrt(mean(e^2)), 3),
      " bias =", round(mean(e), 3), " max|error| =", round(max(abs(e)), 3), "\n")
}
