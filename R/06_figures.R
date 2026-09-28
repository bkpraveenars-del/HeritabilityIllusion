# 06_figures.R -- Figures 1-3 of the manuscript (PNG, 300 dpi)
source("R/00_setup.R")
A <- read.csv("results/agridat_summary.csv"); P <- read.csv("results/cimmyt_per_env.csv")
C <- read.csv("results/cimmyt_summary.csv"); C <- C[C$set == "All environments", ]
L <- read.csv("results/cimmyt_loeo.csv"); X <- read.csv("results/deflation_rule.csv")
col1 <- "#2b6cb0"; col2 <- "#c05621"; grey <- "grey55"
pal4 <- c("#2b6cb0", "#c05621", "#2f855a", "#6b46c1")
dev_png <- function(f, w, h) png(f, width = w, height = h, res = 300, type = if (capabilities("cairo")) "cairo" else "windows")

## Figure 1
dev_png("figures/Fig1.png", 2400, 1150)
layout(matrix(1:2, 1), widths = c(1.15, 1)); par(mar = c(4.2, 7.2, 2.2, 1), cex = 0.8)
A <- A[order(A$H2_single_mean), ]; n <- nrow(A)
plot(NA, xlim = c(0, 1), ylim = c(1, n + 3), yaxt = "n", xlab = expression("Broad-sense heritability (" * H^2 * ")"),
     ylab = "", main = "a  Nineteen public trial series (agridat)", adj = 0)
segments(A$H2_true_1env, 1:n, A$H2_single_mean, 1:n, col = grey)
points(A$H2_single_mean, 1:n, pch = 19, col = col2); points(A$H2_true_1env, 1:n, pch = 19, col = col1)
axis(2, at = 1:n, labels = A$ds, las = 1, cex.axis = 0.6, tick = FALSE, line = -0.6)
legend("topleft", c("Single-environment (mean)", "One environment, GxE removed"), pch = 19, col = c(col2, col1), bty = "n", cex = 0.8)
par(mar = c(4.2, 4.2, 2.2, 1))
hist(P$H2, breaks = seq(0, 1, 0.05), col = "grey85", border = "white", xlab = expression("Single-environment " * H^2),
     main = paste0("b  CIMMYT ESWYT, ", nrow(P), " environments"), adj = 0, xlim = c(0, 1), ylim = c(0, 32))
abline(v = C$H2_true1, col = col1, lwd = 2); abline(v = median(P$H2), col = col2, lwd = 2, lty = 2)
legend(0.63, 32, c("Median single-environment", "True one-environment (4 nurseries)"), col = c(col2, col1), lty = c(2, 1), lwd = 2, bty = "n", cex = 0.75)
dev.off()

## Figure 2
dev_png("figures/Fig2.png", 2400, 1100); par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1), cex = 0.8)
lim <- range(c(L$pred, L$real)); pc <- as.numeric(factor(L$nursery))
plot(L$pred, L$real, pch = 19, col = adjustcolor(pal4[pc], 0.6), xlim = lim, ylim = lim,
     xlab = "Predicted gain in selection environment (t/ha)", ylab = "Realised gain in other environments (t/ha)",
     main = "a  Leave-one-environment-out selection", adj = 0)
abline(0, 1, lty = 2, col = grey); abline(h = 0, col = "grey80")
legend("topleft", levels(factor(L$nursery)), pch = 19, col = pal4, bty = "n", cex = 0.7)
M <- merge(L, P[, c("env", "H2")])
plot(sqrt(M$H2), M$r_cross, pch = 19, col = adjustcolor(col1, 0.5), xlim = c(0, 1), ylim = c(-0.6, 1),
     xlab = expression("Implied accuracy  " * sqrt(H^2) * " (single environment)"),
     ylab = "Observed r with mean of other environments", main = "b  Promised vs observed accuracy", adj = 0)
abline(0, 1, lty = 2, col = grey); abline(h = 0, col = "grey80")
dev.off()

## Figure 3
dev_png("figures/Fig3.png", 1300, 1200); par(mar = c(4.2, 4.2, 1.5, 1), cex = 0.8)
plot(X$H2_true1, X$H2_deflated, pch = ifelse(X$source == "CIMMYT", 17, 19), col = ifelse(X$balanced, col1, col2),
     xlim = c(0, 1), ylim = c(0, 1), xlab = expression("Direct REML one-environment " * H^2),
     ylab = expression(H[single]^2 / (1 + rho)))
abline(0, 1, lty = 2, col = grey)
legend("topleft", c("Balanced (>=80% GxE cells filled)", "Unbalanced", "CIMMYT nursery"),
       pch = c(19, 19, 17), col = c(col1, col2, "grey30"), bty = "n", cex = 0.75)
dev.off()
cat("Figures written to figures/\n")
