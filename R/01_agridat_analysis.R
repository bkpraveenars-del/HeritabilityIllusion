# 01_agridat_analysis.R -- single-environment vs one-environment (G x E removed) H2
# for 19 replicated multi-environment trial series in the agridat package.
source("R/00_setup.R")
set.seed(2026)

# Series retained after screening (see Methods). Excluded: gomez.multilocsplitplot
# (2 genotypes), belamkar.augmented (1 replicated environment), sharma.met
# (documented as possibly simulated).
spec <- read.table(text = "
ds g e r y
acorsi.grayleafspot gen env rep y
barrero.maize gen env rep yield
buntaran.wheat gen loc rep yield
damesa.maize gen site rep yield
dasilva.maize gen env rep yield
edwards.oats gen loc+year block yield
gauch.soy gen env rep yield
george.wheat gen loc+year block yield
kang.peanut gen env rep yield
lee.potatoblight gen year rep y
linder.wheat gen env block yield
lu.stability gen env block yield
omer.sorghum gen env rep yield
onofri.winterwheat gen year block yield
shafii.rapeseed gen loc+year rep yield
shaw.oats gen env+year block yield
tesfaye.millet gen site+year rep yield
vargas.wheat1.traits gen year rep yield
verbyla.lupin gen loc+site+year rep yield", header = TRUE, stringsAsFactors = FALSE)

getds <- function(name) { e <- new.env(); data(list = name, package = "agridat", envir = e); get(name, envir = e) }

prep <- function(s) {
  x  <- getds(s$ds); ec <- strsplit(s$e, "\\+")[[1]]
  d  <- data.frame(gen = factor(x[[s$g]]), env = factor(do.call(paste, c(x[ec], sep = "_"))),
                   rep = factor(x[[s$r]]), y = as.numeric(x[[s$y]]))
  d  <- droplevels(d[is.finite(d$y), ])
  d$rep <- factor(paste(d$env, d$rep))
  ok <- tapply(seq_len(nrow(d)), d$env, function(ix)
          length(unique(d$rep[ix])) >= 2 && length(unique(d$gen[ix])) >= 5)
  droplevels(d[d$env %in% names(ok)[ok], ])
}

fit_met <- function(d) {
  m <- suppressMessages(lmer(y ~ (1|env) + (1|rep) + (1|gen) + (1|gen:env), d, control = CTL))
  v <- vc(m); rb <- nrow(d) / (nlevels(droplevels(d$gen)) * nlevels(droplevels(factor(d$env))))
  list(v = v, rb = rb)
}
h2_true1 <- function(f) f$v[["gen"]] / (f$v[["gen"]] + f$v[["gen:env"]] + f$v[["Residual"]] / f$rb)

res <- list(); per <- list(); boot <- list()
for (i in seq_len(nrow(spec))) {
  s <- spec[i, ]; d <- prep(s)
  H1 <- c()
  for (en in levels(d$env)) {
    de <- droplevels(d[d$env == en, ])
    m  <- tryCatch(suppressMessages(lmer(y ~ rep + (1|gen), de, control = CTL)), error = function(e) NULL)
    if (is.null(m)) next
    v  <- vc(m); rr <- nrow(de) / nlevels(de$gen)
    H1[en] <- v["gen"] / (v["gen"] + v["Residual"] / rr)
    per[[length(per) + 1]] <- data.frame(ds = s$ds, env = en, H2_single = H1[en],
                                         sigma2_g = v["gen"], sigma2_e = v["Residual"], r = rr)
  }
  f  <- fit_met(d); v <- f$v; ne <- nlevels(d$env)
  res[[s$ds]] <- data.frame(ds = s$ds, ngen = nlevels(d$gen), nenv = ne, r = round(f$rb, 1),
    H2_single_mean = mean(H1, na.rm = TRUE), H2_single_median = median(H1, na.rm = TRUE),
    H2_true_1env = h2_true1(f),
    H2_MET = v[["gen"]] / (v[["gen"]] + v[["gen:env"]] / ne + v[["Residual"]] / (f$rb * ne)),
    sigma2_G = v[["gen"]], sigma2_GE = v[["gen:env"]], sigma2_e = v[["Residual"]],
    rho = v[["gen:env"]] / v[["gen"]],
    fill = mean(table(d$gen, d$env) > 0))
  # cluster bootstrap over environments
  envs <- levels(d$env); nb <- if (nrow(d) > 3000) min(40, BOOT_AGRIDAT) else BOOT_AGRIDAT
  infl <- replicate(nb, {
    sm <- sample(envs, replace = TRUE)
    db <- do.call(rbind, lapply(seq_along(sm), function(k) {
      z <- d[d$env == sm[k], ]; z$env <- paste0("b", k); z$rep <- paste0("b", k, z$rep); z }))
    fb <- tryCatch(fit_met(db), error = function(e) NULL)
    if (is.null(fb)) NA else mean(H1[sm], na.rm = TRUE) - h2_true1(fb) })
  boot[[s$ds]] <- data.frame(ds = s$ds, infl_lo = quantile(infl, .025, na.rm = TRUE),
                             infl_hi = quantile(infl, .975, na.rm = TRUE), nboot = sum(is.finite(infl)))
  message(s$ds, " done")
}
R <- do.call(rbind, res); rownames(R) <- NULL
R$inflation <- R$H2_single_mean - R$H2_true_1env
R <- merge(R, do.call(rbind, boot), by = "ds")
write.csv(R, "results/agridat_summary.csv", row.names = FALSE)
write.csv(do.call(rbind, per), "results/agridat_per_env.csv", row.names = FALSE)

cat("\nagridat: series =", nrow(R), " environments =", sum(R$nenv), "\n")
cat("median single-env H2 =", round(median(R$H2_single_mean), 2),
    " median true one-env H2 =", round(median(R$H2_true_1env), 2),
    " median inflation =", round(median(R$inflation), 2), "\n")
cat("bootstrap CI excludes 0:", sum(R$infl_lo > 0), "of", nrow(R), "\n")
