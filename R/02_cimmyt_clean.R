# 02_cimmyt_clean.R -- read and clean CIMMYT ESWYT plot-level grain yield data.
# Place the four RawData files downloaded from the CIMMYT Research Data repository
# (data.cimmyt.org) in data/cimmyt/ (see README). They are tab-delimited text
# despite the .xls extension.
source("R/00_setup.R")

files <- c("26TH ESWYT_RawData.xls", "27TH ESWYT_RawData.xls",
           "29ESWYT_RawData.xls",    "38TH ESWYT_RawData.xls")
paths <- file.path("data/cimmyt", files)
if (!all(file.exists(paths))) stop("Missing CIMMYT files in data/cimmyt/: ",
                                    paste(files[!file.exists(paths)], collapse = ", "))

D <- do.call(rbind, lapply(seq_along(paths), function(i) {
  x <- read.delim(paths[i], stringsAsFactors = FALSE, check.names = FALSE, quote = "", fileEncoding = "latin1")
  x <- x[x[["Trait name"]] == "GRAIN_YIELD", ]
  data.frame(nursery = sub("29ESWYT", "29TH ESWYT", sub("_RawData.xls", "", files[i])),
             occ = x$Occ, country = x$Country, loc = x$Loc_desc, cycle = x$Cycle,
             gno = x$Gen_no, gname = x$Gen_name, rep = x$Rep, blk = x$Sub_block, plot = x$Plot,
             y = suppressWarnings(as.numeric(x$Value)), unit = x$Unit)
}))

n0  <- nrow(D)
dup <- duplicated(D[, c("nursery", "occ", "gno", "rep", "blk", "plot", "y")])
cat("Grain-yield plot records:", n0, "\n")
cat("Exact duplicate records removed:", sum(dup), "in",
    length(unique(paste(D$nursery, D$occ)[dup])), "environments\n")
D <- D[!dup, ]
cat("Missing yield:", sum(!is.finite(D$y)), "\n")

# Entry 1 is a local check that differs among sites -> excluded (49 common entries)
D <- D[D$gno != 1 & is.finite(D$y) & D$y > 0, ]
D$env <- paste(D$nursery, D$occ)
D$gen <- factor(paste(D$nursery, D$gno))
saveRDS(D, "results/cimmyt_clean.rds")
cat("Environments before filtering:", length(unique(D$env)), "\n")
