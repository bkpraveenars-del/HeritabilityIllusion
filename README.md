# Heritability Illusion — code and results

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23013121.svg)](https://doi.org/10.5281/zenodo.23013121)

Reproducible R code for the manuscript:

> Kumar P. *Single-environment heritability overstates selection response: a 695-environment calibration across ten crops.* (to submit to the journal)

## Zenodo

Concept DOI: [10.5281/zenodo.23013121](https://doi.org/10.5281/zenodo.23013121)  
Published version: [10.5281/zenodo.23013120](https://doi.org/10.5281/zenodo.23013120)

The code reanalyses 19 replicated multi-environment trial series from the R package **agridat** and plot-level grain yield from four **CIMMYT Elite Spring Wheat Yield Trials** (26th, 27th, 29th and 38th ESWYT). For each trial series it computes:

- single-environment broad-sense heritability;
- the heritability of a one-environment test once genotype × environment variance is removed;
- leave-one-environment-out realised versus predicted selection gain;
- a check of the deflation rule H²/(1 + ρ).

## Folder layout

| Path | Contents |
| --- | --- |
| `run_all.R` | Runs every step in order |
| `R/00_setup.R` | Package installation, shared settings (bobyqa optimizer, bootstrap size) |
| `R/01_agridat_analysis.R` | agridat series: single-environment and multi-environment heritability, bootstrap |
| `R/02_cimmyt_clean.R` | Reads and cleans CIMMYT ESWYT plot data (removes exact duplicates and the local check) |
| `R/03_cimmyt_analysis.R` | CIMMYT heritability, leave-one-environment-out test, robustness subsets (`std` argument = yield standardised within environment) |
| `R/04_cimmyt_bootstrap.R` | Cluster bootstrap over environments for the inflation |
| `R/05_deflation_rule.R` | Tests H²/(1 + ρ) against direct estimates, by trial balance |
| `R/06_figures.R` | Figures 1–3 |
| `results/` | Output tables (CSV) as used in the manuscript |
| `figures/` | Figures 1–3 (PNG, 300 dpi) |
| `data/cimmyt/` | Place the CIMMYT raw data files here (not redistributed) |

## Data

**agridat** data load automatically from the package (Wright K, agridat, https://doi.org/10.32614/CRAN.package.agridat). The published analysis used development version 1.27 from https://github.com/kwstat/agridat; CRAN version 1.26 contains the same datasets.

**CIMMYT ESWYT** data are not redistributed here. Download them from the CIMMYT Research Data repository (https://data.cimmyt.org) and put these four files in `data/cimmyt/`:

- `26TH ESWYT_RawData.xls`
- `27TH ESWYT_RawData.xls`
- `29ESWYT_RawData.xls`
- `38TH ESWYT_RawData.xls` (38th ESWYT dataset: hdl:11529/10548343)

The files are tab-delimited text despite the `.xls` extension.

## How to run

```r
# R >= 4.3; packages lme4 and agridat are installed automatically
Rscript run_all.R                 # full run (about 1–2 h, mostly bootstraps)
HI_BOOT=5 Rscript run_all.R       # quick test with 5 bootstrap resamples (Linux/macOS)
```

On Windows, open R in this folder and run `source("run_all.R")`. For a quick test there, first run `Sys.setenv(HI_BOOT = 5)`.

## Notes on the analysis

- Every mixed model uses the `bobyqa` optimizer. The default optimizer converged to a spurious optimum in one fit.
- 1,300 exact duplicate plot records (13 environments of the 38th ESWYT) are removed before analysis.
- Entry 1 in each ESWYT is a local check that differs among sites, so it is excluded, leaving 49 common entries.
- Bootstrap intervals for ρ are biased low, because resampled duplicate environments are identical copies. Only intervals for the inflation are reported.
- Bootstrap results differ slightly between runs and machines. The seed is fixed at 2026.

## Licence

Code is released under the MIT licence (see `LICENSE`). Results tables and figures are released under CC BY 4.0.

## Citation

See `CITATION.cff`. Please also cite agridat, the original sources of the agridat datasets, and the CIMMYT ESWYT datasets.
