# Risk, Crop Yields, and Weather Index Insurance in Village India: Replication Code
This repository contains the data and R (v.4.1.2) and WinBUGS (v.14) programs required to replicate the tables and figures in Michler, J.D., Viens, F.G., and Shively, G.E. (2022), as well as the associated appendices. The programs generate all summary tables and figures. In order to construct the final regression tables, some manual editing (and copying/pasting across files) is required, so we do not reproduce them here. However, all relevant numbers are in the regression output generated by the programs files listed below.

[![DOI](https://zenodo.org/badge/484186941.svg)](https://zenodo.org/badge/latestdoi/484186941)


The programs are as follows:
- `rci_analysis.R` replicates all tables and figures in the paper and appendices, including running subsidiary `.bug` programs for the Bayesian analysis.
- `india.0.bug` through `india.3.bug` which are subsidiary files called from `rci_analysis.R`. These files simply need to be openned in WinBUGS (v.14) in order to be called from the R script.

The data are included in `rci_data.dta`, a Stata file imported by `rci_analysis.R`.

Note that the code requires installing the correct version of R and WinBUGS as well as the R packages detailed at the top of `rci_analysis.R`.
