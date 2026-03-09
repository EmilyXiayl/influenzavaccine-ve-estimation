# Influenza Vaccine Effectiveness (VE) Estimation

Bayesian hierarchical meta-analysis for influenza vaccine effectiveness: data preparation, JAGS fitting, country-level VE estimates, and LOCO / K-fold / LOSO cross-validation.

## Requirements

- **R** (4.x or later)
- **JAGS**: install from [JAGS](https://sourceforge.net/projects/mcmc-jags/)
- R packages: `HDInterval`, `ggplot2`, `matrixStats`, `dplyr`, `magrittr`, `readxl`, `stringr`, `tidyr`, `metafor`, `tidyverse`, `corrplot`, `rjags`, `coda`

Example:

```r
install.packages(c("HDInterval", "ggplot2", "matrixStats", "dplyr", "magrittr", 
                   "readxl", "stringr", "tidyr", "metafor", "tidyverse", "corrplot", "rjags", "coda"))
```

## Data and run

- Set working directory to the script folder (RStudio or `Rscript --file=` does this automatically).
- Place `model_input_NEWGDP.xlsx` in the same directory. Required sheets: `country`, `ve2`, `global`, `veNTD`.
- Optional: `VE_results_WB.xlsx`, validation `.rds` outputs (see script).

## Run order

1. Run the main model up to country-level VE; check for errors.
2. Run LOCO / K-fold / LOSO as needed (set `parallel = TRUE` in LOCO to speed up).

## License

If you use this code, please cite accordingly. Do not upload sensitive or unpublished data.
