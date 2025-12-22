# MA-GEV

This repository provides R code for implementing a **Model Averaging (MA) approach for the Generalized Extreme Value (GEV) model**, including parameter estimation, model weighting, and sensitivity analysis.

The code is organized in a modular way. All core functions are defined in individual R scripts and are loaded using `source()` rather than installing an R package.

---

## How to run the analysis

1. Clone or download this repository.

2. Set the working directory to the project root (`MA-gev-main`).

```r
setwd("path/MA-gev-main")
```

3. Run the main script
```r
source("R/run_ma_gev.R")
```
