## =================================================
## Main script to run MA-GEV analysis
## =================================================

## -------------------------------------------------
## 1. Set working directory (if needed)
## -------------------------------------------------
# setwd("path/to/project")

## -------------------------------------------------
## 2. Load data
## -------------------------------------------------
load("data/haenam.RData")

## -------------------------------------------------
## 3. Source required R scripts
## -------------------------------------------------
source("R/cand.xi.new_05Dec25.R")
source("R/cand.xi.paper_12Dec25.R")

source("R/MA_gev_main_16Dec25.R")
source("R/mixed.gev.new_14Dec25.R")
source("R/weight.com.BMA_16Dec25.R")
source("R/K_sensitivity_plot_MA_GEV.R")

## -------------------------------------------------
## 4. Run analysis
## -------------------------------------------------

# Apply the ma.gev.open function
## ma.gev.open(data=haenam, quant=c(.98, .99, .995), weight='like', bma=TRUE)
