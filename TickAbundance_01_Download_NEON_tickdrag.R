### ABOUT
## Script to download NEON Tick abundance data for EDA
## 30 Jan 2020
## WEM

library(neonUtilities)
### downloading NEON data #####
# Tick Abundance data are in DP1.10093.001
zipsByProduct(dpID = "DP1.10093.001", site = "all", package = "basic",
              savepath = "data_raw/") # downloads from NEON as zip
stackByTable("data_raw/filesToStack10093/") # merges zips into csv

tickField <- read.csv("data_raw/filesToStack10093/stackedFiles/tck_fielddata.csv")
head(tickField)
