### ABOUT
## Script to download NEON Tick data for EDA
## 30 Jan 2020
## WEM
library(neonUtilities)
### downloading NEON data #####
# Tick Pathogen data are in DP1.10092.001
zipsByProduct(dpID = "DP1.10092.001", site = "all", package = "basic",
              savepath = "data_raw/") # downloads from NEON as zip
stackByTable("./data_raw/filesToStack10092/") # merges zips into csv

tickPathogen <- read.csv("./data_raw/filesToStack10092/stackedFiles/tck_pathogen.csv")
head(tickPathogen)

# This is melissa's version