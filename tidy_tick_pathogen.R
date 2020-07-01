# Download and clean tick data 
# Authors:  Melissa Chen, Wynne Moss, Brendan Hobart, Matt Bitters
# Date: 6/30/20

library(tidyverse)
library(neonUtilities)
library(lubridate)


if(!file.exists("data_raw/tck_sitexspecies.Rdata") |
   !file.exists("data_raw/tck_sitexspecies_env.Rdata") |
   !file.exists("data_raw/tck_longform.Rdata")){
  loadByProduct(dpID = "DP1.10092.001", site = "all", package = "expanded", check.size = F) # downloads from NEON and loads to env

  zipsByProduct(dpID="DP1.10092.001", site="all", package = "expanded", check.size = F, savepath = "./data_raw/")
  stackByTable(filepath = "./data_raw/filesToStack10092/", savepath = "./data_raw/filesToStack10092/")
}

tick_path <- read.csv("data_raw/filesToStack10092/stackedFiles/tck_pathogen.csv", na.strings = c("","NA","na","n/a","N/A"))

# Checking uid is unique
nrow(tick_path) == length(unique(tick_path$uid))
# Check there isn't any missing data in location or time
any(is.na(tick_path[,c("domainID","siteID","plotID", "decimalLatitude", "decimalLongitude", "plotType", "nlcdClass", "elevation", "collectDate")]))

# Look at sample conditions
table(tick_path$sampleCondition) # We'll get rid of all "non-okay" samples

# How does named location differ from just plotID

# How many test results are NA?

# Filter NA test results out

# testPathogenName-- remove "HardTick DNA Quality"
# Any testPathogenName=="HardTick DNA Quality", testResults=="negative"?

# Check that all testingIDs are found in the original dataset





View(tick_path)
