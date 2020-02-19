### About
## Script to download small mammal box trapping data from NEON
## 6 Feb 2020
## BKH (largely copied from WEM)
library(dplyr)
library(neonUtilities)

### downloading NEON data #####
# Small mammal trapping data are in DP1.10072.001
zipsByProduct(dpID = "DP1.10072.001", site = "all", package = "basic",
              savepath = "data_raw/") # downloads from NEON as zip
stackByTable("data_raw/filesToStack10072/") # merges zips into csv

# read in full dataset (takes ~60s)
SMamTraps <- read.csv("data_raw/filesToStack10072/stackedFiles/mam_pertrapnight.csv")

# keep only necessary columns
SMamTraps <- SMamTraps %>% select(domainID,
                                  siteID,
                                  plotID,
                                  nlcdClass,
                                  decimalLatitude,
                                  decimalLongitude,
                                  geodeticDatum,
                                  elevation,
                                  trapStatus,
                                  trapType,
                                  collectDate,
                                  tagID,
                                  taxonID,
                                  scientificName,
                                  taxonRank,
                                  nativeStatusCode,
                                  recapture,
                                  lifeStage,
                                  larvalTicksAttached,
                                  nymphalTicksAttached,
                                  adultTicksAttached
)

head(SMamTraps)
names(SMamTraps)

write.csv(SMamTraps, "data_derived/Small_Mammal_Trapping.csv", row.names=FALSE)

#### TEST BETTER DOWNLOAD WAY ####
test <- loadByProduct()
?loadByProduct
