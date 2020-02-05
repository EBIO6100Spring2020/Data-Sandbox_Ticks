### ABOUT
## Script to download NEON Tick abundance data for EDA
## 30 Jan 2020
## WEM
library(dplyr)
library(neonUtilities)
### downloading NEON data #####
# Tick Abundance data are in DP1.10093.001
zipsByProduct(dpID = "DP1.10093.001", site = "all", package = "basic",
              savepath = "data_raw/") # downloads from NEON as zip
stackByTable("data_raw/filesToStack10093/") # merges zips into csv

tickField <- read.csv("data_raw/filesToStack10093/stackedFiles/tck_fielddata.csv")
head(tickField)
tickTaxon <- read.csv("data_raw/filesToStack10093/stackedFiles/tck_taxonomyProcessed.csv")
# tickTaxon %>% filter(family == "") %>%View()
# head(tickTaxon)
tickTaxon <- tickTaxon %>% select(namedLocation, plotID, collectDate, sampleID, subsampleID, scientificName,
                                  acceptedTaxonID, genus, family, subfamily, sexOrAge, individualCount)
write.csv(tickTaxon, "data_derived/Tick_Abundance_Taxonomy.csv", row.names = FALSE)
write.csv(tickField, "data_derived/Tick_Abundance_Counts.csv", row.names = FALSE)
