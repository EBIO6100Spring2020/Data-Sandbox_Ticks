# Download and clean tick data 
# Authors: Wynne Moss, Melissa Chen, Brendan Hobart, Matt Bitters
# Date: 6/22/20

### to do list
# add total identified column
# add total tick column
# ecocomp
# fix the unassigned code and simplify 
# fix comments
devtools::install_github('NEONScience/NEON-utilities/neonUtilities')

library(tidyverse)
library(neonUtilities)
library(lubridate)

###########################################
#  LOAD TICK DATA FROM NEON OR FILE SOURCE 
###########################################

# Tick Abundance and field data data are in DP1.10093.001

# modify this code for GitHub repo structure (e.g. data raw folder name)

if(file.exists("data_raw/tck_fielddata.Rdata") & file.exists("data_raw/tck_taxonomyProcessed.Rdata")){
  load("data_raw/tck_fielddata.Rdata")
  load("data_raw/tck_taxonomyProcessed.Rdata")
} else {
  Tick_all <- loadByProduct(dpID = "DP1.10093.001",
                            package = "basic", check.size = F) # downloads from NEON and loads to env
  Tick_all <- list(tck_fielddata = 1, tck_taxonomyProcessed = 1, tck_taxonomyRaw = 1)
  tck_fielddata <- Tick_all$tck_fielddata
  tck_tax <- Tick_all$tck_taxonomyProcessed
  tck_tax_raw <- Tick_all$tck_taxonomyRaw
}

# tck_taxonomyProcessed and tck_fielddata are the two datasets we want 


##### Optional: fixing loadByProduct errors ####
  # #### If you get an error with regexp in list.files, and the above doesn't work:
  # zipsByProduct(dpID="DP1.10093.001", site="all", package = "expanded", check.size = F, savepath = ".//data_raw")
  # stackByTable(filepath = "./data_raw/filesToStack10093/", savepath = "./data_raw/filesToStack10093/") # in R-4.0.0, MacOS Mojave v10.14.6, it will stack but won't delete raw files.
  # #### manually delete "raw" files
  # fileNames <- list.files(path="./data_raw/filesToStack10093/", pattern = "NEON[.]", include.dirs = TRUE, full.names = TRUE)
  # unlink(fileNames, recursive = TRUE)
  # tck_fielddata <- read.csv("./data_raw/filesToStack10093/stackedFiles/tck_fielddata.csv", header=TRUE, stringsAsFactors = FALSE)
  # tck_tax <- read.csv("./data_raw/filesToStack10093/stackedFiles/tck_taxonomyProcessed.csv", header=TRUE, stringsAsFactors = FALSE)
  # tck_tax_raw <- read.csv("./data_raw/filesToStack10093/stackedFiles/tck_taxonomyRaw.csv", header=TRUE, stringsAsFactors = FALSE)
  # Tick_all <- list(tck_fielddata=(tck_fielddata),tck_taxonomyProcessed=(tck_tax),tck_taxonomyRaw=(tck_tax_raw))





###########################################
# CLEAN TICK FIELD DATA #
###########################################

#### NAs
# replace empty characters with NAs instead of ""
repl_na <- function(x) ifelse(x=="", NA, x)

tck_fielddata %>% mutate_if(.predicate = is.character, .funs = repl_na) -> tck_fielddata

# check which fields have NAs
# tck_fielddata %>% select(everything()) %>% summarise_all(funs(sum(is.na(.))))
tck_fielddata %>% select(everything()) %>% summarise_all(~sum(is.na(.))) # this is "not" depreciated MYC

#### Quality Flags
# remove samples that had logistical issues
# keep only those with NA in samplingImpractical 
tck_fielddata %>% filter(is.na(samplingImpractical)) -> tck_fielddata
# tck_fielddata %>% select(everything()) %>% summarise_all(funs(sum(is.na(.))))
tck_fielddata %>% select(everything()) %>% summarise_all(~sum(is.na(.))) # MYC

table(tck_fielddata$sampleCondition) # none of these are major issues
table(tck_fielddata$dataQF) # 551 are legacy data (they look OK but let's remove anyways)

tck_fielddata %>% filter(dataQF == "legacyData") %>% mutate(year = year(collectDate)) %>% pull(year) %>% table()
tck_fielddata %>% filter(is.na(dataQF)) -> tck_fielddata

# scan through remarks
# tck_fielddata %>% filter(is.na(dataQF)) %>% pull(remarks) %>% unique() 
# many of these explain why the tax count might be smaller than field count (individuals dropped or lost, etc)

##### Sample IDs
# check that all of the rows WITHOUT a sample ID have no ticks
tck_fielddata %>% filter(is.na(sampleID)) %>% pull(targetTaxaPresent) %>% table()

# check that all the rows WITH a sample ID have ticks
tck_fielddata %>% filter(!is.na(sampleID)) %>% pull(targetTaxaPresent) %>% table()

# sampleID only assigned when there were ticks present; all missing sample IDs are for drags with no ticks (good)
# option: filter out all of the sample events with no ticks 

# make sure sample IDs are unique
tck_fielddata %>% group_by(sampleID) %>% summarise(n = n()) %>% filter(n > 1) # no
sum(duplicated(na.omit(tck_fielddata$sampleID)))

# check that drags with ticks present have a sample ID
tck_fielddata %>% filter(targetTaxaPresent == "Y" & is.na(sampleID)) # none are missing S.ID

### Check Counts vs. NAs 
# some samples have ticks present and a sample ID but no counts 
tck_fielddata %>% filter(targetTaxaPresent == "Y") %>%
  filter(is.na(adultCount) & is.na(nymphCount) & is.na(larvaCount)) %>% pull(sampleID) -> missing.counts

# are they in the tax data?
tck_tax %>% filter(sampleID %in% missing.counts)
tck_tax_raw %>% filter(sampleID %in% missing.counts)

# these aren't in the tax table either...get rid of these records
tck_fielddata %>% filter(!(sampleID %in% missing.counts)) -> tck_fielddata
rm(missing.counts)

# some samples have no ticks present but counts are NA instead of 0
# fill these in with 0s

# tck_fielddata %>% mutate(adultCount = if_else(targetTaxaPresent == "N" & totalSampledArea > 0 & is.na(adultCount), 0, adultCount),
#                          nymphCount = if_else(targetTaxaPresent=="N" & totalSampledArea > 0 & is.na(nymphCount), 0, nymphCount),
#                          larvaCount = if_else(targetTaxaPresent == "N" & totalSampledArea >0 & is.na(larvaCount),0, larvaCount)) -> tck_fielddata

tck_fielddata %>% mutate(adultCount = ifelse(targetTaxaPresent == "N" & totalSampledArea > 0 & is.na(adultCount), 0, adultCount),
                         nymphCount = ifelse(targetTaxaPresent=="N" & totalSampledArea > 0 & is.na(nymphCount), 0, nymphCount),
                         larvaCount = ifelse(targetTaxaPresent == "N" & totalSampledArea >0 & is.na(larvaCount),0, larvaCount)) -> tck_fielddata
## Using base R ifelse instead of tidyverse if_else because if_else gives dbl/int errors for counts MYC


# are there any 0 counts where there should be > 0?
tck_fielddata %>% filter(targetTaxaPresent == "Y") %>% mutate(totalCount = adultCount + nymphCount + larvaCount) %>% filter(totalCount == 0) # no

### Check that all the field samples are in the tax table
tck_fielddata %>% filter(!sampleID %in% tck_tax$sampleID, !is.na(sampleID)) %>% nrow() # some are missing
# are there any in the tax table not in field table?
tck_tax %>% filter(!sampleID %in% tck_fielddata$sampleID, !is.na(sampleID)) %>% nrow() # some are missing (legacy)

# get rid of field samples that have no taxonomic info (not obvious why they weren't id'd)
tck_fielddata %>% filter(sampleID %in% tck_tax$sampleID | is.na(sampleID)) -> tck_fielddata 

# get rid of the tax samples that have no field data 
tck_tax %>% filter(sampleID %in% tck_fielddata$sampleID) -> tck_tax

### Confirm no missing data
# list of fields that shouldn't have NAs
req_cols <- c("siteID", "plotID", "collectDate", "adultCount", "nymphCount", "larvaCount", "totalSampledArea")

# check for NAs again
# tck_fielddata %>% select(req_cols) %>% summarise_all(funs(sum(is.na(.))))
tck_fielddata %>% select(req_cols) %>% summarise_all(~sum(is.na(.))) # MYC


####################################
# CLEAN TICK TAXONOMY DATA 
####################################

### replace "" with NA
tck_tax %>% mutate_if(.predicate = is.character, .funs = repl_na) -> tck_tax

### check that sample IDs match across tick field and tick tax

# are all the tax sample IDs in the field table?
tck_tax %>% filter(!sampleID %in% tck_fielddata$sampleID) %>% nrow() #yes

### check sample quality flags
table(tck_tax$sampleCondition) 
# none of these are deal breakers but just keep those deemed OK
tck_tax %>% filter(sampleCondition == "OK") -> tck_tax

### check for NAs in fields 
tck_tax %>% select(everything()) %>% summarise_all(funs(sum(is.na(.))))
tck_tax %>% select(everything()) %>% summarise_all(~sum(is.na(.))) # MYC 

tck_tax %>% filter(!is.na(identifiedDate)) %>% filter(!is.na(acceptedTaxonID)) -> tck_tax

# epithet are all NAs (don't need these columns)
table(tck_tax$infraspecificEpithet)

# qualifier not sure how to interpret this
table(tck_tax$identificationQualifier)

# legacy data is the only flag (OK)
table(tck_tax$dataQF) 

# minor remarks
unique(tck_tax$remarks)

### check that the IDs all make sense
table(tck_tax$acceptedTaxonID)
table(tck_tax$taxonRank)
table(tck_tax$family)

tck_tax %>% filter(taxonRank == "order") %>% pull(acceptedTaxonID) %>% table()# some were just ID'd to order and all are IXOSP2
tck_tax %>% filter(taxonRank == "family") %>% pull(acceptedTaxonID) %>% table()
tck_tax %>% filter(taxonRank == "species") %>% pull(acceptedTaxonID) %>% table()

### sex or age columns
table(tck_tax$sexOrAge)
tck_tax %>% filter(is.na(sexOrAge)) # all have sex or age
tck_tax %>% filter(individualCount==0 | is.na(individualCount)) %>% nrow() # all have counts

# create a lifestage column so tax counts can be compared to field data
tck_tax %>% mutate(lifeStage = case_when(sexOrAge == "Male" | sexOrAge == "Female" ~ "Adult",
                                         sexOrAge == "Larva" ~ "Larva",
                                         sexOrAge == "Nymph" ~ "Nymph")) -> tck_tax
# this assumes that any ticks that were sexed must be adults (I couldn't confirm this)

#########################################
#  MERGE FIELD AND TAXONOMY DATA # 
#########################################

# there are a few options here:
# 1) left join tick_field to tick_tax: this will keep all the 0s but requires some wide/long manipulation
# 2) left join tick_tax to tick_field: this will only keep the records where ticks were ID'd. could add in 0s after depending on preference

# We will left join tick_field to tick_tax, so that 0s are preserved.
# Doing so will require us to resolve count discrepancies (this may be overkill for those interested in P/A)
# For simplicity we will also collapse M/F into all adults (also optional)


# first check for dup colnames
intersect(colnames(tck_fielddata), colnames(tck_tax))

# rename columns containing data unique to each dataset
colnames(tck_fielddata)[colnames(tck_fielddata)=="uid"] <- "uid_field"
colnames(tck_fielddata)[colnames(tck_fielddata)=="sampleCondition"] <- "sampleCondition_field"
colnames(tck_fielddata)[colnames(tck_fielddata)=="remarks"] <- "remarks_field"
colnames(tck_fielddata)[colnames(tck_fielddata)=="dataQF"] <- "dataQF_field"
colnames(tck_fielddata)[colnames(tck_fielddata)=="publicationDate"] <- "publicationDate_field"

colnames(tck_tax)[colnames(tck_tax)=="uid"] <- "uid_tax"
colnames(tck_tax)[colnames(tck_tax)=="sampleCondition"] <- "sampleCondition_tax"
colnames(tck_tax)[colnames(tck_tax)=="remarks"] <- "remarks_tax"
colnames(tck_tax)[colnames(tck_tax)=="dataQF"] <- "dataQF_tax"
colnames(tck_tax)[colnames(tck_tax)=="publicationDate"] <- "publicationDate_tax"

# in order to merge counts with field data we will first combine counts from male/female in the tax column into one adult class
# if male/female were really of interest you could skip this (but the table will just be much wider since all will have a M-F option)

tck_tax %>% group_by(sampleID, acceptedTaxonID, lifeStage) %>% summarise(individualCount = sum(individualCount, na.rm = TRUE)) -> tck_tax_lifestage

# make tick taxonomy wide; now only one row per sample id
tck_tax_lifestage %>% pivot_wider(id_cols = sampleID,  names_from = c(acceptedTaxonID, lifeStage), values_from = individualCount) -> tck_tax_wide

# replace NAs with 0s
tck_tax_wide %>% replace(is.na(.), 0) -> tck_tax_wide
# each row is now be a unique sample id

sum(duplicated(tck_tax_wide$sampleID)) # none are duplicated
length(unique(tck_tax_wide$sampleID)) ==length(unique(na.omit(tck_fielddata$sampleID)))  # some are missing from either side

tck_fielddata %>% filter(sampleID %in% tck_tax_wide$sampleID | is.na(sampleID)) -> tck_fielddata
tck_tax_wide %>% filter(sampleID %in% tck_fielddata$sampleID) -> tck_tax_wide

tck_fielddata %>% left_join(tck_tax_wide, by = "sampleID") -> tck_merged

### clean up 0s
# if ticks were not found, counts should be 0
tck_merged %>% select(contains(c("_Nymph", "_Adult" , "_Larva"))) %>% colnames() -> taxon.cols # columns containing counts per taxon

for(i in 1:length(taxon.cols)){
  tck_merged[which(tck_merged$targetTaxaPresent =="N"), which(colnames(tck_merged) == taxon.cols[i])] = 0
}

# check for NAs in the count columns

# tck_merged %>% select_if(is_numeric)%>% summarise_all(funs(sum(is.na(.))))
tck_merged %>% select_if(is_numeric)%>% summarise_all(~sum(is.na(.))) ### MYC

#########################################################
#FIX COUNT DISCREPANCIES BETWEEN FIELD AND LAB #
#########################################################

#### FIX COUNT 1.1 FLAG DISCREPANCIES #####
# first flag those where the counts might not match up
tck_merged %>% select(contains("_Adult")) %>% colnames() -> adult_cols
tck_merged %>% select(contains("_Nymph")) %>% colnames() -> nymph_cols
tck_merged %>% select(contains("_Larva")) %>% colnames() -> larva_cols

tck_merged %>% mutate(
  totalAdult_tax = rowSums(.[adult_cols]),
  totalNymph_tax = rowSums(.[nymph_cols]),
  totalLarva_tax = rowSums(.[larva_cols]),
  totalCount_tax = rowSums(.[c(adult_cols, nymph_cols, larva_cols)]),
  totalCount_field = rowSums(.[c("adultCount", "nymphCount", "larvaCount")])
) -> tck_merged

# get a list of columns with counts
tck_merged %>% select(contains(c("adult", "nymph","larva", "total", "Adult",  "Nymph",  "Larva", "Count"), ignore.case = FALSE)) %>% select(-totalSampledArea) %>% colnames() -> countCols

# create a count flag column
tck_merged %>% mutate(CountFlag = NA) -> tck_merged

# flag 1: the total count matches, but the life stage columns don't (error 1)
tck_merged$CountFlag[tck_merged$totalCount_field==tck_merged$totalCount_tax] <- 1

# no flag:  all the counts match (no flag) -- some of the 1s will now be  0s
tck_merged$CountFlag[tck_merged$nymphCount==tck_merged$totalNymph_tax &
                       tck_merged$adultCount == tck_merged$totalAdult_tax &
                       tck_merged$larvaCount == tck_merged$totalLarva_tax] <- 0

# flag 2: the field total is > than tax total
tck_merged$CountFlag[tck_merged$totalCount_field > tck_merged$totalCount_tax] <- 2

# flag 3: the field total is < than the tax total
tck_merged$CountFlag[tck_merged$totalCount_field < tck_merged$totalCount_tax] <- 3 

### reconcile count discrepancies
table(tck_merged$CountFlag)
sum(is.na(tck_merged$CountFlag))

# the vast majority of counts match. 

# add columns to lump unidentified counts
tck_merged %>% mutate(UNIDENTIFIED_Larva = 0, UNIDENTIFIED_Nymph = 0, UNIDENTIFIED_Adult = 0) -> tck_merged

#### FIX COUNT 2.1 RECONCILE DISCREPANCIES WHERE TOTALS MATCH #####

## Flag Type 1: the total count matches
tck_merged %>% filter(CountFlag == 1) %>% select(countCols) 

# many of these are adults misidentified as nymphs or vice versa
# trust the taxonomy counts here (ignore nymphCount, adultCount, larvaCount from field data)

#### FIX COUNT 2.2 RECONCILE DISCREPANCIES WHERE FIELD COUNT > TAX COUNT ####
# more ticks found in the field than in the lab
# in some cases, there were too many larvae and not all were ID'd
# in other cases, there may be issues with the counts and some ended up not being ticks
# if we can't rectify, lump the excess field ticks into 'unidentified'

### first look at those where just the larval count is off (e.g. more larvae in field than in lab, but other counts match)
tck_merged %>% filter(CountFlag==2) %>% filter(larvaCount!=totalLarva_tax,
                                               adultCount == totalAdult_tax,
                                               nymphCount ==  totalNymph_tax) %>% select(countCols)

# how many cases are there?
tck_merged %>% filter(CountFlag==2) %>% filter(larvaCount!=totalLarva_tax,
                                               adultCount == totalAdult_tax,
                                               nymphCount ==  totalNymph_tax) %>% pull(sampleID) -> larva.fix
length(larva.fix) 

# unidentified larvae: larvaCount - larvaCount_tax
tck_merged$UNIDENTIFIED_Larva[which(tck_merged$sampleID %in% larva.fix)] <- tck_merged$larvaCount[which(tck_merged$sampleID %in% larva.fix)]  - tck_merged$totalLarva_tax[which(tck_merged$sampleID %in% larva.fix)] 
tck_merged$CountFlag[which(tck_merged$sampleID %in% larva.fix)] <- "2f" # remove flag
rm(larva.fix)

### look at those where just the nymph count is off (e.g. more nymph in field than lab but other counts match)
tck_merged %>% filter(CountFlag==2) %>% filter(larvaCount==totalLarva_tax,
                                               adultCount == totalAdult_tax,
                                               nymphCount !=  totalNymph_tax) %>% select(countCols)

# how many cases are there?
tck_merged %>% filter(CountFlag==2) %>% filter(larvaCount==totalLarva_tax,
                                               adultCount == totalAdult_tax,
                                               nymphCount !=  totalNymph_tax) %>% pull(sampleID) -> nymph.fix
length(nymph.fix) 

# unidentified nymph: nymph count - nymph count tax
tck_merged$UNIDENTIFIED_Nymph[which(tck_merged$sampleID %in% nymph.fix)] <- tck_merged$nymphCount[which(tck_merged$sampleID %in% nymph.fix)]  - 
  tck_merged$totalNymph_tax[which(tck_merged$sampleID %in% nymph.fix)] 

tck_merged$CountFlag[which(tck_merged$sampleID %in% nymph.fix)] <- "2f" # remove flag
rm(nymph.fix)

### look at those where just the adult count is off (e.g. more adult in field than lab but other counts match)
tck_merged %>% filter(CountFlag==2) %>% filter(larvaCount==totalLarva_tax,
                                               adultCount != totalAdult_tax,
                                               nymphCount ==  totalNymph_tax) %>% select(countCols)

# how many cases are there?
tck_merged %>% filter(CountFlag==2) %>% filter(larvaCount==totalLarva_tax,
                                               adultCount != totalAdult_tax,
                                               nymphCount ==  totalNymph_tax) %>% pull(sampleID) -> adult.fix
length(adult.fix) 

# unidentified adult: adult count - adult count tax
tck_merged$UNIDENTIFIED_Adult[which(tck_merged$sampleID %in% adult.fix)] <- tck_merged$adultCount[which(tck_merged$sampleID %in% adult.fix)]  - 
  tck_merged$totalAdult_tax[which(tck_merged$sampleID %in% adult.fix)] 

tck_merged$CountFlag[which(tck_merged$sampleID %in% adult.fix)] <- "2f" # remove flag
rm(adult.fix)

### cases where adults might be mistaken for nymphs
tck_merged %>% filter(CountFlag==2) %>% filter(adultCount > totalAdult_tax,
                                               nymphCount < totalNymph_tax,
                                               larvaCount == totalLarva_tax) 
# no cases but in here, the excess adults that are not in nymph tax would be classified as unidentified adults


### cases where some larvae are actually nymphs
tck_merged %>% filter(CountFlag==2) %>% filter(adultCount == totalAdult_tax,
                                               nymphCount < totalNymph_tax,
                                               larvaCount > totalLarva_tax)  %>% pull(sampleID) -> larvae.fix
# unidentified larvae : larvaCount (field) - totalLarva_tax + (nymphCount - totalNymph_tax)
tck_merged$UNIDENTIFIED_Larva[which(tck_merged$sampleID %in% larvae.fix)] <- (tck_merged$larvaCount[which(tck_merged$sampleID %in% larvae.fix)]  - tck_merged$totalLarva_tax[which(tck_merged$sampleID %in% larvae.fix)]) +
  (tck_merged$nymphCount[which(tck_merged$sampleID %in% larvae.fix)] - tck_merged$totalNymph_tax[which(tck_merged$sampleID %in% larvae.fix)])

tck_merged$CountFlag[which(tck_merged$sampleID %in% larvae.fix)] <- "2f" # remove flag
rm(larvae.fix)

### cases where some nymphs are actually larvae
tck_merged %>% filter(CountFlag==2) %>% filter(adultCount == totalAdult_tax,
                                               nymphCount > totalNymph_tax,
                                               larvaCount < totalLarva_tax)  %>% pull(sampleID) -> nymph.fix
# unidentified nymph : nymphCount (field) - totalnymph_tax + (larvaeCount - larvae_tax)
# number of nymphs not counted minus the excess larvae
tck_merged$UNIDENTIFIED_Nymph[which(tck_merged$sampleID %in% nymph.fix)] <- (tck_merged$nymphCount[which(tck_merged$sampleID %in% nymph.fix)]  - tck_merged$totalNymph_tax[which(tck_merged$sampleID %in% nymph.fix)]) +
  (tck_merged$larvaCount[which(tck_merged$sampleID %in% nymph.fix)] - tck_merged$totalLarva_tax[which(tck_merged$sampleID %in% nymph.fix)])

tck_merged$CountFlag[which(tck_merged$sampleID %in% nymph.fix)] <- "2f" # remove flag
rm(nymph.fix)

### cases where some nymphs are actually adults
tck_merged %>% filter(CountFlag==2) %>% filter(adultCount < totalAdult_tax,
                                               nymphCount > totalNymph_tax,
                                               larvaCount == totalLarva_tax)  %>% pull(sampleID) -> nymph.fix
# unidentified nymph : nymphCount (field) - totalnymph_tax + (adultCount - totalAdult_tax)
# number of nymphs not counted minus the excess adults
tck_merged$UNIDENTIFIED_Nymph[which(tck_merged$sampleID %in% nymph.fix)] <- (tck_merged$nymphCount[which(tck_merged$sampleID %in% nymph.fix)]  - tck_merged$totalNymph_tax[which(tck_merged$sampleID %in% nymph.fix)]) +
  (tck_merged$adultCount[which(tck_merged$sampleID %in% nymph.fix)] - tck_merged$totalAdult_tax[which(tck_merged$sampleID %in% nymph.fix)])

tck_merged$CountFlag[which(tck_merged$sampleID %in% nymph.fix)] <- "2f" # remove flag
rm(nymph.fix)


### there are under-counts in adults and none are in the nymph/larval stage
tck_merged %>% filter(CountFlag==2) %>% filter(adultCount > totalAdult_tax, # there are unidentified adults but none are nymph or larvae
                                               nymphCount == totalNymph_tax | nymphCount > totalNymph_tax, larvaCount == totalLarva_tax | larvaCount > totalLarva_tax) %>%pull(sampleID) -> adult.fix

tck_merged$UNIDENTIFIED_Adult[which(tck_merged$sampleID %in% adult.fix)] <- (tck_merged$adultCount[which(tck_merged$sampleID %in% adult.fix)]  - tck_merged$totalAdult_tax[which(tck_merged$sampleID %in% adult.fix)]) 

tck_merged$CountFlag[which(tck_merged$sampleID %in% adult.fix)] <- "2f" # remove flag
rm(adult.fix)

### there are under-counts in nymphs and none are in the adult/larval stage
tck_merged %>% filter(CountFlag==2) %>% filter(nymphCount > totalNymph_tax, # there are unidentified adults but none are nymph or larvae
                                               adultCount == totalAdult_tax | adultCount > totalAdult_tax, larvaCount == totalLarva_tax | larvaCount > totalLarva_tax) %>%pull(sampleID) -> nymph.fix

tck_merged$UNIDENTIFIED_Nymph[which(tck_merged$sampleID %in% nymph.fix)] <- (tck_merged$nymphCount[which(tck_merged$sampleID %in% nymph.fix)]  - tck_merged$totalNymph_tax[which(tck_merged$sampleID %in% nymph.fix)]) 

tck_merged$CountFlag[which(tck_merged$sampleID %in% nymph.fix)] <- "2f" # remove flag
rm(nymph.fix)


### there are under-counts in larvae and none are in the adult/nymph stage
tck_merged %>% filter(CountFlag==2) %>% filter(larvaCount > totalLarva_tax, # there are unidentified adults but none are nymph or larvae
                                               adultCount == totalAdult_tax | adultCount > totalAdult_tax, nymphCount == totalNymph_tax | nymphCount > totalNymph_tax) %>%pull(sampleID) -> larva.fix

tck_merged$UNIDENTIFIED_Larva[which(tck_merged$sampleID %in% larva.fix)] <- (tck_merged$larvaCount[which(tck_merged$sampleID %in% larva.fix)]  - tck_merged$totalLarva_tax[which(tck_merged$sampleID %in% larva.fix)]) 

tck_merged$CountFlag[which(tck_merged$sampleID %in% larva.fix)] <- "2f" # remove flag
rm(nymph.fix)



### cases where some larvae are actually nymphs or adults
tck_merged %>% filter(CountFlag==2) %>% filter(adultCount < totalAdult_tax,
                                               nymphCount < totalNymph_tax,
                                               larvaCount > totalLarva_tax)  %>% pull(sampleID) -> larvae.fix
# unidentified larvae : larvaCount (field) - total_tax
tck_merged$UNIDENTIFIED_Larva[which(tck_merged$sampleID %in% larvae.fix)] <- (tck_merged$larvaCount[which(tck_merged$sampleID %in% larvae.fix)]  - tck_merged$totalCount_tax[which(tck_merged$sampleID %in% larvae.fix)])
tck_merged$CountFlag[which(tck_merged$sampleID %in% larvae.fix)] <- "2f" # remove flag
rm(adult.fix)


table(tck_merged$CountFlag)

#### FIX COUNT 2.3 RECONCILE DISCREPANCIES WHERE FIELD COUNT < TAX COUNT

## case 1) counts are off by minor numbers (trust tax count) 
# if the discrepancy is less than 20% of the total count or only 1 or two individuals
tck_merged %>% filter(CountFlag == 3) %>% filter(abs(totalCount_tax-totalCount_field)/totalCount_tax < 0.2 |
                                                   (totalCount_tax - totalCount_field)<=1) %>% pull(sampleID) -> ignore.diff

tck_merged$CountFlag[which(tck_merged$sampleID %in% ignore.diff)] <- "3f"
rm(ignore.diff)

tck_merged %>% filter(CountFlag == 3) # in these cases, fine to also trust the tax counts (higher)

table(tck_merged$CountFlag)

# all issues are resolved.

# note that some of these could have a "most likely" ID assigned based on what the majority of larval IDs are but the user can decide what to do with unidentified.


#####################################################
# RESHAPE FINAL DATASET
######################################################

# will depend on the end user but creating two options
# note that for calculation of things like richness, might want to be aware of the level of taxonomic resolution

### option 1) long form data including 0s

# get rid of excess columns (user discretion)
tck_merged_final <- tck_merged %>% select(-CountFlag, -totalAdult_tax, -totalNymph_tax, -totalLarva_tax, -totalCount_tax, -totalCount_field,
                                          -nymphCount, -larvaCount, -adultCount, -coordinateUncertainty, -elevationUncertainty,
                                          -samplingImpractical, -measuredBy, -dataQF_field, -publicationDate_field)
# long form data
tck_merged_final %>% pivot_longer(cols = c(taxon.cols, UNIDENTIFIED_Larva, UNIDENTIFIED_Nymph, UNIDENTIFIED_Adult), names_to = "Species_LifeStage", values_to = "IndividualCount") %>%
  separate(Species_LifeStage, into = c("acceptedTaxonID", "LifeStage"), sep = "_", remove = FALSE) -> tck_merged_final


# add in taxonomic resolution
table(tck_tax$acceptedTaxonID, tck_tax$taxonRank) # accepted TaxonID is either at family, order, or species level
tck_tax %>% group_by(acceptedTaxonID) %>% summarise(taxonRank =  first(taxonRank)) -> tax_rank

tck_merged_final %>% left_join(tax_rank, by = "acceptedTaxonID") -> tck_merged_final


### option 2) a sample x species table with only counts, summed by life form
tck_merged_final %>% select(uid_field, acceptedTaxonID, IndividualCount) %>% group_by(uid_field, acceptedTaxonID) %>%
  summarise(IndividualCount =  sum(IndividualCount)) %>% pivot_wider(id_cols = uid_field, names_from = acceptedTaxonID, values_from = IndividualCount) -> tck_site_species_table

tck_fielddata %>% filter(uid_field %in% tck_site_species_table$uid_field) %>% select(-adultCount, -nymphCount, -larvaCount)-> tck_site_species_env


########################################
# EXPORT DATA #
#######################################

saveRDS(tck_site_species_table, "data_raw/tck_sitexspecies.Rdata") # change to GitHub approprriate folder
saveRDS(tck_site_species_env, "data_raw/tck_sitexspecies_env.Rdata")
saveRDS(tck_merged_final, "data_raw/tck_longform.Rdata")

########################################
# NOTES FOR END USERS  #
#######################################

# be aware of higher order taxa ID when looking at richness (e.g. can't just do row sums)
# unidentified could be semi-identified by looking at most likely ID (e.g. if 200 of 700 larvae are IXOSP, very likely that the Unidentified larvae are IXOSP)
# unidentified coudl be changed to IXOSPP across the board 

