### ABOUT
## Script to explore NEON Tick Pathogen data
## 30 Jan 2020
## WEM

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
#### read in data #####
tickPathogen_raw <- read.csv("data_raw/filesToStack10092/stackedFiles/tck_pathogen.csv")
colnames(tickPathogen_raw)

# remove unneeded columns
tp <- tickPathogen_raw %>% select(subsampleID, domainID, siteID, plotID,
                                  nlcdClass, decimalLatitude,
                                  decimalLongitude, elevation,
                                  collectDate, testingID, individualCount,
                                  sampleCondition,testResult, testPathogenName)


#### explore individual columns #####
# for info see variables file in the data directory
# relevant columns:
# testingID (individual tick ID)
# siteID (where tick was collected)
# collectDate (when tick was collected)
# testResult (not tested, negative, positive)
# testPathogenName (which pathogen result corresponds to)

table(tp$domainID)
table(tp$siteID)
length(unique(tp$subsampleID))
length(unique(tp$testingID))
table(tp$testPathogenName, tp$testResult)
table(tp$sampleCondition)
# looks like Borrelia spp have best data; think about using this rather than a specific species 
#### filter low quality data/errors #####
tp <- tp %>% filter(!is.na(testPathogenName), sampleCondition == "OK")


#### re-format and aggregate data #######
# pathogen results: make not tested NAs and this column numeric
levels(tp$testResult) <- c(NA, "0", "1")
tp$testResult <- as.numeric(as.character(tp$testResult))

# remove not tested
tp <- tp %>% filter(!is.na(testResult))

table(tp$testResult)

# get date columns
tp$Date <- as.Date(tp$collectDate)
tp$Year <- year(tp$Date)
tp$Month <- month(tp$Date)

# separate the testingID into other columns
# some of this ID is redundant with other columns so we'll toss them from table
tp <- tp %>% separate(col = subsampleID, sep = "\\.", into = c("throw", "throw2", "Taxonomy", "Lifestage"), remove = FALSE) %>% select(-"throw", -"throw2")



# make tick data in wide format, with each row a single tick, each column a pathogen
# note new dplyr function pivot_wider() which replaces spread()


# first aggregate the non-pathogen test columns (other covariates)
tp_wide_cols <- tp %>% group_by(testingID) %>% summarise_at(vars(subsampleID:sampleCondition, Date:Month), first) 

# now re-arrange so that each row is a tick and there are multiple columns for each pathogen ("wide" format), join with other covariates
tp_wide <- tp %>% pivot_wider(id_cols = testingID, names_from = testPathogenName,values_from = testResult) %>% left_join(tp_wide_cols, by = "testingID")

rm(tp_wide_cols)
tp_wide <- tp_wide %>% select(-`HardTick DNA Quality`)


# aggregate by population
tp_wide %>% group_by(siteID, Date, Taxonomy) %>% 
  summarise(Borrelia.Prev = sum(`Borrelia sp.`)/n(),
            nlcdClass = first(nlcdClass),
            elevation = first(elevation),
            Year = first(Year),
            Month = first(Month)) -> tp_agg

#### visualize data ######
ggplot(data = tp_agg, aes(x=Date, y = Borrelia.Prev, group = siteID))+
  geom_line(aes(color = siteID))

ggplot(data = tp_agg, aes(x=Year, y = Borrelia.Prev)) + 
  geom_jitter(aes(group=siteID, color = siteID), alpha = .7)

ggplot(data = tp_agg, aes(x=Month, y = Borrelia.Prev)) + 
  geom_jitter(aes(group=siteID), alpha = .2)

ggplot(data = tp_agg, aes(x=nlcdClass, y = Borrelia.Prev)) + 
  geom_violin()
ggplot(data = tp_agg, aes(x=siteID, y = Borrelia.Prev)) + 
  geom_violin()
