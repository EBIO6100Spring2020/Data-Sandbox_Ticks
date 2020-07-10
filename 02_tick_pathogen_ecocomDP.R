# title: "Tick Pathogen ecocomDP"
# author: Melissa Chen, Wynne Moss, Brendan Hobart, Matt Bitters
# date: 9 July 2020

library("tidyverse")
library("dataCleanr")

#### Load ####
tick_path <- read.csv("data/tick_pathogen.txt")

### Required columns #### 
req_col_observation <- c("observation_id"
                         , "event_id", "package_id", "location_id","taxon_id"
                         , "observation_datetime"
                         , "variable_name","value","unit")
req_col_location <- c("location_id"
                      , "location_name"
                      , "latitude", "longitude", "elevation")

req_col_taxon <- c("taxon_id"
                   , "taxon_rank","taxon_name"
                   , "authority_system","authority_taxon_id")

req_col_location_ancillary <- c("location_ancillary_id", "location_id"
                                , "variable_name","value") # for nlcd class

req_col_observation_ancillary <- c("observation_ancillary_id", "event_id"
                                   , "variable_name", "value") # for host species

#### Observation table ####
colnames(tick_path)
dpID <- "DP1.10092.001" # package ID

c("observation_id"
  , "event_id", "package_id", "location_id","taxon_id"
  , "observation_datetime"
  , "variable_name","value","unit")

tick_path$observation_id <- tick_path$uid
tick_path$event_id <- tick_path$testingID
tick_path$package_id <- dpID
tick_path$location_id <- tick_path$namedLocation
tick_path$taxon_id <- tick_path$testPathogenName
tick_path$observation_datetime <- tick_path$collectDate
?iso8601_convert()## Need to check whether I need to convert for time zones
