## ABOUT
## Script to explore NEON small mammal trapping data
## 12 Feb 2020
## BKH

# need to do:
# - how to deal w/ mult. rows when mult. capture events happen
# - look into taxonID, how to deal

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
rm(list=ls())

#### load "derived" data, which just has unneccesary columns removed ####
mamdat <- read.csv("data_derived/Small_Mammal_Trapping.csv")

#### about the data ####
# each row is a trap night
# re-trapped individuals are known and indicated

#### important columns to note: ####
# siteID: NEON name for sampling location
# trapStatus: categorical description
#  0 = no data
#  1 = trap not set
#  2 = trap disturbed/door closed but empty
#  3 = trap door open or closed with spoor left
#  4 = >1 capture in a given night
#  5 = capture
#  6 = trap set and empty
# collectDate: date of trapping
# taxonID: species code (see cheat sheet below)
# recapture: indicates if capture was recap; helpful if doing abundance
# larval/nymphal/adultTicksAttached columns: need to explore

#### make taxonID cheat sheet ####
mam_taxon <- mamdat[13:14] %>% distinct() %>% arrange(taxonID)

#### get date/month/year columns ####
mamdat$Date <- as.Date(mamdat$collectDate)
mamdat$Month <- month(mamdat$Date)
mamdat$Year <- year(mamdat$Date)

#### explore columns ####
table(mamdat$siteID)
table(mamdat$trapStatus) # can remove "" and "1"s
table(mamdat$Month) # most observations May thru Oct
table(mamdat$Year) # 2013 a little thin; 2014 - 2019 good
table(mamdat$taxonID) # wide range among species
table(mamdat$larvalTicksAttached) # more "U"s than "Y"s
table(mamdat$nymphalTicksAttached) # more "U"s than "Y"s
table(mamdat$adultTicksAttached) # more "U"s than "Y"s

#### add effort column ####
# here, effort == # traps/night
# trap being set corresponds to trapStatus 2:6
mamdat$trapset <- "1"
mamdat$trapset[mamdat$trapStatus==""] <- "0"
mamdat$trapset[mamdat$trapStatus=="1 - trap not set"] <- "0"
mamdat$trapset <- as.numeric(mamdat$trapset)
mamdat$plot_date <- paste(mamdat$plotID, mamdat$Date) # add plot-date column
mamdat$plot_date <- as.factor(mamdat$plot_date)
# mutate to add column with # trap-nights per trapping event per plot
mamdat <- mamdat %>% group_by(plot_date) %>%
  mutate(sum.traps = sum(trapset)) %>%
  ungroup()
hist(mamdat$sum.traps)
max(mamdat$sum.traps) # 1102 seems too high
#View(subset(mamdat, plot_date=="HARV_032 2014-09-19"))
#View(subset(mamdat, siteID=="HARV"))
# well, leaving it for now. Maybe they did a huge effort here for a side project
# *** issue: when trap code = 4, meaning that there were >1 smam in a single trap,
#         there are multiple rows of data for a single trap (1 per smam individual)
#         so # rows per plot_night is NOT entirely accurate estimate of # of traps set

# note: for now, going to ignore CPUE until I can deal with the above issue


#### add total smam abundance column ####
# add column for "capture"
mamdat$capture <- "0"
mamdat$capture[mamdat$trapStatus=="4 - more than 1 capture in one trap"] <- "1"
mamdat$capture[mamdat$trapStatus=="5 - capture"] <- "1"
mamdat$capture <- as.numeric(mamdat$capture)
mamdat <- mamdat %>% 
  group_by(plot_date) %>%
  mutate(sum_captures = sum(capture)) %>%
  ungroup()

#### aggregate by site, plot, year, and month ####
mamdat %>%
  group_by(siteID, plotID, Year, Month) %>%
  summarise(nlcdClass = first(nlcdClass),
            lat = first(decimalLatitude),
            lon = first(decimalLongitude),
            elevation = first(elevation),
            ntraps = sum(trapset),
            ncaps = sum(capture),
            richness = length(unique(taxonID))) %>%
  ungroup() -> mamdat_agg

#### aggregate to site ####
mamdat %>%
  group_by(siteID) %>%
  summarise(nlcdClass = first(nlcdClass),
            lat = first(decimalLatitude),
            lon = first(decimalLongitude),
            elevation = first(elevation),
            ntraps = sum(trapset),
            ncaps = sum(capture),
            richness = length(unique(taxonID))) %>%
  ungroup() -> mamdat_site_agg

#### PLOTS ####

# plot trapping *ROUGH* effort by date, coloured by site
ggplot(data=mamdat_agg, aes(x=Month, y=ntraps, group=siteID)) +
  geom_jitter(aes(colour=siteID), alpha=0.3, width=.2) +
  theme(legend.position="none") +
  facet_wrap(~Year) +
  scale_x_continuous(breaks=seq(0, 12, by=2))

# plot captures by coloured by site
ggplot(data=mamdat_agg, aes(x=Month, y=ncaps, group=siteID)) +
  geom_jitter(aes(colour=siteID), alpha=0.3, width=.2) +
  theme(legend.position="none") +
  facet_wrap(~Year) +
  scale_x_continuous(breaks=seq(0, 12, by=2))

# plot *ROUGH* CPUE by month x year
ggplot(data=mamdat_agg, aes(x=Month, y=ncaps/ntraps, group=siteID)) +
  geom_jitter(aes(colour=siteID), alpha=0.3, width=.2) +
  theme(legend.position="none") +
  facet_wrap(~Year) +
  scale_x_continuous(breaks=seq(0, 12, by=2))

# plot richness, month x year
ggplot(data=mamdat_agg, aes(x=Month, y=richness, group=siteID)) +
  geom_jitter(aes(colour=siteID), alpha=0.3, width=.2) +
  theme(legend.position="none") +
  facet_wrap(~Year) +
  scale_x_continuous(breaks=seq(0, 12, by=2))

#### MAPS ####
library(ggalt)  # for custom map projections

# Get the shape of North America
usa <- map_data("world", region = c("USA"))
# Exclude Hawaii if you want to
usa <- usa[!(usa$subregion %in% "Hawaii"),]
usa <- usa[!(usa$subregion %in% "Alaska"),]

# plot mean? effort by site
ggplot() +
  geom_map(map=usa, data=usa,
           aes(long, lat, map_id=region),
           colour="gray80", fill="gray80") +
  geom_point(data=mamdat_site_agg,
             aes(x=lon, y=lat, fill=richness),
             alpha=0.5, size=4, colour="grey30",
             shape=21) +
  scale_fill_viridis_c(option = "magma", direction = -1, begin = 0.2)



