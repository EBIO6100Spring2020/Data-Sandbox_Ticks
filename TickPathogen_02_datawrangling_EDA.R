##### Exploring ######

# Packages
library("tidyverse")
library("lubridate") # for dates
library("stringr") # for adding padded zeros
# library('gtable') if we want to do nested facets
# library('grid') if we want to do nested facets

# First, download all RAW files
# Tick field data
tck_field <- read.csv("data_raw/filesToStack10093/stackedFiles/tck_fielddata.csv", as.is = TRUE)
# Tick taxonomy data
tck_tax <- read.csv("data_raw/filesToStack10093/stackedFiles/tck_taxonomyProcessed.csv", as.is = TRUE)
# Tick pathogen data
tck_path <- read.csv("data_raw/filesToStack10092/stackedFiles/tck_pathogen.csv", as.is = TRUE)

dir.create("./output/EDA_plots")


###### Tick Field Data ########
# Check nested-ness of domain, site, and plot

dir.create("./output/EDA_plots/data_summary")
cat("SUMMARY INFORMATION ",""
    , file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", sep = "\n")
# Number of domains:
length(unique(tck_field$domainID))
cat("","NUMBER OF DOMAINS: ", length(unique(tck_field$domainID))
    , file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", sep = "\n", append=TRUE)
# Number of sites:
length(unique(tck_field$siteID))
cat("","NUMBER OF SITES: ", length(unique(tck_field$siteID))
    , file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", sep = "\n", append=TRUE)
# Number of plots
length(unique(tck_field$plotID))
cat("","NUMBER OF PLOTS: ", length(unique(tck_field$plotID))
    , file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", sep = "\n", append=TRUE)
# Number of sampling days
length(unique(as.Date(tck_field$collectDate)))
cat("","NUMBER OF SAMPLING DAYS: ", length(unique(as.Date(tck_field$collectDate)))
    , file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", sep = "\n", append=TRUE)
# Number of sampling years
length(unique(year(as.Date(tck_field$collectDate))))
cat("","SAMPLING YEARS: ", unique(year(as.Date(tck_field$collectDate)))
    , file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", sep = "\n", append=TRUE)

# Table of domain:site
tck_field %>%
  mutate(domainID=factor(domainID)) %>%
  mutate(siteID= factor(siteID, levels=unique(siteID))) %>%
  select(domainID, siteID) %>%
  table() %>%
  write.table(file="./output/EDA_plots/data_summary/tick_field_domain_site_matrix.txt", sep = "\t", quote = FALSE)

# Table of site:plot
tck_field %>%
  select(plotID, siteID) %>%
  table() %>%
  write.table(file="./output/EDA_plots/data_summary/tick_field_plot_site_matrix.txt", sep = "\t", quote = FALSE)

## Looking at coordinate and elevation uncertainty
ggsave(file="./output/EDA_plots/data_summary/tick_field_coordinateUncertainty.pdf",
  tck_field %>%
    ggplot() +
    geom_histogram(aes(x=coordinateUncertainty))
)
ggsave(file="./output/EDA_plots/data_summary/tick_field_elevationUncertainty.pdf",
tck_field %>%
  ggplot() +
  geom_histogram(aes(x=elevationUncertainty))
)
# I think we can ignore this uncertainty for our model for now.

## Looking at distribution of sampling areas
ggsave(file="./output/EDA_plots/data_summary/tick_field_totalSampledArea.pdf",
tck_field %>%
  ggplot() +
  geom_histogram(aes(x=totalSampledArea))
)

## What different types of plots are there?
cat("","TYPES OF PLOTS:", file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", sep="\n", append=TRUE)
tck_field %>%
  select(plotType) %>%
  table() %>%
  capture.output(file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", append=TRUE)
# They are all distributed; can get rid of this

## What different types of environments are there?
cat("","TYPES OF ENVIRONMENTS:", file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", sep="\n", append=TRUE)
tck_field %>%
  select(nlcdClass) %>%
  table() %>%
  capture.output(file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", append=TRUE)

# This is potentially important; let's keep this.

## What sampling conditions are there?
tck_field %>%
  select(sampleCondition) %>%
  table()
# Cold chain broken?
tck_field %>%
  filter(sampleCondition == "Cold chain broken (conditions described in remarks)") %>%
  select(remarks) %>%
  table()
# Distance estimated?
tck_field %>%
  filter(sampleCondition == "Distance estimated (described in remarks)") %>%
  select(remarks) %>%
  table()
# Other?
tck_field %>%
  filter(sampleCondition == "Other (described in remarks)") %>%
  select(remarks) 
# Tick count estimated?
tck_field %>%
  filter(sampleCondition == "Tick count estimated (error and instar specified in remarks)") %>%
  select(remarks) 
# Only thing noteworthy is "Larval count >500; check to see what actual larval count was"
tck_field %>%
  filter(remarks == "Larval count >500. Individuals above 500 were discarded.")
# In this sample, it looks likethey stopped counting after 500.
# This suck sbecause there are other samples where larval count was >500 anyway
tck_field %>%
  filter(larvaCount>500) %>%
  nrow()
# Not sure why they did this-- but I'm going to leave it at 500 for now and hope that they didn't discard TOO many.
# Final verdict: keep all samples.

## What is event ID?
tck_field %>%
  select(eventID) %>%
  table()
# I don't think we'll need this.

## Range in samplingMethods?
cat("","RANGE IN SAMPLING METHODS:", file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", append=TRUE, sep="\n")
tck_field %>%
  select(samplingMethod) %>%
  table() %>%
  capture.output(file="./output/EDA_plots/data_summary/tick_field_summary_info.txt", append=TRUE)
# Include this column.

# Lastly, for some reason sampleID is an incomplete table. I am going to fill this in.
# First, I check that I'm doing it correctly by cross-checking with present sampleIDs
tck_field %>%
  mutate(year= year(as.Date(collectDate))
         , month=str_pad(month(as.Date(collectDate)),2,side="left",pad="0")
         , day=str_pad(day(as.Date(collectDate)),2,side="left",pad="0")
         , spacer="." ) %>%
  unite(plotID,spacer,year,month,day, col = sampleID2, sep="") %>%
  filter(sampleID!="") %>%
  mutate(isSame=ifelse(sampleID==sampleID2, "same", "different")) %>%
  filter(isSame=="different")
# Strangely, there are some sampleIDs that were inputted that DIFFER from collectDate.
# Let me check the sampleIDs of the equivilent samples in the tick_tax set.
subset_tck_field <- tck_field %>%
  filter(sampleID !="")
tck_tax %>%
  filter(plotID %in% subset_tck_field$plotID, collectDate %in% subset_tck_field$collectDate) %>%
  select(plotID, collectDate, sampleID) %>%
  mutate(year= year(as.Date(collectDate))
         , month=str_pad(month(as.Date(collectDate)),2,side="left",pad="0")
         , day=str_pad(day(as.Date(collectDate)),2,side="left",pad="0")
         , spacer="." ) %>%
  unite(plotID,spacer,year,month,day, col = sampleID2, sep="") %>%
  mutate(isSame=ifelse(sampleID==sampleID2, "same", "different")) %>%
  filter(isSame=="different")
# There are also different sampleIDs here! Problem is, tck_tax has less samples than tck_field, so you can't do a simple match-subset.
# To keep things consistent between datasets, I think what I need to do is:
# (1) Match all sampleIDs between tck_tax and tck_field--> input these into tck_field
# (2) Then, manually create sampleIDs for tck_field given the "formula" it's supposed to follow.

# First, get sampleIDs from tck_path
samples_tck_tax <- tck_tax %>%
  select(plotID, collectDate, sampleID) %>%
  mutate(year= year(as.Date(collectDate))
         , month=str_pad(month(as.Date(collectDate)),2,side="left",pad="0")
         , day=str_pad(day(as.Date(collectDate)),2,side="left",pad="0")
         , spacer="." ) %>%
  unite(plotID,spacer,year,month,day, col = sampleID_2, sep="")

# Select only unecessary columns; adjust data as needed
tck_field_filt <- tck_field %>% 
  mutate(year= year(as.Date(collectDate))
         , month=str_pad(month(as.Date(collectDate)),2,side="left",pad="0")
         , day=str_pad(day(as.Date(collectDate)),2,side="left",pad="0")
         , spacer="." ) %>%
  unite(plotID,spacer,year,month,day, col = sampleID2, sep="", remove=FALSE) %>%
  mutate(sampleID_tax = samples_tck_tax[match(sampleID2, samples_tck_tax$sampleID_2), "sampleID"]) %>%
  mutate(sampleID = ifelse(sampleID !="", sampleID, ifelse(!is.na(sampleID_tax), sampleID_tax, sampleID2))) %>% ### Combines all sampleID options in a hierarchy
  select(uid, domainID, siteID, plotID, sampleID
         , nlcdClass
         , decimalLatitude, decimalLongitude, elevation
         , collectDate
         , samplingMethod, totalSampledArea, adultCount, nymphCount, larvaCount
         , samplingProtocolVersion, measuredBy)

##### Tick taxonomy ######

cat("SUMMARY INFO","", file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", sep="\n")
## What conditions are there?
cat("","SAMPLE CONDITIONS:", file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", sep="\n", append=TRUE)
tck_tax %>% 
  select(sampleCondition) %>%
  table() %>%
  capture.output(file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", append=TRUE)
# Just get rid of all non-OKs

## Summary of sex and age
cat("","SEX AND AGE:", file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", sep="\n", append=TRUE)
tck_tax %>% 
  select(sexOrAge) %>%
  table() %>%
  capture.output(file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", append=TRUE)

## What's up with sex and age?
tck_tax %>%
  select(sexOrAge) %>%
  table()
# We could get rid of the 4 that aren't classified properly... but let's keep them for now.
# Are the 3 unnamed ones adults, or truly don't know?
tck_tax %>%
  filter(sexOrAge=="") %>% select(remarks)
# Okay, so we get rid of those 3 because they are not ticks.
# But let's keep the unsexed adult and change it to appropriate columns.

## Check if identification protocol version is different
tck_tax %>%
  select(identificationProtocolVersion) %>%
  table()
# it's all the same

## What are the individual counts?
tck_tax %>%
  ggplot() +
  geom_histogram(aes(individualCount))
# Extremely long right tail; mode at 1

## Let's look at taxonomy
cat("","FAMILY AND GENUS:", file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", sep="\n", append=TRUE)
tck_tax %>%
  select(family, genus) %>%
  table() %>%
  capture.output(file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", append=TRUE)
# Do the individuals with no famiy HAVE no family, or are they simply not added?
tck_tax %>%
  filter(family=="")  %>%
  select(scientificName, acceptedTaxonID, taxonRank, family, subfamily, tribe, subtribe, genus, subgenus, specificEpithet, infraspecificEpithet)
# Ah, I see-- they are only ID'd down to order because I'm assuming larva are difficult to ID.
# Are there any that don't belong to the broader order of Ixodida?
tck_tax %>%
  filter(family=="") %>%
  select(scientificName) %>% table()
# Three; but these are the "not a tick" ones.
tck_tax %>%
  filter(scientificName=="")

## How does the rest of the data look, down to genus?
cat("","GENUS AND SP", file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", sep="\n", append=TRUE)
tck_tax %>%
  select(genus, specificEpithet) %>%
  table() %>% 
  capture.output(file="./output/EDA_plots/data_summary/tick_taxonomy_data_summary_info.txt", append=TRUE)
# Looks good

tck_tax %>%
  select(infraspecificEpithet) %>%
  table()
tck_tax %>%
  select(subgenus, subfamily) %>%
  table()
tck_tax %>%
  select(tribe, subtribe) %>%
  table()
tck_tax %>%
  select(family) %>%
  table()
# We can just get rid of infraspecificEpithet, tribe, subtribe, subfamily and subgenus since they are just empty columns.
# We should probably standardize the ID's a bit so the individuals have a proper nested taxonomy.
# Are all of the other genera also in the family Ixodida?--> Google searched this; yes, they are all in the ORDER Ixodida and FAMILY Ixodidae

## Lastly, check the different types of taxon tanks
tck_tax %>%
  select(taxonRank) %>%
  table()
# Looks like this column is complete; the 3 unamed are the 3 non-ticks.
# It seems there is a group of ticks that are in the family Ixodidae, and another group (141) that can't be ID'd past order
# Therefore, we should probably also make a separate ORDER column.
tck_tax %>%
  filter(taxonRank=="order") 
tck_tax %>%
  filter(taxonRank=="family")

# Below, I create new sampleIDs and subsampleIDs so they are consistent between all datasets.
tck_tax_filt <- tck_tax %>%
  filter(sampleCondition=="OK", sexOrAge!="") %>%
  mutate(lifeStage=ifelse(sexOrAge %in% c("Male","Female"), "Adult", sexOrAge)
         ,sex=ifelse(sexOrAge %in% c("Male","Female"), sexOrAge, NA) ) %>%
  mutate(order="Ixodida"
         , family=ifelse(family=="", NA, family)
         , genus = ifelse(genus == "", NA, genus)
         , specificEpithet = ifelse(specificEpithet == "",NA,specificEpithet)) %>%
  select(plotID, sampleID, subsampleID
         , collectDate, identifiedDate
         , acceptedTaxonID, scientificName, order, family, genus, specificEpithet, taxonRank
         , lifeStage, sex, individualCount
  )

  #### Note: what I don't understand is the difference between sp. and spp. In this scenario, you either do or don't know the ID, so why ever use spp.?



###### Tick Pathogen Data #######
#### Looking at tck pathogen 
### Looking at tck taxonomy
cat("SUMMARY FILE","", file="./output/EDA_plots/data_summary/tick_pathogen_data_summary_info.txt", sep="\n")
# All plots an samples are found within tck_field? If false, then all included.
any(!tck_tax$sampleID %in% tck_field$sampleID)

# Summary of IDs against plots
tck_tax %>%
  select(plotID, acceptedTaxonID) %>%
  table()

## How many tests, and what results?
cat("","TEST RESULTS:", file="./output/EDA_plots/data_summary/tick_pathogen_data_summary_info.txt", sep="\n", append=TRUE)
table(tck_path$testResult) %>%
  capture.output( file="./output/EDA_plots/data_summary/tick_pathogen_data_summary_info.txt", append=TRUE)
# We should get rid of all rows with unknown

## Are all counts just 1?
table(tck_path$individualCount)
# Yes; we can get rid of this column since each row is an individual observation

## What types of conditions are there?
cat("","SAMPLING CONDITION:", file="./output/EDA_plots/data_summary/tick_pathogen_data_summary_info.txt", sep="\n", append=TRUE)
table(tck_path$sampleCondition) %>%
  capture.output(file="./output/EDA_plots/data_summary/tick_pathogen_data_summary_info.txt", append=TRUE)
# Most are OK, so let's get rid of everything else.

## How many diff types of pathogens?
cat("","DIFFERENT TYPES OF PATHOGENS", file="./output/EDA_plots/data_summary/tick_pathogen_data_summary_info.txt", sep="\n", append=TRUE)
table(tck_path$testPathogenName, tck_path$testResult) %>%
  capture.output(file="./output/EDA_plots/data_summary/tick_pathogen_data_summary_info.txt", append=TRUE)
unique(tck_path$testPathogenName)
# What is Borrelia sp. and HardTick DNA Quality?

# First, is HardTick DNA Quality supposed to be a separate column from all other data?
tck_path %>%
  filter(testPathogenName != "HardTick DNA Quality") %>%
  select(subsampleID) %>%
  unique()

# Check if there are any subsamples that are NOT in HardTick DNA quality subsample list
# If false, probably means we can get rid of HardTick DNA Quality rows
any(!(tck_path %>%
  filter(testPathogenName == "HardTick DNA Quality") %>%
  pull(subsampleID)) %in% (tck_path %>%
  filter(testPathogenName != "HardTick DNA Quality") %>%
  pull(subsampleID) %>%
  unique()))
# Okay, so we can get rid of HardTick DNA Quality.

# Check plotType to see if all the same
tck_path %>%
  select(plotType) %>%
  table()
# All same; can get rid of this.

# Check batch ID nad testProtocolVersion
tck_path %>%
  select(batchID, testProtocolVersion) %>%
  table()
# Only one protocol version; get rid of this column.

# Check batch ID with testedBy
tck_path %>%
  arrange(testedBy) %>%
  mutate(batchID=factor(batchID, levels=unique(batchID))) %>%
  select(batchID, testedBy) %>%
  table()
# batchID is nested with testedBy; so keep both of these random effects.

# Remove all results that weren't completed. 
# Filter to relevant columns of data
# change "age or sex" into proper columns.
tck_path_filt <- tck_path %>%
  filter(testResult %in% c("Positive","Negative"), sampleCondition %in% c("OK"), testPathogenName != "HardTick DNA Quality") %>%
  select(uid, domainID, siteID, plotID, subsampleID # nestedness
         , nlcdClass # Enviro
         , decimalLatitude, decimalLongitude, elevation # geographic
         , collectDate, testedDate # time
         , testResult, testPathogenName # test results
         , batchID, testedBy # random effects
         ) 

###### Plotting pathogen data to visualize variables #######

# Now, let's see distribution of pathogens across plots and years
dir.create("./output/EDA_plots/tick_pathogen")

ggsave(filename="./output/EDA_plots/tick_pathogen/pathogen_vs_plotID.pdf", width=12, height=8,
tck_path_filt %>%
  arrange(domainID, siteID) %>%
  mutate(plotID=factor(plotID, levels=unique(plotID))) %>%
  mutate(Year=year(as.Date(collectDate))) %>%
  group_by(plotID, Year, testPathogenName) %>%
  summarise(ProportionDetected = sum(testResult=="Positive")/length(testResult), TotalTests = length(testResult)) %>%
  ggplot() +
  geom_point(aes(x=plotID, y=testPathogenName, col=ProportionDetected, cex=log(TotalTests))) +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Pathogen Name") +
  facet_wrap(.~Year, nrow=4)
)

dir.create("./output/EDA_plots/tick_pathogen/byPathogen")
# Get ordered plotID by domain and site
plot_order <- tck_path_filt %>%
  arrange(domainID, siteID) %>%
  pull(plotID) %>%
  unique()
for ( p in unique(tck_path_filt$testPathogenName)) {
  ggsave(filename = paste0("./output/EDA_plots/tick_pathogen/byPathogen/",p,".pdf"), width=12, height=8,
  tck_path_filt %>%
    filter(testPathogenName==p) %>%
    mutate(Year=year(as.Date(collectDate)), yDay = yday(as.Date(collectDate))) %>%
    group_by(plotID, Year, yDay, testPathogenName) %>%
    summarise(ProportionDetected = sum(testResult=="Positive")/length(testResult), TotalTests = length(testResult)) %>%
    ungroup() %>%
    mutate(ProportionDetected = ifelse(ProportionDetected>0,as.numeric(ProportionDetected),NA)) %>%
    mutate(ProportionDetected = as.numeric(ProportionDetected), plotID = factor(plotID, levels=plot_order)) %>%
    ggplot() +
    geom_point(aes(x=yDay, y=plotID, col=(ProportionDetected), cex=log(TotalTests))) +
    theme(axis.text.x = element_text(angle=90)) +
    ylab("PlotID") + xlab("Day of year") +
    scale_color_gradientn(limits = c(0,1),
                         colours=c("white", "darkred"),
                         breaks=c(0,1), labels=format(c(0,1)))+
    xlim(c(0,365))+
    scale_y_discrete(drop=FALSE)
  )
}
# 
# ggsave(filename="./output/EDA_plots/allPathogens_plotvsdate.pdf",
#   tck_path_filt %>%
#   arrange(domainID, siteID) %>%
#   mutate(plotID=factor(plotID, levels=unique(plotID))) %>%
#   # filter(testPathogenName==p) %>%
#   mutate(Year=year(as.Date(collectDate)), yDay = yday(as.Date(collectDate))) %>%
#   group_by(plotID, Year, yDay, testPathogenName) %>%
#   summarise(ProportionDetected = sum(testResult=="Positive")/length(testResult), TotalTests = length(testResult)) %>%
#     ungroup() %>%
#     mutate(ProportionDetected = ifelse(ProportionDetected>0,as.numeric(ProportionDetected),NA)) %>%
#     mutate(ProportionDetected = as.numeric(ProportionDetected), plotID = factor(plotID, levels=plot_order)) %>%
#     ggplot() +
#   geom_point(aes(x=yDay, y=plotID, col=as.numeric(ProportionDetected), cex=log(TotalTests))) +
#   theme(axis.text.x = element_text(angle=90)) +
#   ylab("PlotID") + xlab("Day of year") +
#   scale_color_gradientn(limits = c(0,1),
#                         colours=c("white", "darkred"),
#                         breaks=c(0,1), labels=format(c(0,1)))+
#   xlim(c(0,365)) +
#   facet_grid(testPathogenName~.)
#   , heigh=40, width=12
# )

dir.create("./output/EDA_plots/tick_pathogen/byPlot")
for ( p in unique(tck_path_filt$plotID)) {
  ggsave(filename=paste0("./output/EDA_plots/tick_pathogen/byPlot/",p,".pdf"), width=12, height=8,
  tck_path_filt %>%
    filter(plotID==p) %>%
    mutate(Year=year(as.Date(collectDate)), yDay = yday(as.Date(collectDate))) %>%
    group_by(plotID, Year, yDay, testPathogenName) %>%
    summarise(ProportionDetected = sum(testResult=="Positive")/length(testResult), TotalTests = length(testResult)) %>%
    ungroup() %>%
    mutate(ProportionDetected = ifelse(ProportionDetected>0,ProportionDetected,NA)) %>%
    mutate(ProportionDetected = as.numeric(ProportionDetected), plotID = factor(plotID, levels=plot_order)) %>%
    ggplot() +
    geom_point(aes(x=yDay, y=testPathogenName, col=ProportionDetected, cex=log(TotalTests))) +
    theme(axis.text.x = element_text(angle=90)) +
    ylab("PathogenName") + xlab("Day of year") +
    scale_color_gradientn(limits = c(0,1),
                          colours=c("white", "darkred"),
                          breaks=c(0,1), labels=format(c(0,1)))+
    xlim(c(0,365)) 
)  
}

dir.create("./output/EDA_plots/tick_pathogen/byDomain_pathogenvsdate")
for ( d in unique(tck_path_filt$domainID)) {
  nplots <- tck_path_filt %>%
    filter(domainID==d) %>%
    pull(plotID) %>%unique() %>%length()
  ggsave(filename=paste0("./output/EDA_plots/tick_pathogen/byDomain_pathogenvsdate/",d,".pdf"), width=12, height=nplots*2+1,
         tck_path_filt %>%
           filter(domainID==d) %>%
           mutate(Year=year(as.Date(collectDate)), yDay = yday(as.Date(collectDate))) %>%
           group_by(siteID, plotID, Year, yDay, testPathogenName) %>%
           summarise(ProportionDetected = sum(testResult=="Positive")/length(testResult), TotalTests = length(testResult)) %>%
           ungroup() %>%
           mutate(ProportionDetected = ifelse(ProportionDetected>0,ProportionDetected,NA)) %>%
           mutate(ProportionDetected = as.numeric(ProportionDetected), plotID = factor(plotID, levels=plot_order)) %>%
           ggplot() +
           geom_point(aes(x=yDay, y=testPathogenName, col=ProportionDetected, cex=log(TotalTests))) +
           theme(axis.text.x = element_text(angle=90)) +
           ylab("PathogenName") + xlab("Day of year") +
           scale_color_gradientn(limits = c(0,1),
                                 colours=c("white", "darkred"),
                                 breaks=c(0,1), labels=format(c(0,1)))+
           xlim(c(0,365)) +
           facet_grid(siteID+ plotID ~ .)
  )
}


#### Merging taxonomy with tick field #### 



#### Merging with tick pathogens #####


# 
# 
# # Combined plot of domain:site:plot
# for ( d in unique(tck_field$domainID)) {
#     d=unique(tck_field$domainID)[1]
#     temp <- tck_field %>%
#       as_tibble() %>%
#       filter(domainID==d)
#     for ( s in unique(temp$siteID)) {
#       s= unique(temp$siteID)[2]
#       temp2 <- temp %>%
#         filter(siteID==s)
#       for ( p in unique(temp$plotID)) {
#         temp2 %>%
#           filter(plotID==p) %>%
#           mutate(Date=as.Date(collectDate)) %>%
#           mutate(Year=year(Date), Month=month(Date), Day=yday(Date)) %>%
#           # mutate(siteID=as.character(siteID)) %>%
#           select(Date, Year, Month, Day, adultCount, nymphCount, larvaCount) %>%
#           gather(adultCount, nymphCount, larvaCount, key=lifestage, value=count) %>%
#           complete() %>%
#           as_tibble() %>%
#           ggplot() +
#           geom_point(aes(x=lifestage, y=count, col=factor(Year)))
#         # facet_grid(.~ siteID)
#       }
#     }
# }


