#!/bin/bash

##### This is the official script for filtering tick field data. Created 20feb2020
library("tidyverse") # data wrangling
library("lubridate") # for dates
library("stringr") # for adding padded zeros

# First, download all RAW files
# Tick field data
tck_field <- read.csv("data_raw/filesToStack10093/stackedFiles/tck_fielddata.csv", as.is = TRUE)
# Tick taxonomy data
tck_tax <- read.csv("data_raw/filesToStack10093/stackedFiles/tck_taxonomyProcessed.csv", as.is = TRUE)
# Tick pathogen data
tck_path <- read.csv("data_raw/filesToStack10092/stackedFiles/tck_pathogen.csv", as.is = TRUE)


###### Tick field data manipulations #######
##### 1. Make useful dates #####

tck_field_1 <- tck_field %>% 
  mutate(year= year(as.Date(collectDate))
         , month=month(as.Date(collectDate))
         , day=day(as.Date(collectDate))
         , dayOfYear=yday(as.Date(collectDate)))

tck_tax_1 <- tck_tax %>%
  mutate(year= year(as.Date(collectDate))
         , month=month(as.Date(collectDate))
         , day=day(as.Date(collectDate))
         , dayOfYear=yday(as.Date(collectDate)))

tck_path_1 <- tck_path %>%
  mutate(year= year(as.Date(collectDate))
         , month=month(as.Date(collectDate))
         , day=day(as.Date(collectDate))
         , dayOfYear=yday(as.Date(collectDate)))

##### 2. SampleID problem ######
# First, get sampleIDs from tck_path
samples_tck_tax <- tck_tax_1 %>%
  group_by(plotID, collectDate, sampleID) %>%
  summarize(check=1) %>%
  select(plotID, collectDate, sampleID) %>%
  rename(sampleID_tax = sampleID)

# Adjust sample IDs
tck_field_2 <- tck_field_1 %>% 
    rename(sampleID_original = sampleID) %>%
    mutate(month_pad=str_pad(month,2,side="left",pad="0")
         , day_pad=str_pad(day,2,side="left",pad="0")
         , spacer="." ) %>%  
    unite(plotID,spacer,year,month_pad,day_pad, col = sampleID2, sep="", remove=FALSE) %>%
  select(-month_pad, -day_pad) %>%
    left_join(samples_tck_tax) %>%
    mutate(sampleID = ifelse(sampleID_original!="", sampleID_original, ifelse(!is.na(sampleID_tax), sampleID_tax, sampleID2))) ### Combines all sampleID options in a hierarchy

##### 3. Remove samples that had flooding, fire, or other "acts of god" (or are incomplete) #####
# Looked at the "sampleCondition" and "remarks" sections to see if there was anything noteworthy.

tck_field_3 <- tck_field_2 %>%
  filter(!(sampleID %in% c("OSBS_022.20160310"
                           ,"LAJA_003.20161214"
                           ,"LAJA_001.20161214"
                           ,"LAJA_002.20161214"
                           ,"LAJA_030.20161215"
                           ,"UKFS_003.20170523"
                           ,"TALL_001.20180522"
                           ,"GRSM_009.20180704")))

##### 4. Investigate non-unique sampleIDs in tck_field #####
# length(unique(tck_field_3$sampleID))
# length(tck_field_3$sampleID)
# tck_field_3 %>%
#   mutate(duplicated=duplicated(sampleID)) %>%
#   filter(duplicated) %>%
#   select(sampleID)%>% pull()
# There are 16 non-unique samples

# Remove all SECOND drags; sort by time, then remove duplicates
tck_field_4 <- tck_field_3 %>%
  arrange(collectDate) %>%
  mutate(dup=duplicated(sampleID)) %>%
  filter(!dup)

# Check it is the same length
# length(unique(tck_field_4$sampleID))
# length(tck_field_4$sampleID)

##### 5. Edit and filter tck_tax, tck_field, tck_path #####
tck_tax_2 <- tck_tax_1 %>%
  filter(sampleCondition=="OK") %>% # get rid of compromised samples (3 of them)
  filter(sexOrAge!="") %>% # Get rid of non-ticks
  mutate(lifeStage=ifelse(sexOrAge %in% c("Male","Female"), "Adult", sexOrAge)
         ,sex=ifelse(sexOrAge %in% c("Male","Female"), sexOrAge, NA) ) %>%
  mutate(order="Ixodida"
         , family=ifelse(family=="", NA, family)
         , genus=ifelse(genus == "", NA, genus)
         , specificEpithet = ifelse(specificEpithet == "",NA,specificEpithet)) %>%
  select(plotID, sampleID, subsampleID
         , collectDate
         , acceptedTaxonID, scientificName, order, family, genus, specificEpithet, taxonRank
         , lifeStage, sex, individualCount
  ) 

tck_field_5 <- tck_field_4 %>%
  filter(!(is.na(adultCount)& is.na(nymphCount)& is.na(larvaCount))) %>% # get rid of all samples with NO tick counts
  mutate(larvaCount = ifelse(is.na(larvaCount),0,larvaCount)) %>%
  select(domainID, siteID, plotID, sampleID # sample identification info
         , nlcdClass, decimalLatitude, decimalLongitude, elevation # Environmental conditions
         , collectDate, year, month, day, dayOfYear
         , totalSampledArea, adultCount, nymphCount, larvaCount
         )

tck_path_2 <- tck_path_1 %>%
  filter(sampleCondition=="OK") %>%
  filter(testPathogenName !="HardTick DNA Quality") %>%
  separate(subsampleID, into=c("site","plot","date","acceptedTaxonID","lifeStage"), remove=FALSE) %>%
  unite(site, plot, sep="_", col="plotID_temp") %>%
  unite(plotID_temp, date, sep=".", col=sampleID) %>%
  select(domainID, siteID, plotID, subsampleID, testingID# sample identification info
         , nlcdClass, decimalLatitude, decimalLongitude, elevation # Environmental conditions
         , collectDate, year, month, day, dayOfYear
         , individualCount, testPathogenName, testResult
  )

##### 6. Matching counts #####

# Okay, make a tck_tax that has JUST counts
MATCHEDCOUNTS_tckfield_tcktax <- tck_tax_2 %>%
  mutate(sexOrAge=ifelse(is.na(sex),lifeStage,sex)) %>%
  select(sampleID, sexOrAge, individualCount, acceptedTaxonID) %>%
  # spread(key=ageOrSex, value=individualCount)
  group_by(sampleID, sexOrAge) %>% # make male and female adults the same line
  summarize(allCounts=sum(individualCount)) %>%
  spread(key=sexOrAge, value=allCounts) %>%
  mutate(adultCount_tax=sum(c(Female, Male), na.rm = TRUE), nymphCount_tax = Nymph, larvaCount_tax = Larva) %>%
  select(sampleID, adultCount_tax, Female, Male, nymphCount_tax, larvaCount_tax)%>%
  right_join(tck_field_5) %>%
  mutate(totalCount_nolarva_tax=sum(c(adultCount_tax, nymphCount_tax), na.rm = TRUE)
         , totalCount_nolarva = sum(c(adultCount, nymphCount), na.rm=TRUE)) %>%
  select(sampleID, totalCount_nolarva, totalCount_nolarva_tax, adultCount, adultCount_tax, Female, Male, nymphCount, nymphCount_tax, larvaCount, larvaCount_tax) %>%
  replace(is.na(.), 0) 

#*********** Error level 0 ************#
# Perfect match; keep
# Get sampleIDs with identical counts.
toKeep <- MATCHEDCOUNTS_tckfield_tcktax %>%
  filter(adultCount == adultCount_tax, nymphCount == nymphCount_tax, larvaCount==larvaCount_tax) %>%
  pull(sampleID) 

tck_field_adj <- tck_field_5 %>%
  filter(sampleID%in%toKeep) 

#*********** Error level 1 ************#
# Larva counts wrong only; trust tck_field.
toKeep1 <- MATCHEDCOUNTS_tckfield_tcktax %>%
  filter(adultCount == adultCount_tax & nymphCount == nymphCount_tax & larvaCount!=larvaCount_tax) %>%
  pull(sampleID) 

tck_field_adj <- tck_field_5 %>%
  filter(sampleID %in% toKeep1) %>%
  full_join(tck_field_adj)

#*********** Error level 2 ************#
# No taxonomy at all; trust field counts.
toKeep2 <- MATCHEDCOUNTS_tckfield_tcktax %>%
  filter(adultCount != adultCount_tax | nymphCount != nymphCount_tax) %>%
  filter(totalCount_nolarva_tax==0) %>%
  pull(sampleID)

tck_field_adj <- tck_field_5 %>%
  filter(sampleID %in% toKeep2) %>%
  full_join(tck_field_adj)

#*********** Error level 3 ************#
# Numbers got shuffled between tick IDs; keep taxonomy counts.
# Get adjusted sample IDs
adj_counts_3 <- MATCHEDCOUNTS_tckfield_tcktax %>%
  filter(adultCount != adultCount_tax | nymphCount != nymphCount_tax) %>%
  filter(totalCount_nolarva_tax!=0) %>%
  mutate(totalCount = adultCount + nymphCount + larvaCount, totalCount_tax = adultCount_tax+ nymphCount_tax + larvaCount_tax) %>%
  filter(totalCount == totalCount_tax) %>%
  select(sampleID, adultCount_tax, nymphCount_tax, larvaCount_tax) 

toAdj3 <- adj_counts_3 %>%
  pull(sampleID)

tck_field_adj <- tck_field_5 %>%
  filter(sampleID %in% toAdj3) %>%
  left_join(adj_counts_3) %>%
  select(-c(adultCount, nymphCount, larvaCount)) %>%
  rename(adultCount=adultCount_tax, nymphCount=nymphCount_tax, larvaCount=larvaCount_tax) %>%
  full_join(tck_field_adj)

#*********** Error level 4 ************#
# Field > tax, and adult/nymph variation is stable.
adj_counts_4 <- MATCHEDCOUNTS_tckfield_tcktax %>%
  filter(adultCount != adultCount_tax | nymphCount != nymphCount_tax) %>%
  filter(totalCount_nolarva_tax!=0) %>%
  mutate(totalCount = adultCount + nymphCount + larvaCount, totalCount_tax = adultCount_tax+ nymphCount_tax + larvaCount_tax) %>%
  filter(totalCount != totalCount_tax) %>%
  mutate(total_diff = totalCount-totalCount_tax, nonLarv_diff = totalCount_nolarva-totalCount_nolarva_tax) %>%
  filter(total_diff>0) %>%
  mutate(larvaCount_est = totalCount-adultCount_tax-nymphCount_tax)

# Quickly check that estimated larva counts are pretty close to field-observed larva counts
adj_counts_4 %>%
  ggplot() +geom_point(aes(x=larvaCount, y=larvaCount_est))
# There's only one point of concern-- manually add this one back in since it is pretty obvious they just didn't ID all 2800 ticks.
exception <- adj_counts_4 %>%
  filter((larvaCount_est>2500&larvaCount==0)) %>%
  select(sampleID, adultCount_tax, nymphCount_tax, larvaCount_est)
adj_counts_4 <- adj_counts_4 %>%
  filter(!(larvaCount_est>2500&larvaCount==0)) %>%
  select(sampleID, adultCount_tax, nymphCount_tax, larvaCount_est)

toAdj4 <- adj_counts_4 %>% pull(sampleID)
toAdj4_exception <- exception %>%pull(sampleID)

tck_field_adj <- tck_field_5 %>%
  filter(sampleID %in% toAdj4) %>%
  left_join(adj_counts_4) %>%
  select(-c(adultCount, nymphCount, larvaCount)) %>%
  rename(adultCount=adultCount_tax, nymphCount=nymphCount_tax, larvaCount=larvaCount_est) %>%
  full_join(tck_field_adj)
# Manually add that one back in; but keep field counts.
tck_field_adj <- tck_field_5 %>%
  filter(sampleID %in% toAdj4_exception) %>%
  left_join(exception) %>%
  select(-c(adultCount_tax, nymphCount_tax, larvaCount_est)) %>%
  full_join(tck_field_adj)

#*********** Error level 5 ************#
# What about tax>field?
adj_counts_5 <- MATCHEDCOUNTS_tckfield_tcktax %>%
  filter(adultCount != adultCount_tax | nymphCount != nymphCount_tax) %>%
  filter(totalCount_nolarva_tax!=0) %>%
  mutate(totalCount = adultCount + nymphCount + larvaCount, totalCount_tax = adultCount_tax+ nymphCount_tax + larvaCount_tax) %>%
  filter(totalCount != totalCount_tax) %>%
  mutate(total_diff = totalCount-totalCount_tax, nonLarv_diff = totalCount_nolarva-totalCount_nolarva_tax) %>%
  filter(total_diff<0, abs(nonLarv_diff) <= 10) 
# Note: a threshold of 10 was chosen because I looked at the tck_path file and found support for
# all of these 
# Quickly check that tax counts are pretty close to field-observed larva counts
adj_counts_5 %>%
  ggplot() +geom_point(aes(x=adultCount, y=adultCount_tax), col="darkred")+
  geom_point(aes(x=nymphCount, y=nymphCount_tax), col="purple", alpha=0.1)  +
  geom_point(aes(x=larvaCount, y=larvaCount_tax), col="green", alpha=0.1)
# Get rid of that one outlier point-- is very strange. After some sleuthing, it's sample "ORNL_007.20150618"
adj_counts_5 <- adj_counts_5 %>%
  mutate(adult_diff = adultCount-adultCount_tax) %>%
  filter(adult_diff<=18) 

adj_counts_5 <- adj_counts_5 %>%
  select(sampleID, adultCount_tax, nymphCount_tax, larvaCount_tax)

toAdj5 <- adj_counts_5 %>% pull(sampleID)

tck_field_adj <- tck_field_5 %>%
  filter(sampleID %in% toAdj5) %>%
  left_join(adj_counts_5) %>%
  select(-c(adultCount, nymphCount, larvaCount)) %>%
  rename(adultCount=adultCount_tax, nymphCount=nymphCount_tax, larvaCount=larvaCount_tax) %>%
  full_join(tck_field_adj)

#*********** Checking to see if counts are consistent now  ************#
# Now, we check to see if counts are consistent. 
tck_tax_2 %>%
  mutate(sexOrAge=ifelse(is.na(sex),lifeStage,sex)) %>%
  select(sampleID, sexOrAge, individualCount, acceptedTaxonID) %>%
  # spread(key=ageOrSex, value=individualCount)
  group_by(sampleID, sexOrAge) %>% # make male and female adults the same line
  summarize(allCounts=sum(individualCount))%>% 
  spread(key=sexOrAge, value=allCounts) %>%
  mutate(adultCount_tax=sum(c(Adult, Female, Male), na.rm = TRUE), nymphCount_tax = Nymph, larvaCount_tax = Larva) %>%
  select(sampleID, adultCount_tax, nymphCount_tax, larvaCount_tax) %>%
  right_join(tck_field_adj) %>%
  mutate(totalCount_nolarva_tax=sum(c(adultCount_tax, nymphCount_tax), na.rm = TRUE)
         , totalCount_nolarva = sum(c(adultCount, nymphCount), na.rm=TRUE)) %>%
  select(sampleID, totalCount_nolarva, totalCount_nolarva_tax, adultCount, adultCount_tax,  nymphCount, nymphCount_tax, larvaCount, larvaCount_tax) %>%
  mutate(diff=totalCount_nolarva!=totalCount_nolarva_tax) %>%
  filter(!is.na(adultCount_tax),!is.na(nymphCount_tax) ) %>%
  filter(sampleID != toAdj4_exception) %>% # There was that one exception where it won't match
  pull(diff) %>% any() 
# All adult/larva are the same now, EXCEPT when it is zero or NA in the tck_path dataset.

##### 7. Merging taxonomy and field counts #####

# Can uncount, but the problem is you also get rid of zeros, which can be important. 
# So, extract zeros first and THEN uncount.
temp_field_zeros <- tck_field_adj %>%
  gather(adultCount, nymphCount, larvaCount, key=lifeStageCat, value=Count) %>%
  filter(Count==0 | is.na(Count))

# Uncount
temp_field_long <- tck_field_adj %>%
  gather(adultCount, nymphCount, larvaCount, key=lifeStageCat, value=Count) %>%
  filter(Count!=0) %>%
  uncount(weights = Count, .id = "unique") %>%
  full_join(temp_field_zeros) %>%
  mutate(Count=ifelse(is.na(Count),1,Count)) %>%
  mutate(lifeStage=ifelse(lifeStageCat=="adultCount","Adult", ifelse(lifeStageCat=="nymphCount","Nymph",ifelse(lifeStageCat=="larvaCount","Larva",NA)))) %>%
  select(-lifeStageCat)

temp_tax_long <- tck_tax_2 %>%
  uncount(weights=individualCount) 
# Manually make unique ID; don't want unique by species or sex; only by life stage and sample.
temp_tax_long$unique <- with(temp_tax_long, ave(rep(1, nrow(temp_tax_long)), sampleID, lifeStage, FUN = seq_along))

### Finally, add them together!!!
tck_field_tax_joined <-  temp_field_long %>%
  left_join(temp_tax_long) %>%
  select(domainID, siteID, plotID, sampleID, subsampleID
         , nlcdClass, decimalLatitude, decimalLongitude, elevation
         , collectDate, year, month, day, dayOfYear
         , totalSampledArea
         , lifeStage, sex, taxonRank, acceptedTaxonID, order, family, genus, specificEpithet, Count)

##### 8. Merging tick counts/taxonomy with tick pathogen #####

# are all subsamples in tck_path also in other dataset?
any(!(tck_path_2$subsampleID %in% tck_field_tax_joined$subsampleID))
# Convert to short format for pathogen testing
tck_path_2_short <- tck_path_2 %>%
  mutate(testPathogenName = gsub(testPathogenName,pattern = " ", replacement = "_"))%>%
  spread(key=testPathogenName, value=testResult)  %>%
  select(-individualCount) 

# Merge with other data
tck_all_merged <- left_join(tck_field_tax_joined, tck_path_2_short) %>%
  replace(.=="",NA)

write_csv(tck_all_merged, path="data_derived/MASTER_all_tck_data_merged.csv")

##### 9. Creating subsets and appropriate summaries #####
dir.create("data_derived/subset_data")
dir.create("data_derived/subset_data/ticks_taxonomy_collapsed")
dir.create("data_derived/subset_data/ticks_taxonomy_separated")
dir.create("data_derived/subset_data/ticks_IXOSCA")
dir.create("data_derived/subset_data/ticks_AMBAME")


#### Dataset 1: Create dataset of ALL tick abundances, taxonomy collapsed and separated ####
all_ticks_taxonomy_collapsed <- tck_all_merged %>%
  group_by(domainID, siteID, plotID, sampleID
           , nlcdClass, decimalLatitude, decimalLongitude, elevation
           , collectDate, year, month, day, dayOfYear, totalSampledArea 
           , lifeStage) %>%
  summarize(Count=sum(Count)) 

# We can check how well the altered counts and non-altered counts compare to each other
tck_field_5 %>%
  gather(adultCount, nymphCount, larvaCount, key=lifeStage, value=Original_Count) %>%
  mutate(lifeStage=ifelse(lifeStage=="adultCount", "Adult",ifelse(lifeStage=="nymphCount","Nymph",ifelse(lifeStage=="larvaCount","Larva",NA)))) %>%
  select(sampleID, lifeStage, Original_Count) %>%
  right_join(all_ticks_taxonomy_collapsed) %>%
  ggplot() +geom_point(aes(x=log(Original_Count), y=log(Count), col=lifeStage), alpha=0.25) 

# The problem with grouping by sample and species then summing this is that sometimes the number of identified taxa are NOT equal to number of counted taxa.
# Therefore, we should perhaps use identifications as a proxy of PROPORTION of species? This is potentially important when there is 1 larva ID'd, but 500 larva in the sample.

# Count how many ID's were made for each sample and life stage
tck_countsIDed_only <- tck_all_merged %>%
  mutate(identified = !is.na(acceptedTaxonID)) %>%
  # filter(!is.na(acceptedTaxonID)) %>%
  group_by(sampleID, lifeStage, identified) %>%
  summarize(IDsMadeInSample=sum(identified)) %>%
  mutate(IDsMadeInSample = ifelse(IDsMadeInSample==0,NA,IDsMadeInSample))
# Now, standardize actual counts for each species by the proportion of ticks that were ID'd
all_ticks_taxonomy_separated_countsAdjusted_intermediate <- tck_all_merged %>%
  mutate(identified = !is.na(acceptedTaxonID)) %>%
  # filter(!is.na(acceptedTaxonID)) %>%
  group_by(domainID, siteID, plotID, sampleID
           , nlcdClass, decimalLatitude, decimalLongitude, elevation
           , collectDate, year, month, day, dayOfYear, totalSampledArea 
           , lifeStage, identified, acceptedTaxonID) %>%
  summarize(rawCount=sum(Count)) %>%
  ungroup() %>%
  left_join(all_ticks_taxonomy_collapsed) %>%
  left_join(tck_countsIDed_only) %>%
  rename(TotalTickCountsInSample=Count)%>%
  mutate(estimated_species_proportion = rawCount/IDsMadeInSample) %>%
  mutate(estimated_species_proportion = ifelse(is.na(acceptedTaxonID) & rawCount==TotalTickCountsInSample, 1, estimated_species_proportion)) %>%  
  mutate(estimatedCount = estimated_species_proportion*TotalTickCountsInSample) %>%
  select(-identified)

all_ticks_taxonomy_separated <- all_ticks_taxonomy_separated_countsAdjusted_intermediate %>%
  select(-c(TotalTickCountsInSample, IDsMadeInSample, estimated_species_proportion))

# We can vizualize how well raw and adjusted counts map to each other:
all_ticks_taxonomy_separated %>%
  mutate(acceptedTaxonID=factor(ifelse(is.na(acceptedTaxonID),"Unknown",acceptedTaxonID))) %>%
  ggplot() +geom_point(aes(x=log(rawCount), y=log(estimatedCount), col=acceptedTaxonID, pch=lifeStage))+
  scale_color_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                             "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                             "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                             "#8A7C64", "#599861"))


#### Dataset 3: Create dataset of only IXOSCA only ####

# We want to keep zeros in (instead of omitting them) because that will tell us where we DIDN'T sample versus where they were just NOT found.
# Will use "complete" to make sure all real samples are included, regardless of whether there are ticks of that species.

# How do you differentiate from non-taxa ID'ed and other?
# We can see if that particular sample is in tck_tax; if it is in it, it means SOMETHING was ID'ed.
# BUT if the count is zero in tck_field, then it should also be properly zero since we KNOW there were no ticks of any kind there.
all_ticks_IXOSCA <- all_ticks_taxonomy_separated %>%
  filter(acceptedTaxonID=="IXOSCA") %>%
  right_join(all_ticks_taxonomy_collapsed) %>%
  mutate(acceptedTaxonID=ifelse(is.na(acceptedTaxonID), "IXOSCA", acceptedTaxonID) 
         , rawCount=ifelse(is.na(rawCount) & sampleID %in% tck_tax_2$sampleID, 0, rawCount)
         , estimatedCount=ifelse(is.na(estimatedCount) & sampleID %in% tck_tax_2$sampleID, 0, estimatedCount)) %>%
  mutate(rawCount=ifelse(Count==0, 0, rawCount)
         , estimatedCount=ifelse(Count==0, 0, estimatedCount)) %>%
  rename(allSpeciesTotal=Count) 


#### Dataset 4: Create dataset of only AMBAME only ####

# We want to keep zeros in (instead of omitting them) because that will tell us where we DIDN'T sample versus where they were just NOT found.
# Will use "complete" to make sure all real samples are included, regardless of whether there are ticks of that species.

# How do you differentiate from non-taxa ID'ed and other?
# We can see if that particular sample is in tck_tax; if it is in it, it means SOMETHING was ID'ed.
# BUT if the count is zero in tck_field, then it should also be properly zero since we KNOW there were no ticks of any kind there.
all_ticks_AMBAME <- all_ticks_taxonomy_separated %>%
  filter(acceptedTaxonID=="AMBAME") %>%
  right_join(all_ticks_taxonomy_collapsed) %>%
  mutate(acceptedTaxonID=ifelse(is.na(acceptedTaxonID), "AMBAME", acceptedTaxonID)
         , rawCount=ifelse(is.na(rawCount) & sampleID %in% tck_tax_2$sampleID, 0, rawCount)
         , estimatedCount=ifelse(is.na(estimatedCount) &sampleID %in% tck_tax_2$sampleID, 0, estimatedCount)) %>%
  mutate(rawCount=ifelse(Count==0, 0, rawCount)
         , estimatedCount=ifelse(Count==0, 0, estimatedCount)) %>%
  rename(allSpeciesTotal=Count) 

#### WORKING BELOW #####
for ( v in c("taxonomy_collapsed", "taxonomy_separated", "IXOSCA", "AMBAME")) {
  current_data <- get(paste0("all_ticks_",v))

  #### Validation set: months May June July in 2018.
  time_validation_set <- current_data %>%
    filter(month %in% c(5,6,7), year==2018)
  time_test_set <- current_data %>%
    filter(!(sampleID %in% time_validation_set$sampleID))

  #### Validation set: the entire domain, D07
  domain_validation_set <- current_data %>%
    filter(domainID=="D07")
  domain_test_set <- current_data %>%
    filter(domainID!="D07")

  #### Validation set: Sites SCBI, DSNY, STEI, DELA, DEJU
  site_validation_set <- current_data %>%
    filter(siteID%in% c("SCBI","DSNY","STEI","DELA","DEJU"))
  site_test_set <- current_data %>%
    filter(!(sampleID %in% site_validation_set$sampleID))

  #### Validation set: Completely random
  nSample <- length(unique(current_data$sampleID))
  set.seed(23022022)
  validationset <- sample(1:nSample, size=round(nSample*0.30), replace=FALSE)
  random_validation_set <- current_data[validationset,]
  random_test_set <- current_data[-validationset,]

  write_csv(current_data, path=paste0("./data_derived/subset_data/ticks_",v,"/complete_data_",v,".csv"))
  write_csv(time_validation_set, path=paste0("./data_derived/subset_data/ticks_",v,"/time_validationset_",v,".csv"))
  write_csv(time_test_set, path=paste0("./data_derived/subset_data/ticks_",v,"/time_testset_",v,".csv"))
  write_csv(domain_validation_set, path=paste0("./data_derived/subset_data/ticks_",v,"/domain_validationset_",v,".csv"))
  write_csv(domain_test_set, path=paste0("./data_derived/subset_data/ticks_",v,"/domain_testset_",v,".csv"))
  write_csv(random_validation_set, path=paste0("./data_derived/subset_data/ticks_",v,"/random_validationset_",v,".csv"))
  write_csv(random_test_set, path=paste0("./data_derived/subset_data/ticks_",v,"/random_testset_",v,".csv"))
  
}
