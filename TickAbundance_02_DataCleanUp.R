# export tick abundance training and testing datasets for gams
library(dplyr)
library(tidyr)
library(forcats)
# read in training dataset
tick_abun <- read.csv("data_derived/subset_data/ticks_IXOSCA/time_testset_IXOSCA.csv")

tick_abun %>% filter(lifeStage == "Nymph") %>% droplevels() %>% # get only nymphs
 mutate( # add other columns and scale
  collectDate = as.Date(collectDate),
  nymph_presence = case_when(estimatedCount > 0 ~ 1, estimatedCount == 0 ~ 0),
  season = case_when(month %in% c(5,6,7,8,9)~1,!month %in% c(5,6,7,8,9)~0 ),
  density = estimatedCount/totalSampledArea,
  logSampledArea = log(totalSampledArea),
  sDecimalLongitude = scale(decimalLongitude),
  sDecimalLatitude = scale(decimalLatitude),
  sDayofYear = scale(dayOfYear),
  estimatedCount = as.integer(estimatedCount),
  fMonth = as.factor(month),
  fYear = as.factor(year),
  yearPlot = as.factor(paste(plotID, year, sep = "_"))
) %>% filter(!is.na(density)) -> nymph_all

# get total ticks counted at each plot/site/domain over all visits, add to data
nymph_all %>% group_by(plotID) %>% mutate(totalTicks_plot = sum(estimatedCount, na.rm = TRUE)) %>% ungroup() -> nymph_all
nymph_all %>% group_by(domainID) %>% mutate(totalTicks_domain = sum(estimatedCount, na.rm = TRUE)) %>% ungroup() -> nymph_all
nymph_all %>% group_by(siteID) %>% mutate(totalTicks_site = sum(estimatedCount, na.rm = TRUE)) %>% ungroup() -> nymph_all
nymph_all %>% group_by(nlcdClass) %>% mutate(totalTicks_nlcd = sum(estimatedCount, na.rm = TRUE)) %>% ungroup() -> nymph_all

# get rid of domains and nlcd classes that have no ticks (but keep sites/plots that never have ticks)
nymph_all %>% filter(totalTicks_domain>0) %>% droplevels() -> nymph_all
nymph_all %>% filter(totalTicks_nlcd>0) %>% droplevels() -> nymph_all


# summarise total ticks per plot
nymph_all %>% group_by(plotID) %>% summarise(totalTicks_plot = sum(estimatedCount, na.rm = TRUE))  %>% ungroup() -> total_ticks

# reorder factors
nymph_all$nlcdClass <- fct_reorder(nymph_all$nlcdClass, nymph_all$totalTicks_nlcd)
nymph_all$domainID <- fct_reorder(nymph_all$domainID, nymph_all$totalTicks_domain)
nymph_all$plotID <- fct_reorder(nymph_all$plotID, nymph_all$totalTicks_plot)

# export data
saveRDS(nymph_all, "data_derived/nymph_abun_all.Rdata")



#### REPEAT WITH TEST DATASET
# read in/reformat the test dataset
tick_abun_test <- read.csv("data_derived/subset_data/ticks_IXOSCA/time_validationset_IXOSCA.csv")

tick_abun_test %>% filter(lifeStage == "Nymph") %>% droplevels() %>% # get only nymphs
  mutate( # add other columns and scale
    collectDate = as.Date(collectDate),
    nymph_presence = case_when(estimatedCount > 0 ~ 1, estimatedCount == 0 ~ 0),
    season = case_when(month %in% c(5,6,7,8,9)~1,!month %in% c(5,6,7,8,9)~0 ),
    density = estimatedCount/totalSampledArea,
    logSampledArea = log(totalSampledArea),
    sDecimalLongitude = scale(decimalLongitude, center = -91.2, scale = 15.1 ), # scale by same amt as training dataset!
    sDecimalLatitude = scale(decimalLatitude, center = 38.1, scale = 7.26), # scale by same amt as training
    sDayofYear = scale(dayOfYear, center = 199, scale = 68.1),
    estimatedCount = as.integer(estimatedCount),
    fMonth = as.factor(month),
    fYear = as.factor(year),
    yearPlot = as.factor(paste(plotID, year, sep = "_"))
  ) %>% filter(!is.na(density)) -> nymph_test


# get total ticks counted at each plot/site/domain over all visits, add to data
nymph_test %>% group_by(plotID) %>% mutate(totalTicks_plot = sum(estimatedCount, na.rm = TRUE)) %>% ungroup() -> nymph_test
nymph_test %>% group_by(domainID) %>% mutate(totalTicks_domain = sum(estimatedCount, na.rm = TRUE)) %>% ungroup() -> nymph_test
nymph_test %>% group_by(siteID) %>% mutate(totalTicks_site = sum(estimatedCount, na.rm = TRUE)) %>% ungroup() -> nymph_test
nymph_test %>% group_by(nlcdClass) %>% mutate(totalTicks_nlcd = sum(estimatedCount, na.rm = TRUE)) %>% ungroup() -> nymph_test

# get rid of domains and nlcd classes that have no ticks (but keep sites/plots that never have ticks)
nymph_test %>% filter(totalTicks_domain>0) %>% droplevels() -> nymph_test
nymph_test %>% filter(totalTicks_nlcd>0) %>% droplevels() -> nymph_test




# reorder factors
nymph_test$nlcdClass <- fct_reorder(nymph_test$nlcdClass, nymph_test$totalTicks_nlcd)
nymph_test$domainID <- fct_reorder(nymph_test$domainID, nymph_test$totalTicks_domain)
nymph_test$plotID <- fct_reorder(nymph_test$plotID, nymph_test$totalTicks_plot)

# export data
saveRDS(nymph_test, "data_derived/nymph_abun_test.Rdata")

