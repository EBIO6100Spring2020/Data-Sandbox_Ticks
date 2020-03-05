## Machine learning models on NEON tick abundance data
## BKH; 3/4/2020

library(gbm)
library(dismo) # for gbm.step
library(randomForest) # for randomForest
library(ranger) # for ranger; can give uncertainty in predictions
library(dplyr)

rm(list=ls())

## load data
# complete dataset
data.all <- read.csv("data_derived/subset_data/ticks_taxonomy_collapsed/complete_data_taxonomy_collapsed.csv")
# domain-based training & test datasets
data.domain.train <- read.csv("data_derived/subset_data/ticks_taxonomy_collapsed/domain_testset_taxonomy_collapsed.csv")
data.domain.test <- read.csv("data_derived/subset_data/ticks_taxonomy_collapsed/domain_validationset_taxonomy_collapsed.csv")
# time-based training & test datasets
data.time.train <- read.csv("data_derived/subset_data/ticks_taxonomy_collapsed/time_testset_taxonomy_collapsed.csv")
data.time.test <- read.csv("data_derived/subset_data/ticks_taxonomy_collapsed/time_validationset_taxonomy_collapsed.csv")
# random training & test datasets
data.rand.train <- read.csv("data_derived/subset_data/ticks_taxonomy_collapsed/random_testset_taxonomy_collapsed.csv")
data.rand.test <- read.csv("data_derived/subset_data/ticks_taxonomy_collapsed/random_validationset_taxonomy_collapsed.csv")

## get list of sites that have *some* tick observations
# table with all site counts
data.all %>% group_by(siteID) %>% summarise(siteTotal = sum(Count)) %>%
  ungroup() -> site.table
# table with sites that have ticks
site.table %>% subset(siteTotal>0) %>% pull(siteID) -> tickSites

#View(data.time.test)
#View(filter(data.time.test, siteID %in% tickSites))

## code for BRT and RF taken from supplementary appendix 2 of 
##   Oppel et al 2012, Biological Conservation
##   paper doi: https://doi.org/10.1016/j.biocon.2011.11.013
##   appendix url: https://ars.els-cdn.com/content/image/1-s2.0-S0006320711004319-mmc2.txt

##############################
##### Nymph stage models #####
##############################

predictorIndex <- c(1,2,3,5,10,11) # define columns in data that are predictors
responseIndex <- c(16)


##### DOMAIN-based split #####

# I don't like splitting train/test by domain; removed this code
# That's where "BRT1" & "RF1" went to



##### TIME-based split #####

## Boosted regression trees ##
# Run model on training data
?gbm.step
BRT2 <- gbm.step(data=subset(data.time.train, lifeStage=="Nymph" & siteID %in% tickSites),
                 gbm.x = predictorIndex, # indexes or names predictor variables
                 gbm.y = responseIndex, # indexes or names response variable
                 family = "poisson", # no negative binomial option
                 tree.complexity = 8, # this is what paper used; explore
                 learning.rate = 0.001, # this is what paper used
                 bag.fraction = 0.64 # this is what paper used
)
summary(BRT2) # plotID, month, year, & domainID important (in decreasing order)
# Predict to test data
data.time.test.n <- subset(data.time.test,
                           lifeStage=="Nymph" & siteID %in% tickSites)
data.time.test.n$BRT_pred <- predict.gbm(BRT2, 
                                     data.time.test.n,
                                     n.trees = BRT2$gbm.call$best.trees,
                                     type = "response")
## plot predicted versus observed
plot(log(data.time.test.n$Count+1), log(data.time.test.n$BRT_pred+1))


## Random Forests ##
?randomForest
## run model to training set
# ***randomForest cannot handle categorical predictors with >53 categories,
# so removing plotID from model
RF2 <- randomForest(Count ~ domainID + siteID + nlcdClass + year + month,
                    data=data.time.train,
                    importance=T, # tell it you want imp. of predictors
                    replace=T, # sampling with replacement
                    mtry=5, # #of var rand sampled as candidates at each split
                    ntree=1500, # #of trees to grow 
                    na.action=na.omit) # other settings too, this is what paper did
RF2$importance
# interp. help: https://stats.stackexchange.com/questions/162465/in-a-random-forest-is-larger-incmse-better-or-worse
# %IncMSE: higher number --> more important
# IncNodePurity: higher number --> better variable
# Decreasing Importance: siteID > month > year > nlcdClass > domainID
## predict to test set
data.time.test.n$RF_pred <- predict(RF2, data.time.test.n)
## plot observed versus estimated counts
plot(log(data.time.test.n$Count+1), log(data.time.test.n$RF_pred+1))
# over-estimating when obs count == 0, but pretty decent



##### RANDOM-based split #####

## Boosted regression trees ##

## Run model on training data
BRT3 <- gbm.step(data=subset(data.rand.train, lifeStage=="Nymph" & siteID %in% tickSites),
                 gbm.x = predictorIndex, # indexes or names predictor variables
                 gbm.y = responseIndex, # indexes or names response variable
                 family = "poisson", # no negative binomial option
                 tree.complexity = 8, # this is what paper used; explore
                 learning.rate = 0.001, # this is what paper used
                 bag.fraction = 0.64 # this is what paper used
)
summary(BRT3) # plotID, year, & month most important
## Predict to test data
data.rand.test.n <- subset(data.rand.test,
                           lifeStage=="Nymph" & siteID %in% tickSites)
?predict.gbm
data.rand.test.n$BRT_pred <- predict.gbm(BRT3, 
                                   data.rand.test.n,
                                   n.trees = BRT3$gbm.call$best.trees,
                                   type = "response")
## plot predicted versus observed
plot(log(data.rand.test.n$Count+1), log(data.rand.test.n$BRT_pred+1))


## Random forests ##
?randomForest
## run model to training set
# ***randomForest cannot handle categorical predictors with >53 categories,
# so removing plotID from model
RF3 <- randomForest(Count ~ domainID + siteID + nlcdClass + year + month,
                    data=data.rand.train,
                    importance=T, # tell it you want imp. of predictors
                    replace=T, # sampling with replacement
                    mtry=5, # #of var rand sampled as candidates at each split
                    ntree=1500, # #of trees to grow 
                    na.action=na.omit) # other settings too, this is what paper did
RF3$importance
# Decreasing Importance: month > siteID > nlcdClass > year > domain
## predict to test set
data.rand.test.n$RF_pred <- predict(RF3, data.rand.test.n)
# ERROR: Type of predictors in new data do not match that of the training data.
# not sure why it threw here but not above...
# need to look into it ***

## plot observed counts versus predicted
plot(log(data.rand.test.n$Count+1), log(data.rand.test.n$RF_pred+1))



## BKH only: save workspace (models take f**king forever to run on my laptop)
# *** note to self: add below file to gitIgnore
#save.image(file="machineLearnSpace.RData")
#load.image(file="machineLearnSpace.RData")