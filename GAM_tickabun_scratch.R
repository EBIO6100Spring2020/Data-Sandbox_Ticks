###### GAMs with Tick Data
## Author: Wynne Moss
## Date: March 17 2020
## About: 
# trying to use generalized additive models to predict tick abundance across time/space
# analysis uses NEON's tick datasets, cleaned previously by M. Chen
library(dplyr)
library(ggplot2)
library(tidyr)
# library(padr)
library(rjags)
library(mgcv)

#### Read in and clean up data ####
# read in training set
tick_abun <- read.csv("data_derived/subset_data/ticks_IXOSCA/time_testset_IXOSCA.csv")
# read in test set
tick_abun_test <- read.csv("data_derived/subset_data/ticks_IXOSCA/time_validationset_IXOSCA.csv")

# get only nymphs from plots that have non zero counts
tick_abun %>% group_by(plotID) %>% filter(lifeStage == "Nymph") %>% mutate(totalTicks = sum(estimatedCount, na.rm=TRUE)) %>% filter(totalTicks > 0) %>% ungroup() %>% droplevels() -> nymph

# convert counts to density using the drag length
nymph <-nymph %>% mutate(density = estimatedCount/totalSampledArea) %>% filter(!is.na(density))
# fix date
nymph$collectDate <- as.Date(nymph$collectDate)

#### Run GAMs #####
# gaussian
g.g <- gam(log(density+1) ~ nlcdClass + s(month) + s(decimalLongitude) , data = nymph, method = "REML")
plot(g.g, residuals = TRUE, pch = 1, all.terms = TRUE)
summary(g.g)
gam.check(g.g)
# log-density still not gaussian, consider another family

# try a negative binomial
g.nb <- gam(estimatedCount ~ nlcdClass + s(month) , offset = log(totalSampledArea), data = nymph, family = nb, method = "REML")
plot(g.nb, residuals = TRUE, pch = 1, all.terms = TRUE)
# this looks super weird...


# zero inflated Poisson
g.zip <- gam(estimatedCount ~ nlcdClass + s(month) , offset = log(totalSampledArea), data = nymph, family = ziP, method = "REML")
plot(g.zip, residuals = TRUE, pch = 1, all.terms = TRUE)

# a two stage Poisson (one linear predictor controls prob of presence, and the other controls the mean)

g.ziplss <- gam(list(estimatedCount ~ nlcdClass + s(month)), data = nymph, family = ziplss(), method = "REML")
