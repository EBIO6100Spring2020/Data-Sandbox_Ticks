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
library(brms)
library(schoenberg)

#### Read in and clean up data ####
# read in training set
tick_abun <- read.csv("data_derived/subset_data/ticks_IXOSCA/time_testset_IXOSCA.csv")
# read in test (validation) set
tick_abun_test <- read.csv("data_derived/subset_data/ticks_IXOSCA/time_validationset_IXOSCA.csv")

# filter only plots that have non zero counts of ticks, and only nymphs
tick_abun %>% group_by(plotID) %>% filter(lifeStage == "Nymph") %>% mutate(totalTicks = sum(estimatedCount, na.rm=TRUE)) %>% filter(totalTicks > 0) %>% ungroup() %>% droplevels() -> nymph

# convert counts to density using the drag length
nymph <-nymph %>% mutate(density = estimatedCount/totalSampledArea) %>% filter(!is.na(density))

# fix date
nymph$collectDate <- as.Date(nymph$collectDate)

# add p/a column
nymph %>% mutate(nymph_presence = case_when(estimatedCount > 0 ~ 1, estimatedCount == 0 ~ 0)) -> nymph

# add scaled/logged variables
nymph$logSampledArea <- log(nymph$totalSampledArea)
nymph$sDecimalLongitude <- scale(nymph$decimalLongitude)
nymph$sDayofYear <- scale(nymph$dayOfYear)

#### Run simple GAMs #####
# gaussian
gam_dens <- gam(density ~ nlcdClass + siteID + s(month) + s(decimalLongitude),
             data = nymph, method = "REML")
summary(gam_dens)
plot(gam_dens, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
# this doesn't look great, likely due to distribution of density:
hist(nymph$density)
# does log density perform any better?
hist(log(nymph$density+1))
# not really, and still very zero-inflated
gam_log_dens <- gam(log(density+1) ~ nlcdClass + + siteID+ s(month) + s(decimalLongitude) , data = nymph, method = "REML")
plot(gam_log_dens, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
summary(gam_log_dens)
gam.check(gam_log_dens)
# log-density still not gaussian, consider another family

#### Try GAM for presence/absence ####
table(nymph$nymph_presence, nymph$domainID)
table(nymph$nymph_presence, nymph$plotID)
table(nymph$nymph_presence, nymph$nlcdClass)

gam_binom <- gam(nymph_presence ~ s())
#### Run GAM with negative binomial family ####
# try a negative binomial
gam_count_nb <- gam(estimatedCount ~ s(month) + s(plotID, bs = "re"), offset = log(totalSampledArea), data = nymph, family = nb, method = "REML")

summary(gam_count_nb)

plot(gam_count_nb, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
# this looks super weird...

gam_count_nb_2 <- gam(list(
  estimatedCount ~ s(month) + nlcdClass + offset(log(totalSampledArea))+s(plotID, bs = "re"),
  ~1),
  data = nymph, family = ziplss())

# zero inflated Poisson
gam_count_zip <- gam(estimatedCount ~ nlcdClass + s(month) , offset = log(totalSampledArea), data = nymph, family = ziP, method = "REML")
plot(gam_count_zip, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)

gam_count_zip_2 <- gam(estimatedCount ~ nlcdClass + siteID + s(month) + s(decimalLongitude) , offset = log(totalSampledArea), data = nymph, family = ziP, method = "REML")
plot(gam_count_zip_2, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
# what is going on with month?


gam_count_zip_3 <- gam(estimatedCount ~ nlcdClass + siteID + s(sDayofYear) + s(decimalLongitude) , offset = logSampledArea, data = nymph, family = ziP, method = "REML")
plot(gam_count_zip_3, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)


# a two stage Poisson (one linear predictor controls prob of presence, and the other controls the mean)

gam_count_ziplss <- gam(list(estimatedCount ~ s(sDayofYear),
                     ~ nlcdClass + siteID + s(month) + s(sDecimalLongitude)),
            offset = logSampledArea, data = nymph, family = ziplss())
plot(gam_count_ziplss, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
summary(gam_count_ziplss)


# try an interaction
gam_count_ziplss_2 <- gam(list(estimatedCount ~ te(sDayofYear, sDecimalLongitude), # controls the number
                             ~ nlcdClass + siteID + s(month)), # controls presence/absence
                        offset = logSampledArea, data = nymph, family = ziplss())
plot(gam_count_ziplss_2, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
par(mfrow=c(1,1))
vis.gam(gam_count_ziplss_2, view = c("sDayofYear", "sDecimalLongitude"), plot.type = "contour", scheme = TRUE, pch = 1, residuals = TRUE, cex = 1, cex.axis = 1)

summary(gam_count_ziplss_2)
