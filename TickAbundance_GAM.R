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
library(mgcv) # fitting GAMS
library(brms) # fitting Bayesian models
library(schoenberg) # for brms plots (not used yet)
library(itsadug) # for plotting RE with GAM

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
# add season (dummy code: 0 = not tick season, 1 = tick season)
nymph %>% mutate(season = case_when(month %in% c(5,6,7,8,9)~1,!month %in% c(5,6,7,8,9)~0 )) -> nymph
# add scaled/logged variables
nymph$logSampledArea <- log(nymph$totalSampledArea)
nymph$sDecimalLongitude <- scale(nymph$decimalLongitude)
nymph$sDecimalLatitude <- scale(nymph$decimalLatitude)
nymph$sDayofYear <- scale(nymph$dayOfYear)
nymph$estimatedCount <- as.integer(nymph$estimatedCount)
nymph$fMonth <- as.factor(nymph$month)
nymph$fYear <- as.factor(nymph$month)

#### Make some exploratory plots ####
nymph %>% ggplot(aes(x=dayofYear, y = density)
#### Run simple GAMs #####
# gaussian
gam_dens <- gam(density ~ nlcdClass + siteID + s(sDayofYear) + s(sDecimalLongitude),
             data = nymph, method = "REML")
summary(gam_dens)
plot(gam_dens, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
# this doesn't look great, likely due to distribution of density:
hist(nymph$density)
# does log density perform any better?
hist(log(nymph$density+1))
# not really, and still very zero-inflated
gam_log_dens <- gam(log(density+1) ~ nlcdClass + siteID+ s(sDayofYear) + s(sDecimalLongitude) , data = nymph, method = "REML")
plot(gam_log_dens, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
summary(gam_log_dens)
gam.check(gam_log_dens)
# log-density still not gaussian, consider another family

#### Binomial GAM with presence/absence ####

gam_binom <- gam(nymph_presence ~ nlcdClass + siteID + s(sDayofYear) + s(sDecimalLongitude),
                 data = nymph, family = binomial)
summary(gam_binom)
plot(gam_binom, pch = 1, all.terms = TRUE, pages = 1, trans = plogis)
termplot(gam_binom, pch = 1, rug = FALSE,
         se = TRUE, terms = "nlcdClass", 
         las =2)
termplot(gam_binom, pch = 1, rug = FALSE,
         se = TRUE, terms = "siteID", 
         las =2)
# kept all samples that observation said there were ticks even if model said there weren't supposed to be ticks
# observed zeros -- probability of that particular sample being a 0, if the pr>0.5 it's a "poisson 0". 
# things that weren't a hurdle 0 those data are used to fit the next model


#### Add Random Effects to Binomial GAM  ####
gam_binom_re <- gam(nymph_presence ~ nlcdClass+ s(sDayofYear) + te(sDecimalLongitude, sDecimalLatitude) +
                       s(siteID, bs = "re")+s(plotID, bs = "re"), data = nymph, family = binomial)
summary(gam_binom_re)
# site ID doesn't matter if plotID is in the model
plot(gam_binom_re, pages = 1)
plot_smooth(gam_binom_re, view = "sDayofYear", cond=list(nlcdClass = "deciduousForest"),
            rm.ranef=TRUE, rug = FALSE, transform = plogis)

plot_smooth(gam_binom_re, view = "sDayofYear", cond=list(nlcdClass = "deciduousForest", plotID = "TALL_008"),
            col = "blue", rug = FALSE, transform = plogis)
plot_smooth(gam_binom_re, view = "sDayofYear", cond=list(nlcdClass = "deciduousForest", plotID = "TREE_017"),
            col = "red", rug = FALSE,  transform = plogis, add = TRUE)
# looks like plotID >> siteID

# random slopes 
gam_binom_rs <- gam(nymph_presence ~ nlcdClass+ s(sDayofYear) + te(sDecimalLongitude, sDecimalLatitude) +
                      +s(plotID, bs = "re") + s(plotID, sDayofYear, bs = "re"), data = nymph, family = binomial)
plot_smooth(gam_binom_rs, view = "sDayofYear", cond=list(nlcdClass = "deciduousForest"),
            rm.ranef=TRUE, rug = FALSE, transform = plogis)

plot_smooth(gam_binom_rs, view = "sDayofYear", cond=list(nlcdClass = "deciduousForest", plotID = "TALL_008"),
            col = "blue", rug = FALSE, transform = plogis)
plot_smooth(gam_binom_rs, view = "sDayofYear", cond=list(nlcdClass = "deciduousForest", plotID = "TREE_017"),
            col = "red", rug = FALSE,  transform = plogis, add = TRUE)
plot_smooth(gam_binom_rs, view = "sDayofYear", cond=list(nlcdClass = "deciduousForest", plotID = "OSBS_001"),
            col = "green", rug = FALSE,  transform = plogis, add = TRUE)
plot_smooth(gam_binom_rs, view = "sDayofYear", cond=list(nlcdClass = "deciduousForest", plotID = "UNDE_019"),
            col = "purple", rug = FALSE,  transform = plogis, add = TRUE)

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

gam_count_ziplss <- gam(list(estimatedCount ~ s(sDayofYear) + s(plotID, bs = "re"), # controls mean
                     ~ nlcdClass +  s(month) + te(sDecimalLatitude, sDecimalLongitude) + s(plotID, bs = "re")), # controls probability of presence
            offset = logSampledArea, data = nymph, family = ziplss())
start_event(data = nymph, column = "collectDate", event = c("plotID"))
plot(gam_count_ziplss, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
summary(gam_count_ziplss)

# how good is this model?
acf_resid(gam_count_ziplss)

# treat explicitly as a time series?
nymph$Time <- as.integer(nymph$collectDate - min(nymph$collectDate))
nymph <- data.frame(nymph)
nymph <- start_event(data = nymph,  event = c("plotID", "siteID"))
valRho  <- acf(resid(gam_count_ziplss), plot=FALSE)$acf[2]

gam_count_ziplss_ts <- gam(list(estimatedCount ~ s(sDayofYear) + s(plotID, bs = "re") + s(Time, by = plotID), # controls mean
                             ~ nlcdClass +  s(month) + te(sDecimalLatitude, sDecimalLongitude) + s(plotID, bs = "re")), # controls probability of presence
                        AR.start = nymph$start.event, rho = valRho,
                        offset = logSampledArea, data = nymph, family = ziplss())
check_resid(gam_count_ziplss_ts)
# plotIDdaysSince <- paste(nymph$daysSince, nymph$plotID)
# length(unique(plotIDdaysSince))
# day of year vs. month


# try an interaction
gam_count_ziplss_2 <- gam(list(estimatedCount ~ te(sDayofYear, sDecimalLongitude), # controls the number
                             ~ nlcdClass + siteID + s(month)), # controls presence/absence
                        offset = logSampledArea, data = nymph, family = ziplss())
plot(gam_count_ziplss_2, residuals = TRUE, pch = 1, all.terms = TRUE, pages = 1)
par(mfrow=c(1,1))
vis.gam(gam_count_ziplss_2, view = c("sDayofYear", "sDecimalLongitude"), plot.type = "contour", scheme = TRUE, pch = 1, residuals = TRUE, cex = 1, cex.axis = 1)

summary(gam_count_ziplss_2)

### Using brms ####
# notes: can't use tensor smooths in brm
# need to use tensor product smooths
# smooths are basically random effects
# "wiggly" parts of spline are random effects
# variance parameter controls degree of wigliness
# smooth parts are the fixed effects

# try a very simple model
gam_count_zinb1 <- brm(estimatedCount ~ nlcdClass + as.factor(month) + offset(logSampledArea), data = nymph, family = zero_inflated_poisson())
summary(gam_count_zinb1)
# zero inflation is high! 0.60
# factor month is probably not legit

# add in the random effect for plot
gam_count_zinb2<- brm(estimatedCount ~ nlcdClass + fMonth +(1|plotID)+ offset(logSampledArea) , data = nymph, family = zero_inflated_poisson(), iter =1000, chains = 2)
summary(gam_count_zinb2)
# add in some predictors for zero inflation. zi ~ things that predict 0s
gam_count_zinb3 <- brm(bf(
  estimatedCount ~ nlcdClass + fMonth+ (1|plotID) + offset(logSampledArea),
  zi ~ nlcdClass + (season|plotID)),
  data = nymph, family = zero_inflated_poisson(), iter =1000, chains = 2)
summary(gam_count_zinb3)
plot(conditional_effects(gam_count_zinb3), ask = FALSE)
# add in a smooth function for month and random effect for year nested within plot
gam_count_zinb4 <- brm(bf(
  estimatedCount ~ nlcdClass + s(sDayofYear)+ (1|plotID/fYear) + offset(logSampledArea),
  zi ~ nlcdClass + (season|plotID)),
  data = nymph, family = zero_inflated_poisson(), iter =1000, chains = 2)
summary(gam_count_zinb4)

newdata <- expand.grid(
  "sDayofYear" = seq(-2,2, by = .1),
  "nlcdClass" = levels(nymph$nlcdClass),
  "plotID" = "new",
  "fYear" = "new",
  "logSampledArea" = 0,
  "season" = 1
)
ranef(gam_count_zinb4)
preds.zinb4 <- predict(gam_count_zinb4, newdata = newdata, re_formula = NA
                       )
head(preds.zinb4)
newdata <- cbind(newdata, preds.zinb4)
head(newdata)
ggplot(newdata, aes(x=sDayofYear, y = Estimate, group = nlcdClass))+
  geom_line(aes(color = nlcdClass))

  Estimate~sDayofYear, newdata)
#### To do list ####
# refine models using diagnostics
# model selection
# test on other dataset and quantify error
# extract uncertainties on predictions 
# the extreme events: what do they have in common? 