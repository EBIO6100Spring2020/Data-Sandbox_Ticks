######## GAM model with manual hurdle component ##########
library(tidyverse)
library(pscl) # hurdle models
library(lme4) # mixed effect models
library(boot) # logit/inv logit
library(car) # Anova (type III)
library(mgcv) # for GAMs
library(randomForest) # for random forest

tck <- read.csv("data_derived/MASTER_all_tck_data_merged.csv")

#### Filter tick dataset and adjust ####

## How many plots were tested, and how many times?
tck %>%
  mutate(tested=ifelse(is.na(testingID),"no","yes")) %>%
  select(plotID, tested) %>% table(useNA="ifany")

## Filter tck to only include samples that tested Borrelia_sp.; count number of different life stages; create proportiontested/positive tc
tck_allsamples_borr <- tck %>%
  mutate(tested=ifelse(is.na(Borrelia_sp.),0,1) # Was each tick  ever tested for borrellia?
         ,isPositive=ifelse(Borrelia_sp.=="Positive", 1,0) # And if it was tested, was it positive or negative?
  ) %>%
  group_by(domainID, siteID, nlcdClass, plotID, elevation, collectDate, dayOfYear, year, month) %>% # collapse by sample-- summing all tck counts together
  summarize(numberTested=sum(tested, na.rm=TRUE) # number of tested ticks in that sample
            ,n=n() # total ticks in that sample
            , numberPositive=sum(isPositive,na.rm = TRUE) # number of positive ticks in that sample
            , nAdult=sum(lifeStage=="Adult")
            , nNymph=sum(lifeStage=="Nymph")
            , nLarva=sum(lifeStage=="Larva")) %>%
  ungroup() %>%
  mutate(proportionTested=numberTested/nNymph # proportion of all nymph ticks tested-- only nymphs were ever tested.
         , proportionPositive=numberPositive/numberTested
         , tested = ifelse(numberTested > 0 , TRUE, FALSE) # new true/false tested, which is summed across ticks
         , testingStatus = ifelse(numberTested > 0, "Tested", ifelse(nNymph>0, "Nymphs present, not tested", "No nymphs"))
         ) %>%
  mutate(nlcdClass=factor(nlcdClass, levels=c("emergentHerbaceousWetlands","cultivatedCrops","pastureHay","grasslandHerbaceous"
                                              ,"dwarfScrub","shrubScrub","sedgeHerbaceous"
                                              ,"woodyWetlands","deciduousForest","evergreenForest","mixedForest"))) %>%
  mutate(borrPresent = ifelse(numberPositive>0,1,0)
         , domainID = factor(as.character(domainID))
         , siteID = factor(as.character(siteID))
         , plotID = factor(as.character(plotID))
         , year = factor(year)
         , tckDensity = sum(c(nLarva, nNymph, nAdult))/n
         , ntckDensity = nNymph/n
         , NLtckDensity = sum(c(nNymph, nAdult), na.rm = TRUE)/n
         , atckDensity = nAdult/n)  %>%
  mutate(lognymphDensity=log(ntckDensity)
         , logNLtckDensity=log(NLtckDensity)
         , logtckDensity = log(tckDensity)
         , logadultDensity=log(atckDensity))

## Plot over time, by year, by plot, to see how zero-inflated it is.
tck_allsamples_borr %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=plotID, fill=proportionPositive, col=testingStatus), pch=21) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_manual(values=c(Tested="blue", `Nymphs present, not tested`="grey", `No nymphs`="white")) +
  facet_grid(.~year)

## Let's filter by plots that NEVER had any positive borrelia
infectedPlots <- tck %>%
  filter(Borrelia_sp.=="Positive") %>%
  select(plotID) %>%pull() %>% unique()

## Filtering same data as above; only now removing non-infected plots.
tck_borrelia_positivePlots <- tck_allsamples_borr %>%
  filter(plotID %in% as.character(infectedPlots))

tck_borrelia_positivePlots %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=plotID, fill=proportionPositive, col=testingStatus), pch=21) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_manual(values=c(Tested="blue", `Nymphs present, not tested`="grey", `No nymphs`="white")) +
  facet_grid(.~year)

## There are STILL a lot of zeros-- zeros where there are nymphs but they were never tested.
## Let's remove those.
tck_borrelia <- tck_borrelia_positivePlots %>%
  mutate(numberPositive=ifelse(numberTested==0,NA,numberPositive)) %>% # Make sure that numberPositive is not artificially zero-- if there were no tests, it should be NA
  filter(numberTested!=0) # get rid of all samples where the didn't actually test.

tck_borrelia %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=plotID, fill=proportionPositive, col=testingStatus), pch=21) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_manual(values=c(Tested="blue", `Nymphs present, not tested`="grey", `No nymphs`="white")) +
  facet_grid(.~year)

nrow(tck_borrelia)
# There are only 311 samples left-- still zero inflated, but slightly better.

#### Preliminary Plotting ####

# First, what is the distribution of Borrelia prevalence?
tck_borrelia %>%
  ggplot() +geom_histogram(aes(x=log(numberPositive+1)), bins=100)
# histogram without any zeros
tck_borrelia %>%
  filter(numberPositive>0) %>%
  ggplot() +geom_histogram(aes(x=(numberPositive)), bins=100)
# histogram log abundance, no zeros
tck_borrelia %>%
  filter(numberPositive>0) %>%
  ggplot() +geom_histogram(aes(x=log(numberPositive)), bins=50)

# What is relationship between proportion positive and number positive
tck_borrelia %>%
  ggplot() +geom_point(aes(x=log(numberPositive), y=log(proportionPositive)))
# Now, what is the relationship between number tested and nNymph? Might be a binning artifact?
tck_borrelia %>%
  ggplot() +geom_point(aes(x=log(nNymph), y=log(numberTested)))
# Just get rid of all data where they didn't test most nymphs-- see how many that is
tck_borrelia %>%
  filter(abs(nNymph-numberTested)>0) %>% select(plotID, collectDate, nNymph, numberTested)
# Most of these samples are not a problem, except HARV_004 (2015) and maybe SCI_002 (2017) and SERC_001 (2017-07). 
# Let's get rid of only HARV_004; I think all the rest are fine.
tck_borrelia_adj <- tck_borrelia %>%
  filter(abs(nNymph-numberTested)<3) # HARV_004 is the only site that differs nNymph and numberTested by more than 3


#### Fitting GAMs with a hurdle model component ####

# First, model a binomial component with dayOfYear, nlcdClass, and non-larval tick density (NLtckDensity) as predictors of borrellia presence/absence.
# Plot is a random effect, which I include as a varying-intercept, varying-slope model.
mod.gambin_all <- gam(borrPresent ~ s(dayOfYear) + nlcdClass + 
                    s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + 
                    s(year, bs="re") + s(dayOfYear, year, bs="re") +
                    s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                    s(elevation) + 
                    offset(log(numberTested)) 
                # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                , data=tck_borrelia_adj
                , method="REML"
                , family=binomial)
summary(mod.gambin_all)

mod.gambin_noanytckDensity <- gam(borrPresent ~ s(dayOfYear) + nlcdClass + 
                        s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + 
                        s(year, bs="re") + s(dayOfYear, year, bs="re") +
                        # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + 
                        s(elevation) + 
                        offset(log(numberTested)) 
                      # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                      , data=tck_borrelia_adj
                      , method="REML"
                      , family=binomial)
summary(mod.gambin_noanytckDensity)
c(mod.gambin_all$aic-mod.gambin_noanytckDensity$aic) # fit is actually better if we DON'T include tck density

mod.gambin_noelevation <- gam(borrPresent ~ s(dayOfYear) + nlcdClass + 
                                s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + 
                                s(year, bs="re") + s(dayOfYear, year, bs="re") +
                                # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                                # s(elevation) + 
                                offset(log(numberTested)) 
                              # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                              , data=tck_borrelia_adj
                              , method="REML"
                              , family=binomial)
c(mod.gambin_noanytckDensity$aic - mod.gambin_noelevation$aic) # fit is better WITH elevation

# Check how fit is if we do the reverse-- include tck density but exclude elevation
mod.gambin_noelevation_wdensity <- gam(borrPresent ~ s(dayOfYear) + nlcdClass + 
                                s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + 
                                s(year, bs="re") + s(dayOfYear, year, bs="re") +
                                s(logNLtckDensity) +
                                # s(elevation) + 
                                offset(log(numberTested)) 
                              # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                              , data=tck_borrelia_adj
                              , method="REML"
                              , family=binomial)
c(mod.gambin_noanytckDensity$aic, mod.gambin_noelevation$aic, mod.gambin_noelevation_wdensity$aic) # fit is better WITH elevation

mod.gambin_noyear <- gam(borrPresent ~ s(dayOfYear) + nlcdClass + 
                                s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + 
                                # s(year, bs="re") + s(dayOfYear, year, bs="re") +
                                # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                                s(elevation) +
                                offset(log(numberTested)) 
                              # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                              , data=tck_borrelia_adj
                              , method="REML"
                              , family=binomial)
c(mod.gambin_noanytckDensity$aic - mod.gambin_noyear$aic) # fit is better WITH year

mod.gambin_noplot <- gam(borrPresent ~ s(dayOfYear) + nlcdClass +
                           # s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
                           s(year, bs="re") + s(dayOfYear, year, bs="re") +
                           # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                           s(elevation) +
                           offset(log(numberTested))
                         # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                         , data=tck_borrelia_adj
                         , method="REML"
                         , family=binomial)
c(mod.gambin_noanytckDensity$aic - mod.gambin_noplot$aic) # fit is better WITH plot


mod.gambin_nonlcdClass <- gam(borrPresent ~ s(dayOfYear) + #nlcdClass +
                           s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
                           s(year, bs="re") + s(dayOfYear, year, bs="re") +
                           # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                           s(elevation) +
                           offset(log(numberTested))
                         # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                         , data=tck_borrelia_adj
                         , method="REML"
                         , family=binomial)
c(mod.gambin_noanytckDensity$aic - mod.gambin_nonlcdClass$aic) # fit is actually better WITHOUT nlcdClass

mod.gambin_nodayOfYear <- gam(borrPresent ~ #s(dayOfYear) + #nlcdClass +
                                s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
                                s(year, bs="re") + s(dayOfYear, year, bs="re") +
                                # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                                s(elevation) +
                                offset(log(numberTested))
                              # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                              , data=tck_borrelia_adj
                              , method="REML"
                              , family=binomial)
c(mod.gambin_nonlcdClass$aic - mod.gambin_nodayOfYear$aic) # fit is  better WITH dayOfYear

# Double check nlcdClass and elevation are not conflated?
mod.gambin_noelev_wnlcd <- gam(borrPresent ~ nlcdClass +# s(dayOfYear)  +
                                s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
                                s(year, bs="re") + s(dayOfYear, year, bs="re") +
                                # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                                # s(elevation) +
                                offset(log(numberTested))
                              # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                              , data=tck_borrelia_adj
                              , method="REML"
                              , family=binomial)
mod.gambin_noelev_wnlcd$aic
mod.gambin_nonlcdClass$aic
# Yup, nlcd Class doesn't contain much info.

# Now check; do I need interactions for random effects?
summary(mod.gambin_nonlcdClass) # best fitting model so far
mod.gambin_noplot <- gam(borrPresent ~ s(dayOfYear) + #nlcdClass +
                                s(dayOfYear,plotID, bs="re") +# s(plotID, bs="re") + 
                                s(year, bs="re") + s(dayOfYear, year, bs="re") +
                                # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                                s(elevation) +
                                offset(log(numberTested))
                              # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                              , data=tck_borrelia_adj
                              , method="REML"
                              , family=binomial)
c(mod.gambin_nonlcdClass$aic - mod.gambin_noplot$aic) # model is slightly worse without plots, but BARELY. Keep in.


summary(mod.gambin_nonlcdClass) # best fitting model so far
mod.gambin_noplotInter <- gam(borrPresent ~ s(dayOfYear) + #nlcdClass +
                           s(plotID, bs="re") + 
                           s(year, bs="re") + s(dayOfYear, year, bs="re") +
                           # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                           s(elevation) +
                           offset(log(numberTested))
                         # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                         , data=tck_borrelia_adj
                         , method="REML"
                         , family=binomial)
c(mod.gambin_nonlcdClass$aic - mod.gambin_noplotInter$aic) # model is  worse without plot interaction

summary(mod.gambin_nonlcdClass) # best fitting model so far
mod.gambin_noyearInter <- gam(borrPresent ~ s(dayOfYear) + #nlcdClass +
                           s(dayOfYear,plotID, bs="re") + s(plotID, bs="re") + 
                           s(year, bs="re") +# s(dayOfYear, year, bs="re") +
                           # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                           s(elevation) +
                           offset(log(numberTested))
                         # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                         , data=tck_borrelia_adj
                         , method="REML"
                         , family=binomial)
c(mod.gambin_nonlcdClass$aic - mod.gambin_noyearInter$aic) # model is slightly better with interaction of year, dayofyear;Keep in.

# See if domain has an effect
mod.gambin_domain <- gam(borrPresent ~ s(dayOfYear) + #nlcdClass +
                                s(dayOfYear,plotID, bs="re") + s(plotID, bs="re") + 
                                s(year, bs="re") +# s(dayOfYear, year, bs="re") +
                                # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                                s(elevation) + s(domainID, bs="re") + s(domainID, dayOfYear, bs="re") +
                                offset(log(numberTested))
                              # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                              , data=tck_borrelia_adj
                              , method="REML"
                              , family=binomial)
c(mod.gambin_nonlcdClass$aic - mod.gambin_domain$aic)

# What if we took elevation away now?
mod.gambin_domainnoelev <- gam(borrPresent ~ s(dayOfYear) + #nlcdClass +
                           s(dayOfYear,plotID, bs="re") + s(plotID, bs="re") + 
                           s(year, bs="re") +# s(dayOfYear, year, bs="re") +
                           # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                           s(domainID, bs="re") + s(domainID, dayOfYear, bs="re") +
                           offset(log(numberTested))
                         # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                         , data=tck_borrelia_adj
                         , method="REML"
                         , family=binomial)
c(mod.gambin_nonlcdClass$aic - mod.gambin_domainnoelev$aic) # elevation still tells you stuff

mod.gambin_fixedDomain <- gam(borrPresent ~ s(dayOfYear) + #nlcdClass +
                               s(dayOfYear,plotID, bs="re") + s(plotID, bs="re") + 
                               s(year, bs="re") +# s(dayOfYear, year, bs="re") +
                               # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                               s(elevation) + domainID + 
                               offset(log(numberTested))
                             # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                             , data=tck_borrelia_adj
                             , method="REML"
                             , family=binomial)
c(mod.gambin_withDomain$aic - mod.gambin_fixedDomain$aic)

mod.gambin_fixeddomainnoelev <- gam(borrPresent ~ s(dayOfYear) + #nlcdClass +
                                     s(dayOfYear,plotID, bs="re") + s(plotID, bs="re") + 
                                     s(year, bs="re") +# s(dayOfYear, year, bs="re") +
                                     # s(logNLtckDensity) + s(lognymphDensity) + s(logadultDensity) + s(logtckDensity) +
                                     domainID +
                                     offset(log(numberTested))
                                   # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                                   , data=tck_borrelia_adj
                                   , method="REML"
                                   , family=binomial)
c(mod.gambin_fixedDomain$aic, mod.gambin_fixeddomainnoelev$aic, mod.gambin_withDomain$aic)

## FINAL MODEL:
mod.gambin<- gam(borrPresent ~ s(dayOfYear)  +
                                 s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
                                 s(year, bs="re") + s(dayOfYear, year, bs="re") +
                                 # s(elevation) +
                                domainID +
                                 offset(log(numberTested))
                               # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                               , data=tck_borrelia_adj
                               , method="REML"
                               , family=binomial)
gam.check(mod.gambin)
summary(mod.gambin)
plot(mod.gambin, scale=0, pages=1)
vis.gam(mod.gambin, view = c("elevation","dayOfYear"), theta=-45)
vis.gam(mod.gambin, view = c("year","dayOfYear"), theta=-45)
vis.gam(mod.gambin, view = c("plotID","dayOfYear"), theta=-45)

# Get the "predicted" presence/absence of borrellia, and see how "correct" it was?
plot(NULL,xlim=c(0,1),ylim=c(0,1), xlab=c("Probability Threshold"), ylab=c("Percentage Preditions correct"))
for ( p in seq(0,1, by = 0.05)) {
  points(x=p, y=mean((predict(mod.gambin, type = "response")>p) == tck_borrelia_adj$borrPresent))
}
# Proportion of correct predictions
mean((predict(mod.gambin, type = "response")>0.5) == tck_borrelia_adj$borrPresent)

# Of those where they predicted presence/absence incorrectly, how "far" were they off? 
# (i.e. were there are "negative" predictions that actually had a LOT of borrellia?)
tck_borrelia_adj %>%
  mutate(pred =as.numeric(predict(mod.gambin, type="response"))
         ,correctPredBin = ((pred>0.5) == tck_borrelia_adj$borrPresent) ) %>%
  ggplot() +geom_jitter(aes(x=log(numberPositive+1), y=correctPredBin), height=0.25, width=0) + ylab("Correctly predicted")

tck_borrelia_adj %>%
  mutate(pred =as.numeric(plogis(predict(mod.gambin)))
         ,correctPredBin =((pred>0.5) == tck_borrelia_adj$borrPresent) ) %>%
  ggplot() +geom_point(aes(x=log(numberPositive+1), y=pred, col=correctPredBin))

# Ideally, you have a positive, linear correlation between proportion positive and predicted probability. "Errors" should be correlated with small sample size.
tck_borrelia_adj %>%
  mutate(pred =as.numeric(plogis(predict(mod.gambin)))
         ,correctPredBin =((pred>0.5) == tck_borrelia_adj$borrPresent) ) %>%
  ggplot() +geom_point(aes(x=proportionPositive, y=pred, col=correctPredBin, cex=log(numberTested)))

# What are the groups that they least accurately predicted?
tck_borrelia_adj %>%
  mutate(predictedProbability =as.numeric(plogis(predict(mod.gambin)))
         ,correctPredBin =((predictedProbability>0.5) == tck_borrelia_adj$borrPresent) 
         # , proportionPositive = ifelse(proportionPositive==0, NA, proportionPositive)
         ) %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=plotID, fill=proportionPositive, col=predictedProbability), pch=21,cex=3) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_gradient(low="white", high="darkred") +
  facet_grid(.~year)
# Is there a bias for predicting negative or positive results?
tck_borrelia_adj %>%
  mutate(predictedProbability =as.numeric(plogis(predict(mod.gambin)))
         , predictedPA = (predictedProbability>0.5)
         ,correctPredBin =((predictedProbability>0.5) == tck_borrelia_adj$borrPresent) ) %>%
  select(borrPresent,predictedPA) %>% table()
16/(190+16) # False positive rate
33/(69+33) # False negative rate
# The false negative rate is actually a lot higher than the false positive rate. 
# This means, on average, more samples are positive then you'd expect, given the model.

## So finally, we need to filter the dataset by the predicted "absent" and fit a poisson or negative binomial model
# To include the maximum amount of samples as possible, I'm going to filter the dataset sort of un-usually--
# First, I am going to keep all positive borrelia results. Then, for each "negative" borrelia result, 
# I will look at the predicted probability in our model and determine whether it is a "binomial/hurdle" negative, or a "poisson/NB" negative.

tck_borrelia_filtBin <- tck_borrelia_adj %>%
  mutate(pred=predict(mod.gambin, type="response")) %>%
  filter((borrPresent>0 | pred > 0.5))

# Double check I did the filtering correctly
tck_borrelia_filtBin %>%
  mutate(predictedPA = (pred > 0.5)) %>%
  select(borrPresent,predictedPA) %>% table()

# Check if there is overdispersion
test_glm_filt <- glm(numberPositive ~ factor(month) + elevation + plotID + year
                     , data=tck_borrelia_filtBin
                     , family="quasipoisson"
                     , offset = log(numberTested)
)
summary(test_glm_filt)

# Try a regular poisson, and see how many zeros it finds
test_glm_filt2 <- glm(numberPositive ~ factor(month) + elevation + plotID + year
                     , data=tck_borrelia_filtBin
                     , family="poisson"
                     , offset = log(numberTested)
)

summary(test_glm_filt2)
mu <- predict(test_glm_filt2, type="response") # estimated mean
exp <- sum(dpois(x=0, lambda=mu)) # get probability of zero, then add those up to get total zeros expected?
round(exp)
sum(tck_borrelia_filtBin$numberPositive==0) # real number of zeros
# Yes, there appears to be overdispersion... but the zero-inflation is not nearly as much!

# Now try a negative binomial fit, and see how many zeros it predicts
test_glm_filt3 <- glm(proportionPositive ~ factor(month) + elevation + plotID + year
                      , data=tck_borrelia_filtBin
                      , family=negative.binomial(theta = 1) 
                      , offset = log(numberTested)
)
summary(test_glm_filt3)
mu <- predict(test_glm_filt3, type="response") # estimated mean
sum(round(mu*tck_borrelia_filtBin$numberTested) == 0) # Number of zeros it estimates, given r for each sample
sum(tck_borrelia_filtBin$numberPositive==0) # real number of zeros
# The zero-inflation is significantly reduced!

# make year a factor
tck_borrelia_filtBin <- tck_borrelia_filtBin %>%
  mutate(year=factor(year))



#### Finally, the second GAM model for abundance ####

# We should do the same model selection as above
mod.gampois_all <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                      s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                      s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                      s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                        s(elevation) +
                        s(logtckDensity, sp=2) +s(logNLtckDensity, sp=2) 
                    , data=tck_borrelia_filtBin
                    , method="REML"
                    , family="poisson")
summary(mod.gampois_all)

mod.gampois_notckdensity <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                         s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                         s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                         s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                         s(elevation) 
                         # s(logtckDensity, sp=2) +s(logNLtckDensity, sp=2) 
                       , data=tck_borrelia_filtBin
                       , method="REML"
                       , family="poisson")
c(mod.gampois_all$aic - mod.gampois_notckdensity$aic) # fit is WAY worse when no tck density.

mod.gampois_withtckDensity <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                                  s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                                  s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                  s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                  s(elevation) +
                                s(logtckDensity, sp=2) #+s(logNLtckDensity, sp=2)
                                , data=tck_borrelia_filtBin
                                , method="REML"
                                , family="poisson")
mod.gampois_withNLDensity <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                                    s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                                    s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                    s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                    s(elevation) +
                                    s(logNLtckDensity, sp=2)
                                  , data=tck_borrelia_filtBin
                                  , method="REML"
                                  , family="poisson")
mod.gampois_withnymphDensity <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                                   s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                                   s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                   s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                   s(elevation) +
                                   s(lognymphDensity, sp=2)
                                 , data=tck_borrelia_filtBin
                                 , method="REML"
                                 , family="poisson")
c(mod.gampois_withtckDensity$aic - mod.gampois_notckdensity$aic)  
c(mod.gampois_withNLDensity$aic - mod.gampois_notckdensity$aic) # exact same effect as total tck density; keep NLtckdensity
c(mod.gampois_withnymphDensity$aic - mod.gampois_notckdensity$aic) 

mod.gampois_noelevation <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                                      s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                                      s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                      s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                      # s(elevation) +
                                      s(logNLtckDensity, sp=2)
                                    , data=tck_borrelia_filtBin
                                    , method="REML"
                                    , family="poisson")
c(mod.gampois_withNLDensity$aic - mod.gampois_noelevation$aic)  # slightly better fit without elevation; can remove.

mod.gampois_noyear <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                                 s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                                 s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                 # s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                 # s(elevation) +
                                 s(logNLtckDensity, sp=2)
                               , data=tck_borrelia_filtBin
                               , method="REML"
                               , family="poisson")
c(mod.gampois_noelevation$aic - mod.gampois_noyear$aic) # need year

mod.gampois_noplot <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                            s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                            # s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                            s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                            # s(elevation) +
                            s(logNLtckDensity, sp=2)
                          , data=tck_borrelia_filtBin
                          , method="REML"
                          , family="poisson")
c(mod.gampois_noelevation$aic - mod.gampois_noplot$aic) # need plot

mod.gampois_nonlcdClass <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                            s(dayOfYear, sp=2)   +# nlcdClass +# Main effects: day of year and nymph/adult tick density
                            s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                            s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                            # s(elevation) +
                            s(logNLtckDensity, sp=2)
                          , data=tck_borrelia_filtBin
                          , method="REML"
                          , family="poisson")
c(mod.gampois_noelevation$aic - mod.gampois_nonlcdClass$aic) # need nlcdClass

mod.gampois_nodayOfYear<- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                                 # s(dayOfYear, sp=2)   +
                                 nlcdClass +# Main effects: day of year and nymph/adult tick density
                                 s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                 s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                 # s(elevation) +
                                 s(logNLtckDensity, sp=2)
                               , data=tck_borrelia_filtBin
                               , method="REML"
                               , family="poisson")
c(mod.gampois_noelevation$aic - mod.gampois_nodayOfYear$aic) # need dayOfYear

# Best fit model sofar:
summary(mod.gampois_noelevation)
# Do we need year (no interaction)?
mod.gampois_noyearalone<- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                                s(dayOfYear, sp=2)+nlcdClass +# Main effects: day of year and nymph/adult tick density
                                s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                s(year,dayOfYear, bs="re") + # Random effect of year
                                # s(elevation) +
                                s(logNLtckDensity, sp=2)
                              , data=tck_borrelia_filtBin
                              , method="REML"
                              , family="poisson")
c(mod.gampois_noelevation$aic - mod.gampois_noyearalone$aic)  #no; better if included.

#### FINAL MODEL:

mod.gam_pois <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
                  s(dayOfYear, sp=2) + s(logNLtckDensity, sp=2)  + nlcdClass +# Main effects: day of year and nymph/adult tick density
                  s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                  s(year, bs="re") + s(year,dayOfYear, bs="re") # Random effect of year
                , data=tck_borrelia_filtBin
                , method="REML"
                , family="poisson")
gam.check(mod.gam_pois)
summary(mod.gam_pois)
plot(mod.gam_pois, pages=1, scale=0)
plot(mod.gam_pois$residuals ~ mod.gam_pois$fitted.values, ylim=c(-max(abs(mod.gam_pois$residuals)),max(abs(mod.gam_pois$residuals)) ))
# residuals look sort of wonky
vis.gam(mod.gam_pois, view = c("logNLtckDensity","dayOfYear"), theta=45)
# sanity check to see that this relationship actually makes sense
plot(tck_borrelia_filtBin$proportionPositive ~ tck_borrelia_filtBin$logNLtckDensity)
plot(tck_borrelia_filtBin$proportionPositive ~ tck_borrelia_filtBin$dayOfYear)


### Try the whole model section process again with a binomial
mod.gambin2_all <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                         s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                         s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                         s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                         s(elevation) + domainID +
                         s(logtckDensity, sp=1) +s(logNLtckDensity, sp=1) 
                       , data=tck_borrelia_filtBin
                       , method="REML"
                       , family=binomial)
summary(mod.gambin2_all)

mod.gambin2_nodensity <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                         s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                         s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                         s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                         s(elevation) +domainID
                         # s(logtckDensity, sp=2) +s(logNLtckDensity, sp=2) 
                       , data=tck_borrelia_filtBin
                       , method="REML"
                       , family=binomial)
c(mod.gambin2_nodensity$aic, mod.gambin2_all$aic) # need density
mod.gambin2_logtckdensity <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                               s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                               s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                               s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                               s(elevation) +domainID +
                             s(logtckDensity, sp=1) 
                             , data=tck_borrelia_filtBin
                             , method="REML"
                             , family=binomial)
mod.gambin2_logNLdensity <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                                   s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                                   s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                   s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                   s(elevation) +domainID +
                                 s(logNLtckDensity, sp=1) 
                                 , data=tck_borrelia_filtBin
                                 , method="REML"
                                 , family=binomial)
c(mod.gambin2_all$aic, mod.gambin2_logNLdensity$aic, mod.gambin2_logtckdensity$aic)
plot(mod.gambin2_logNLdensity, select=7)
plot(mod.gambin2_logtckdensity, select=7)

# Let's just keep NLdensity?

mod.gambin2_noelevation <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                               s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                               s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                               s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                 domainID +# s(elevation) + 
                             s(logNLtckDensity, sp=1)
                             , data=tck_borrelia_filtBin
                             , method="REML"
                             , family=binomial)
c(mod.gambin2_noelevation$aic, mod.gambin2_logNLdensity$aic) #  better with elev; should keep

mod.gambin2_noyear <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                                 s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                                 s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                 # s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                 s(elevation) + domainID +
                                 s(logNLtckDensity, sp=1)
                               , data=tck_borrelia_filtBin
                               , method="REML"
                               , family=binomial)
c(mod.gambin2_logNLdensity$aic, mod.gambin2_noyear$aic) # better with year

mod.gambin2_noplot <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                            s(dayOfYear, sp=2)   + nlcdClass +# Main effects: day of year and nymph/adult tick density
                            # s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                            s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                            s(elevation) + domainID +
                            s(logNLtckDensity, sp=1)
                          , data=tck_borrelia_filtBin
                          , method="REML"
                          , family=binomial)
c(mod.gambin2_logNLdensity$aic, mod.gambin2_noplot$aic) # better with plot

mod.gambin2_nonlcdClass <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                            s(dayOfYear, sp=2)   + #nlcdClass +# Main effects: day of year and nymph/adult tick density
                            s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                            s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                            s(elevation) + domainID +
                            s(logNLtckDensity, sp=1)
                          , data=tck_borrelia_filtBin
                          , method="REML"
                          , family=binomial)
c(mod.gambin2_logNLdensity$aic, mod.gambin2_nonlcdClass$aic) # nlcd Class improves a bit

mod.gambin2_nodayofyear <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                                 # s(dayOfYear, sp=2)   + 
                                 nlcdClass +# Main effects: day of year and nymph/adult tick density
                                 s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                 s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                 s(elevation) + domainID +
                                 s(logNLtckDensity, sp=1)
                               , data=tck_borrelia_filtBin
                               , method="REML"
                               , family=binomial)
c(mod.gambin2_logNLdensity$aic, mod.gambin2_nodayofyear$aic) # better with dayofyear

mod.gambin2_nodomain <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                                 # s(dayOfYear, sp=2)   + 
                                 nlcdClass +# Main effects: day of year and nymph/adult tick density
                                 s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                 s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                 s(elevation) + domainID +
                                 s(logNLtckDensity, sp=1)
                               , data=tck_borrelia_filtBin
                               , method="REML"
                               , family=binomial)
c(mod.gambin2_logNLdensity$aic, mod.gambin2_nodayofyear$aic) # better with dayofyear

### BEST MODEL:
mod.gambin_2 <- gam(cbind(numberPositive,numberTested) ~ #offset(log(numberTested)) +  #Offset
                                 s(dayOfYear, sp=2)   +
                                 nlcdClass +# Main effects: day of year and nymph/adult tick density
                                 s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                                 s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
                                 s(elevation) + domainID +s(logNLtckDensity, sp=1)
                               , data=tck_borrelia_filtBin
                               , method="REML"
                               , family=binomial)
gam.check(mod.gambin_2)
summary(mod.gambin_2)
plot(mod.gambin_2, pages=1)
vis.gam(mod.gambin_2, view = c("dayOfYear","logNLtckDensity"), theta=45)
vis.gam(mod.gambin_2, view = c("dayOfYear","elevation"), theta=45)

# 
# #### GAM on non-zero counts ####
# # Filter to only non-zeros
# tck_borrelia_nozeros <- tck_borrelia_adj %>%
#   filter(numberPositive>0) %>%
#   mutate(year=factor(year))
# 
# mod.gam_pois_nozeros <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
#                       s(dayOfYear, sp=2) + s(logNLtckDensity, sp=1)  + nlcdClass +# Main effects: day of year and nymph/adult tick density
#                       s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
#                       s(year, bs="re") + s(year,dayOfYear, bs="re") # Random effect of year
#                     , data=tck_borrelia_nozeros
#                     , method="REML"
#                     , family="poisson")
# gam.check(mod.gam_pois_nozeros)
# summary(mod.gam_pois_nozeros)
# plot(mod.gam_pois_nozeros, pages=1, scale=0)
# plot(mod.gam_pois_nozeros$residuals ~ mod.gam_pois_nozeros$fitted.values, ylim=c(-max(abs(mod.gam_pois_nozeros$residuals)),max(abs(mod.gam_pois_nozeros$residuals)) ))
# # residuals look sort of wonky
# vis.gam(mod.gam_pois_nozeros, view = c("logNLtckDensity","dayOfYear"), theta=45)
# # sanity check to see that this relationship actually makes sense
# plot(tck_borrelia_nozeros$proportionPositive ~ tck_borrelia_nozeros$logNLtckDensity)
# plot(tck_borrelia_nozeros$proportionPositive ~ tck_borrelia_nozeros$dayOfYear)
# 
# ## The poisson model's residuals are wonky, and we've already established that there is likely over-dispersion. 
# # Therefore, let's fit a quasipoisson to see if this improves residual distribution.
# 
# mod.gam_qpois_nozeros <- gam(numberPositive ~ offset(log(numberTested)) +  #Offset
#                        s(dayOfYear, sp=2) + s(logNLtckDensity, sp=1)  + nlcdClass +# Main effects: day of year and nymph/adult tick density
#                        s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
#                        s(year, bs="re") + s(year,dayOfYear, bs="re") # Random effect of year
#                      , data=tck_borrelia_nozeros
#                      , method="REML"
#                      , family="quasipoisson")
# gam.check(mod.gam_qpois_nozeros)
# summary(mod.gam_qpois_nozeros)
# plot(mod.gam_qpois_nozeros, pages=1, scale=0)
# plot(mod.gam_qpois_nozeros$residuals ~ mod.gam_qpois_nozeros$fitted.values, ylim=c(-max(abs(mod.gam_qpois_nozeros$residuals)),max(abs(mod.gam_qpois_nozeros$residuals)) ))
# plot(mod.gam_qpois_nozeros$residuals ~ mod.gam_qpois_nozeros$fitted.values, ylim=c(-2,2) )
# # residuals look sort of wonky still-- and not symmetric...
# vis.gam(mod.gam_qpois_nozeros, view = c("logNLtckDensity","dayOfYear"), theta=45)
# 
# 
# ## Last option: a negative binomial.
# mod.gam_nb_nozeros <- gam(proportionPositive ~ offset(log(numberTested)) +  #Offset
#                     s(dayOfYear) + s(logNLtckDensity)  + nlcdClass +# Main effects: day of year and nymph/adult tick density
#                     s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
#                     s(year, bs="re") + s(year,dayOfYear, bs="re") # Random effect of year
#                   , data=tck_borrelia_nozeros
#                   , method="REML"
#                   , family=nb)
# gam.check(mod.gam_nb_nozeros)
# summary(mod.gam_nb_nozeros)
# plot(mod.gam_nb_nozeros, pages=1, scale=0)
# plot(mod.gam_nb_nozeros$residuals ~ mod.gam_nb_nozeros$fitted.values, ylim=c(-max(abs(mod.gam_nb_nozeros$residuals)),max(abs(mod.gam_nb_nozeros$residuals)) ))
# plot(mod.gam_nb_nozeros$residuals ~ mod.gam_nb_nozeros$fitted.values, ylim=c(-2,2) )
# # residuals look sort of wonky still-- and not symmetric...
# vis.gam(mod.gam_nb_nozeros, view = c("logNLtckDensity","dayOfYear"), theta=45)
# 

#### Can I use this in brms? ####
library("brms")

##### STOP HERE #####
#### TESTING BRMS AND GAMs #####
get_prior(bf(borrPresent ~ s(dayOfYear)  +
               s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
               s(year, bs="re") + s(dayOfYear, year, bs="re") +
               s(elevation) + domainID +
               offset(log(numberTested)))
             , data=tck_borrelia_adj
             , family=bernoulli
)
# 
# mod.brm_test1 <- brm(bf(borrPresent ~ s(dayOfYear)  +
#                           s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
#                           s(year, bs="re") + s(dayOfYear, year, bs="re") +
#                           s(elevation) +
#                           offset(log(numberTested)))
#           , data=tck_borrelia_adj
#           , family=bernoulli
# )
# summary(mod.brm_test1)
# plot(mod.brm_test1)
# conditional_smooths(mod.brm_test1, smooths = "s(dayOfYear)")
# conditional_smooths(mod.brm_test1, smooths = "s(elevation)")


# Hurdle component
# Let's use the best model from regular GAMs
summary(mod.gambin)
gambin_priors <- c(prior(student_t(3,-1.7, 5), class="Intercept" )
                   ,prior(normal(-1.21, 5), class="b", coef="domainIDD02" )
                   ,prior(normal(-5.86, 5), class="b", coef="domainIDD03" )
                   ,prior(normal(-1.27, 5), class="b", coef="domainIDD05" )
                   ,prior(normal(-6.14, 5), class="b", coef="domainIDD06" )
                   ,prior(normal(-5.70, 5), class="b", coef="domainIDD07" )
                   ,prior(normal(-6.35, 5), class="b", coef="domainIDD08" )
                   ,prior(normal(2.98, 5), class="b", coef="sdayOfYear_1" )
                   ,prior(normal(7.85, 5), class="b", coef="selevation_1" )
)

if (FALSE ) {
  brm_bin <- brm(bf(borrPresent ~ s(dayOfYear)  +
                      s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
                      s(year, bs="re") + s(dayOfYear, year, bs="re") +
                      s(elevation) + domainID +
                      offset(log(numberTested)))
                 , seed=48
                 , data=tck_borrelia_adj
                 , family=bernoulli
                 , prior = gambin_priors
                 , control=list(adapt_delta=0.99, max_treedepth=15)
                 )
  # brm_bin_withDomain <- brm(bf(borrPresent ~ s(dayOfYear)  +nlcdClass +
  #                     s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") +
  #                     s(year, bs="re") + s(dayOfYear, year, bs="re") +
  #                     s(elevation) +s(domainID, bs="re") +
  #                     offset(log(numberTested)))
  #                , seed=48
  #                , data=tck_borrelia_adj
  #                , family=bernoulli)

  save(brm_bin, file="brm_bin.RData")
  # save(brm_bin_noelev, file="brm_bin_noelev.RData")
  
} else {
  load("brm_bin.RData")
  # load("brm_bin_noelev.RData")
}

summary(brm_bin)
conditional_effects(brm_bin, effects = "dayOfYear")
conditional_effects(brm_bin, effects = "elevation") # elevation effect super weird

pp_check(brm_bin, nsamples = 100)
pp_check(brm_bin, type = "ecdf_overlay")

# Look at how estimated probability maps to actual data for borrelia
fitted(brm_bin) %>%
  cbind(brm_bin$data) %>%
  ggplot() +geom_point(aes(x=Estimate, y=borrPresent))

# Inspect error rate (false positive and false negative)
predict(brm_bin) %>%
  cbind(brm_bin$data) %>%
  mutate(Pred = ifelse(Estimate>0.5, 1, 0)) %>%
  select(borrPresent, Pred) %>% table()
11/(11+195) # false positive rate
29/(29+73) # false negative rate

# See if there is correlation between probability and how many positives there actually were
fitted(brm_bin) %>%
  cbind(tck_borrelia_adj[,c("proportionPositive","numberPositive", "numberTested")]) %>%
  ggplot() + geom_point(aes(x=numberPositive, y=Estimate, cex=numberTested)) 

## Look at effect of various predictors
# Effect of day of year
fitted(brm_bin) %>%
  cbind(brm_bin$data) %>%
  left_join(tck_borrelia_adj) %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=Estimate, col=nlcdClass)) + geom_point(aes(x=dayOfYear, y=borrPresent), col="red",alpha=0.2) + geom_smooth(aes(x=dayOfYear, y=borrPresent))
conditional_smooths(brm_bin, smooths="s(dayOfYear)")

# Effect of tck density
fitted(brm_bin) %>%
  cbind(brm_bin$data) %>%
  left_join(tck_borrelia_adj) %>%
  ggplot() + geom_point(aes(x=elevation, y=Estimate, col=domainID)) + geom_point(aes(x=elevation, y=borrPresent), col="red",alpha=0.2) + geom_smooth(aes(x=elevation, y=borrPresent)) 
conditional_smooths(brm_bin, smooths="s(elevation)")

# Since each sample is IID, we can include all positive results in poisson component of model, and use fitted probabilities to determine
# whether each zero is a "binomial" zero or a "poisson" zero.

tck_borrelia_filtbin_brm <- tck_borrelia_adj %>%
  cbind(predict(brm_bin)) %>%
  filter(borrPresent>0 | Estimate>0.5)

## Now try to fit a poisson distribution?
tck_borrelia_filtbin_brm %>%
  ggplot() + geom_histogram(aes(x=proportionPositive), bins=20)
tck_borrelia_filtbin_brm %>%
  ggplot() + geom_histogram(aes(x=numberPositive), bins=20)

### Binomial
brm_bin2 <- brm(bf(numberPositive | trials(numberTested) ~ #offset(log(numberTested)) +  #Offset
      s(dayOfYear)   +
      nlcdClass + # Main effects: day of year and nymph/adult tick density
      s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
      s(year, bs="re") + s(year,dayOfYear, bs="re") + # Random effect of year
      s(logNLtckDensity))
    , data=tck_borrelia_filtbin_brm
    , family=binomial()
    , seed=24
)

summary(brm_bin2)
plot(brm_bin2)
conditional_smooths(brm_bin2)
plot()
summary(mod.gambin_2)


## POISSON
# cbind with binomial; also, re-center
tck_borrelia_filtbin_brm <- tck_borrelia_filtbin_brm %>%
  mutate(dayOfYear_cen = (dayOfYear-mean(dayOfYear)/sd(dayOfYear))
         , logNLtckDensity_cen = (logNLtckDensity-mean(logNLtckDensity)/sd(logNLtckDensity))
         , logadultDensity_cen = (logadultDensity-mean(logadultDensity)/sd(logadultDensity))
         , lognymphDensity_cen = (lognymphDensity-mean(lognymphDensity)/sd(lognymphDensity))
         , year=factor(year))
brm_bin2 <- brm(numberPositive | trials(numberTested) ~ 
                  s(dayOfYear_cen) + s(logNLtckDensity_cen)  + nlcdClass +# Main effects: day of year and nymph/adult tick density
                  s(plotID, bs="re") + s(dayOfYear_cen,plotID, bs="re") + # Random Effect of plot
                  s(year, bs="re") + s(year,dayOfYear_cen, bs="re") # Random effect of year
                , data=tck_borrelia_filtbin_brm
                , family=binomial()
                , control = list(adapt_delta=0.99, max_treedepth=15)
                )
# try setting priors to limit how wiggly things can be-- might help with convergence
# play around with wigglyness in GAM and then set into brm
summary(brm_bin2)
conditional_smooths(brm_bin2)

conditional_smooths(brm_bin2, smooths="s(dayOfYear_cen)")
conditional_smooths(brm_bin2, smooths="s(logNLtckDensity_cen)")
posterior_predict()

predict(brm_bin2) %>%
  cbind(tck_borrelia_filtbin_brm[,c("numberPositive","dayOfYear","numberTested","logNLtckDensity","year","plotID")]) %>%
  ggplot() + geom_point(aes(x=log(numberPositive+1), y=log(Estimate+1))) +
  geom_segment(aes(x=log(numberPositive+1), xend=log(numberPositive+1), y=log(Q2.5+1), yend=log(Q97.5+1)), col="red")



brm_pois <- brm(numberPositive ~ offset(log(numberTested)) +  #Offset
                  s(dayOfYear) + s(logNLtckDensity)  + nlcdClass +# Main effects: day of year and nymph/adult tick density
                  s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                  s(year, bs="re") + s(year,dayOfYear, bs="re") # Random effect of year
                , data=tck_borrelia_filtbin_brm
                , family=poisson
                , control = list(adapt_delta=0.99, max_treedepth=15))
save(brm_pois, file="brms_pois.RData")
pp_check(brm_pois)

predict(brm_pois) %>%
  cbind(tck_borrelia_filtbin_brm[,c("numberPositive","dayOfYear","numberTested","logNLtckDensity","year","plotID")]) %>%
  ggplot() + geom_point(aes(x=log(numberPositive+1), y=log(Estimate+1))) +
  geom_segment(aes(x=log(numberPositive+1), xend=log(numberPositive+1), y=log(Q2.5+1), yend=log(Q97.5+1)), col="red")



# Trying to figure out how tin incorporate uncertainty into predictions (of number tested) and also
# trying to figure out which samples to keep vs discard for post-hurdle component of model



### Trying to figure out how to get confidence intervals base don sample size
dbinom(x=0, size=2, prob=0.8)

test_brm_hurdPois <- brm(bf(borrPresent ~ offset(log(numberTested)) +  #Offset
                              s(dayOfYear, sp=2) + s(logNLtckDensity, sp=1)  + nlcdClass +# Main effects: day of year and nymph/adult tick density
                              s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
                              s(year, bs="re") + s(year,dayOfYear, bs="re") )
                         , seed = 029348
                    # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                    , data=tck_borrelia_adj
                    , family=hurdle_poisson()
                    , control=list(adapt_delta=0.99))
## Note: 2665 divergent transitions
# 1335 exceeded maximum treedepth of 10
# Probably not enough data given the number of parameters we're testing...



test_brm <- brm(numberPositive ~ offset(log(numberTested)) +  #Offset
       s(dayOfYear) + s(logNLtckDensity)  + nlcdClass +# Main effects: day of year and nymph/adult tick density
       s(plotID, bs="re") + s(dayOfYear,plotID, bs="re") + # Random Effect of plot
       s(year, bs="re") + s(year,dayOfYear, bs="re") # Random effect of year
     , data=tck_borrelia_filtBin
     , family=poisson)
predict(test_brm) %>%
  cbind(tck_borrelia_filtBin[,c("numberPositive","dayOfYear","numberTested","logNLtckDensity","year","plotID")]) %>%
  ggplot() + geom_point(aes(x=log(numberPositive+1), y=log(Estimate+1))) +
  geom_segment(aes(x=log(numberPositive+1), xend=log(numberPositive+1), y=log(Q2.5+1), yend=log(Q97.5+1)), col="red")


# 
# 
# 
# #### Basic hurdle model #####
# # Now, let's try a simple hurdle model with just month and nNymph.
# library("pscl") # for hurdle
# mod.hurdle <- hurdle(numberPositive ~ 1
#        , data=tck_borrelia
#        , offset=log(numberTested)
#        , dist = "poisson"
#        )
# summary(mod.hurdle)
# predict(mod.hurdle)
# # Why does it predict a constant..? Beause no predictors?




# #### Using basic poisson model #####
# mod.pois1 <- glm(numberPositive ~ factor(month) + nlcdClass + log(NLtckDensity)
#                  , family = poisson (link="log")
#                  ,  data=tck_borrelia_adj
#                  , offset = log(numberTested)
# )
# summary(mod.pois1)
# mu <- exp(predict(mod.pois1, type="response")) # estimated mean
# exp <- sum(dpois(x=0, lambda=mu)) # get probability of zero, then add those up to get total zeros expected?
# round(exp)
# sum(tck_borrelia$numberPositive==0)
# # There's definitely underfitting of zeros.

## It turns out the offset needs to be on the same scale as the distribution! (link="log")
# 
# mod.hurdle1 <- hurdle(numberPositive ~  nlcdClass 
#                      , data=tck_borrelia
#                      , offset=log(numberTested)
#                      , dist = "poisson"
# )
# summary(mod.hurdle1)
# plot(tck_borrelia$numberPositive, predict(mod.hurdle1))
# plot((tck_borrelia$numberPositive-predict(mod.hurdle1))~ predict(mod.hurdle1))
# 
# 
# 
# mod.hurdle2 <- hurdle(numberPositive ~ nNymph + nlcdClass
#                       , data=tck_borrelia
#                       # , offset=numberTested
#                       , dist = "poisson"
# )
# summary(mod.hurdle2)
# plot(tck_borrelia$numberPositive, predict(mod.hurdle2))
# plot((tck_borrelia$numberPositive-predict(mod.hurdle2))~ predict(mod.hurdle2))

# The models keep failing when I try to use "complicated" models. Not sure what it is-- not enough data?
# 
# #### GAM ####
# 
# ### First thing: when alone, does the number of nymphs change the number of positive
# mod.gam <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) 
#     , data=tck_borrelia
#     , method="REML")
# summary(mod.gam)
# plot(mod.gam, pages=1)
# 
# mod.gam1 <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + nNymph
#                , data=tck_borrelia
#                , method="REML"
# )
# summary(mod.gam1)
# plot(mod.gam1, pages=1)
# 
# mod.gam1 <- gam(proportionPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + nNymph
#                 , data=tck_borrelia
#                 , method="REML"
# )
# summary(mod.gam1)
# plot(mod.gam1, pages=1)
# 
# tck_borrelia %>% ggplot() + geom_point(aes(x=log(numberTested), y=log(nNymph)))
# 
# mod.gam2 <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + nNymph 
#                 , data=tck_borrelia
#                 , method="REML")
# summary(mod.gam2)
# plot(mod.gam2, pages=1)
# plot(mod.gam2$residuals ~ mod.gam2$fitted.values)
# plot(mod.gam2$fitted.values ~ tck_borrelia$numberPositive, ylim=c(-200,200))
# # There's one outlier
# 
# mod.gam3 <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + nNymph + s(plotID, bs="re")
#                 , data=tck_borrelia
#                 , method="REML")
# summary(mod.gam3)
# plot(mod.gam3, pages=1)
# plot(mod.gam3$residuals ~ mod.gam3$fitted.values)
# plot(mod.gam3$fitted.values ~ tck_borrelia$numberPositive, ylim=c(-200,200))
# # There's one outlier
# 
# mod.gam3a <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + s(plotID, bs="re")
#                 , data=tck_borrelia
#                 , method="REML")
# summary(mod.gam3a)
# plot(mod.gam3a, pages=1)
# plot(mod.gam3a$residuals ~ mod.gam3a$fitted.values)
# plot(mod.gam3a$fitted.values ~ tck_borrelia$numberPositive, ylim=c(-200,200))
# # There's one outlier
# 
# mod.gam4 <- gam(numberPositive ~ offset(log(numberTested)) + s(dayOfYear) + ti(dayOfYear, nNymph) + s(plotID, bs="re") + nlcdClass
#                 , data=tck_borrelia
#                 , method="REML"
#                 , family=poisson(link="log"))
# gam.check(mod.gam4)
# summary(mod.gam4)
# plot(mod.gam4, pages=1)
# vis.gam(mod.gam4, view = c("nNymph","dayOfYear"))
# 
# mod.gam5 <- gam(proportionPositive ~ offset(log(numberTested)) + s(dayOfYear) + nNymph + s(plotID, bs="re") + nlcdClass
#                 , data=tck_borrelia
#                 , method="REML"
#                 , family=nb(link="log"))
# gam.check(mod.gam5)
# summary(mod.gam5)
# vis.gam(mod.gam5, view = c("nNymph","dayOfYear"), theta=-30)
# plot(mod.gam5, pages=1)

# 
# 
# 
# #### Binomial ######
# ### Hmmmm not working. Let's try a poor-man's hurdle model
# bin.intercept <- glm(borrPresent ~ 1
#     , data=tck_borrelia
#     , family="binomial")
# inv.logit(bin.intercept$coefficients)
# 
# bin.mod1 <- glm(borrPresent ~ -1 +factor(month) + nlcdClass + nNymph
#     , data=tck_borrelia
#     # , offset=numberTested
#     , family="binomial"
#     )
# summary(bin.mod1)
# Anova(bin.mod1, type = 3)
# 
# 
# # Now, take the non-zeros and predict
# tck_borrelia_nozeros <- tck_borrelia %>%
#   mutate(month=factor(month)) %>%
#   filter(numberPositive>0)
# # qpois.mod2 <- glm(numberPositive ~ -1 + factor(month) + nlcdClass + nNymph
# #                 , data=tck_borrelia_nozeros
# #                 , family="quasipoisson"
# #                 , offset=numberTested
# #                 ) #### Doesn't work-- dispersion is like 1e14!!
# qpois.mod2 <- glm(numberPositive ~ -1 + factor(month) + nlcdClass + nNymph
#                   , data=tck_borrelia_nozeros
#                   , family="quasipoisson"
#                   # , offset=numberTested
# )
# summary(qpois.mod2)
# exp(coef(qpois.mod2))
# Anova(qpois.mod2, type=3)
# 
# 
# tck_borrelia_nozeros %>%
#   ggplot() +geom_histogram(aes(x=log(numberPositive)))
# tck_borrelia_nozeros %>%
#   ggplot() +geom_histogram(aes(x=proportionPositive))
# 
# # See if the distribution looks right
# 
# cbind(tck_borrelia_nozeros, pred=predict(qpois.mod2)) %>%
#   ggplot() + geom_point(aes(x=(numberPositive), y=exp(pred), col=numberTested)) +
#   geom_abline(aes(intercept=0, slope=1))
# cbind(tck_borrelia_nozeros, pred=predict(qpois.mod2)) %>%
#   ggplot() + geom_histogram(aes(x=numberPositive))+geom_histogram(aes(x=exp(pred)), col="red", alpha=0.5)
# 
# 
# ### Adding random effects ####
# 
# ggnorm.mod3 <- lmer(log(numberPositive) ~ -1 + factor(month) + nlcdClass + (1|domainID)
#                 , data=tck_borrelia_nozeros
#                 # , family="poisson"
#                 # , offset=numberTested
# )
# summary(ggnorm.mod3)
# Anova(ggnorm.mod3)
# ranef(ggnorm.mod3)
# 
# 
# #### Trying GAM ####
# 
# gam1 <- gam(numberPositive ~ s(dayOfYear)+ s(nNymph)+ nlcdClass
#     ,data=tck_borrelia
#     , method = "REML"
#     , family = "quasipoisson")
# gam.check(gam1)
# summary(gam1)
# plot(gam1, shift=coef(gam1)[1], pages=1, seWithMean = TRUE, scheme = 2)
# # 
# # gam2 <- gam(numberPositive ~ s(dayOfYear) + nNymph+ ti(dayOfYear,nNymph)+nlcdClass
# #             ,data=tck_borrelia_fullzeroes
# #             , method = "REML"
# #             )
# # gam.check(gam2)
# # summary(gam2)
# # plot(gam2, shift=coef(gam2)[1], seWithMean = TRUE, scheme=1, pages=1)
# # vis.gam(gam2, view=c("dayOfYear","nNymph"), too.far = 0.05, theta=45)
# 
# gam3 <- gam(numberPositive ~ s(dayOfYear) + s(nNymph)+ ti(dayOfYear,nNymph) + nlcdClass
#             ,data=tck_borrelia_fullzeroes
#             , method = "REML"
# )
# gam.check(gam3)
# summary(gam3)
# plot(gam3, shift=coef(gam3)[1], seWithMean = TRUE, scheme=1, pages=1)
# vis.gam(gam3, view=c("dayOfYear","nNymph"))
# 
# 
# gam4 <- gam(proportionPositive ~ ti(dayOfYear,nNymph) + nlcdClass + s(dayOfYear) + nNymph
#             ,data=tck_borrelia_fullzeroes
#             , method = "REML"
# )
# summary(gam4)
# plot(gam4, shift=coef(gam4)[1], seWithMean = TRUE, scheme=1, pages=1)
# vis.gam(gam4, view=c("dayOfYear","nNymph"))
# 
# gam5 <- gam(proportionPositive ~ ti(dayOfYear,nNymph) + nlcdClass + s(dayOfYear) + s(nNymph)
#             ,data=tck_borrelia_fullzeroes
#             , method = "REML"
# )
# summary(gam5)
# plot(gam5, shift=coef(gam5)[1], seWithMean = TRUE, scheme=1, pages=1)
# vis.gam(gam5, view=c("dayOfYear","nNymph"))
# 
# 
# tck_borrelia_nozeros %>%
#   ggplot(aes(x=dayOfYear)) +geom_point(aes(y=numberPositive)) + geom_point(aes(y=proportionPositive*400), col="red") +
#   scale_y_continuous(name="Number positive for Borrelia", sec.axis=sec_axis(~./400, name="Proportion Positive")) +
#   geom_point(aes(y=log(nNymph)*50), col="blue") +
#   xlim(0,365)
# 
# # Let's try to simulate real data?
# rqpois <- function(n, mu, theta) {
#   rnbinom(n = n, mu = mu, size = mu/(theta-1))
# }
# 
# nsim <- 1000
# bin_coeff <- bin.mod1$coefficients
# bin_pred <- predict(bin.mod1, newdata=tck_borrelia_nozeros,type="response")
# pois_coeff <- c(0,bin.mod2$coefficients, 0, 0, 0)
# 
# var_overdisp <- 37.73
# # x.pred <- tck_borrelia %>%
# #   select(month)
# x.pred <- predict(bin.mod2, type="response")
# pred.mat <- matrix(ncol=length(x.pred), nrow=nsim, dimnames = list(seq(1:nsim), x.pred))
# for ( n in 1:nsim ) {
#   # Binomial
#   # temp.bin <- rbinom(n=length(x.pred), size=1, prob= exp(bin_coeff))
#   # temp.X <- pois_coeff[as.character(x.pred)]
#   # temp.pois <- rnorm(n=length(temp.X), mean=temp.X, sd=temp.X*sqrt(var_overdisp))
#   # temp.Y <- temp.bin*exp(temp.pois)
#   # pred.mat[n,] <- temp.Y
#   # Try neg binomial?
#   pred.mat[n,] <- rqpois(n=length(x.pred), mu=x.pred, theta=var_overdisp)
# }
# 
# pred.mat.long <- pred.mat %>% 
#   as_tibble() %>%
#   gather(key="Predicted", value="Simulated") %>%
#   mutate(Predicted=as.numeric(Predicted), Observed= rep(tck_borrelia_nozeros$numberPositive, each=nsim))
# 
# pred.mat.long %>%
#   ggplot() + geom_point(aes(x=log(Observed), y=log(Simulated)), alpha=0.5)
# 
# 
# cbind(tck_borrelia_nozeros, pred = inv.logit(predict(bin.mod2))) %>%
#   ggplot() +
#   geom_line(aes(x=month, y=predPropPos)) +
#   geom_point(aes(x=month, y=propPositive))
  