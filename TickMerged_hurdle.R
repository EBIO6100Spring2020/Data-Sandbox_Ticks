######## GAM model with manual hurdle component ##########
library(tidyverse)
library(mgcv) # for GAMs
library(randomForest) # for random forest

tck <- read.csv("data_derived/MASTER_all_tck_data_merged.csv")

#### Filter tick dataset and adjust ####

## How many plots were tested, and how many times?
tck %>%
  mutate(tested=ifelse(is.na(testingID),"no","yes")) %>%
  select(plotID, tested) %>% table(useNA="ifany")
## What type of ticks did they test?
tck %>%
  mutate(tested=ifelse(is.na(testingID),"no","yes")) %>%
  select(lifeStage, tested) %>% table(useNA="ifany") 

## Aggregate by plotID:dayOfYear:year. Count number of different life stages; create proportiontested/positive etc
tck_allsamples_borr <- tck %>%
  mutate(tested=ifelse(is.na(Borrelia_sp.),0,1) # Was each tick  ever tested for borrellia?
         ,isPositive=ifelse(Borrelia_sp.=="Positive", 1,0) # And if it was tested, was it positive or negative?
  ) %>%
  group_by(domainID, siteID, nlcdClass, plotID, elevation, collectDate, dayOfYear, year, month, totalSampledArea) %>% # collapse by sample-- summing all tck counts together
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
         , tckDensity = sum(c(nLarva, nNymph, nAdult))/totalSampledArea
         , ntckDensity = nNymph/totalSampledArea
         , NLtckDensity = sum(c(nNymph, nAdult), na.rm = TRUE)/totalSampledArea
         , atckDensity = nAdult/totalSampledArea)  %>%
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

## Let's filter out plots that NEVER had any positive borrelia
infectedPlots <- tck %>%
  filter(Borrelia_sp.=="Positive") %>%
  select(plotID) %>%pull() %>% unique()

## Filtering same data as above; only now removing non-infected plots.
tck_borrelia_positivePlots <- tck_allsamples_borr %>%
  filter(plotID %in% as.character(infectedPlots))

## Re-assess the apperance of plots
tck_borrelia_positivePlots %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=plotID, fill=proportionPositive, col=testingStatus), pch=21) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_manual(values=c(Tested="blue", `Nymphs present, not tested`="grey", `No nymphs`="white")) +
  facet_grid(.~year)

## There are STILL a lot of zeros-- zeros where there are nymphs but they were never tested.
# Let's remove those.
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

## Inspecting possible predictors of borrelia
# dayOfYear vs borrelia
tck_borrelia_adj %>%
  ggplot() +geom_point(aes(x=dayOfYear, y=numberPositive))+
  geom_point(aes(x=dayOfYear, y=proportionPositive*max(tck_borrelia_adj$numberPositive)), col="red", alpha=0.75) +
  scale_y_continuous(sec.axis=sec_axis(trans=~./max(tck_borrelia_adj$numberPositive), name="Proportion positive"))+
  theme(axis.title.y.right = element_text(color="red")) + ylab("Number positive") 
# elevation vs borrelia
tck_borrelia_adj %>%
  ggplot() +geom_point(aes(x=elevation, y=proportionPositive, col=domainID))

# nymph density vs borrelia
tck_borrelia_adj %>%
  ggplot() +geom_point(aes(x=lognymphDensity, y=proportionPositive)) +
  geom_point(aes(x=lognymphDensity, y=borrPresent), col="red")

# adult density vs borrelia
tck_borrelia_adj %>%
  ggplot() +geom_point(aes(x=logadultDensity, y=proportionPositive)) +
  geom_point(aes(x=logadultDensity, y=borrPresent), col="red")

# nymph and adult density vs borrelia-- first, are they correlated?
tck_borrelia_adj %>%
  ggplot() + geom_density2d(aes(x=logadultDensity, y=lognymphDensity))
tck_borrelia_adj %>%
  ggplot() +geom_point(aes(x=logNLtckDensity, y=proportionPositive)) +
  geom_point(aes(x=logNLtckDensity, y=borrPresent), col="red")

#### Fitting a two-step GAM ####

stepGAM <- function(dat, predictors, response, constant_pred=NULL, family, ignore.combos) {
  ## For troubleshooting
  # dat <- tck_borrelia_adj
  # predictors = allPred
  # response = "borrPresent"
  # constant_pred <- "offset(log(numberTested))"
  # family=binomial
  # ignore.combos=list(c("domainID","s(domainID, bs='re')"), c("s(logNLtckDensity)","s(logadultDensity)","s(lognymphDensity)")
  #                    , c("s(logNLtckDensity)","s(logadultDensity)"),c("s(logNLtckDensity)","s(lognymphDensity)"))
  n_rows <- sum(sapply(1:length(predictors),FUN=function(x){ncol(combn(predictors,m=x))} ))
  allAIC <- setNames(data.frame(matrix(nrow=n_rows, ncol=(4+length(predictors)))),nm = c("AIC","REML","Dev.expl","formula",predictors))
  i <- 1
  pb <- txtProgressBar(title = "progress bar", min = 0,
                       max = n_rows, style=3)
  for ( p in 1:length(predictors) ) {
    combos <- combn(predictors, m=p)
    for ( c in 1:ncol(combos) ) {
      # skip any combinations that are meant to be ignored
      if (any(sapply(ignore.combos, FUN=function(ic){sum(ic %in% combos[,c]) == length(ic)}))) {
        i = i + 1
        next
      }
      if ( length(constant_pred)>0 ) {
        frml.temp <- paste(c(paste0(response," ~ ",constant_pred),combos[,c]), collapse="+")
      } else {
        frml.temp <- paste0(paste0(response," ~ "),paste(combos[,c], collapse=" + "), sep="")
      }
      # If more coefficients than data
      res <- try(gam(as.formula(frml.temp)
              , data=dat
              , method="REML"
              , family=family))
      if ( inherits(res, "try-error") ) {
        next
      } else {
        gam.temp <- gam(as.formula(frml.temp)
                        , data=dat
                        , method="REML"
                        , family=family)
        
        pred.temp <- (colnames(allAIC) %in% combos[,c])[-c(1,2,3,4)]
        allAIC[i,c("AIC","REML","Dev.expl","formula")]  <- as.vector(c(gam.temp$aic,gam.temp$gcv.ubre[1],summary(gam.temp)$dev.expl,frml.temp))
        allAIC[i, predictors] <- pred.temp
        
      }
      
      i <- i+1
      Sys.sleep(0.1)
      setTxtProgressBar(pb, i, label=paste( round(i/total*100, 0),
                                            "% done"))
    }
  }
  close(pb)
  return(allAIC)
}
#### Part I: binomial hurdle component ####
# First, model a binomial component with dayOfYear, nlcdClass, and non-larval tick density (NLtckDensity) as predictors of borrellia presence/absence.
# Plot is a random effect, which I include as a varying-intercept, varying-slope model.
allPred <- c("s(dayOfYear)","nlcdClass","s(elevation)"
             ,"s(logNLtckDensity)","s(logadultDensity)","s(lognymphDensity)"
             ,"s(plotID, bs='re')", "s(plotID, dayOfYear, bs='re')"
             ,"s(year, bs='re')", "s(year, dayOfYear, bs='re')"
             , "domainID", "s(domainID, bs='re')"
)


if (FALSE) {
  allAIC <- stepGAM(dat = tck_borrelia_adj, predictors = allPred, response = "borrPresent", constant_pred = "offset(log(numberTested))", family = binomial, ignore.combos=list(c("domainID","s(domainID, bs='re')"), c("s(logNLtckDensity)","s(logadultDensity)","s(lognymphDensity)")                                                                                                                                                                 , c("s(logNLtckDensity)","s(logadultDensity)"),c("s(logNLtckDensity)","s(lognymphDensity)")))
  allAIC_filt <- allAIC %>% filter(!is.na(AIC)) %>% arrange(AIC) %>%mutate(rank=seq(1:length(AIC))) %>% 
    mutate(AIC=as.numeric(AIC)
           , Dev.expl = as.numeric(Dev.expl)
           , REML = as.numeric(REML)) %>%
    filter(AIC<500)
  save(allAIC, file="allAIC.RData")
  save(allAIC_filt, file="allAIC.RData")
  
} else {
  load("allAIC.RData")
  load("allAIC_filt.RData")
  
}
allAIC_filt
# allAIC %>%
#   as_tibble() %>%
#   arrange(AIC) %>% View()

nrow(allAIC_filt)
# Let's see distribution of AIC values
allAIC_filt %>% 
  ggplot() +geom_line(aes(x=rank, y=AIC))+
  geom_line(aes(x=rank, y=Dev.expl*(max(allAIC_filt$AIC))), col="red", alpha=0.2) + scale_y_continuous(sec.axis=sec_axis(~./(max(allAIC_filt$AIC)), name="Deviance Explained")) +
  theme(axis.title.y.right = element_text(colour = "red"))

# Compare AIC for models with and without each predictor
allAIC_filt %>% gather(-c(AIC, REML, Dev.expl, formula,rank), key=Predictor, value=WithPredictor) %>% 
  select(AIC, REML, Dev.expl, Predictor, WithPredictor) %>%
  mutate(WithPredictor = ifelse(WithPredictor==0, FALSE, TRUE)) %>%
  ggplot() + geom_violin(aes(x=WithPredictor, y=AIC))  +  facet_wrap(.~Predictor) 

# Compare deviance explained for models with and without each predictor
allAIC_filt %>% gather(-c(AIC, REML, Dev.expl, formula,rank), key=Predictor, value=WithPredictor) %>% 
  select(AIC, REML, Dev.expl, Predictor, WithPredictor) %>%
  mutate(WithPredictor = ifelse(WithPredictor==0, FALSE, TRUE)) %>%
  ggplot() + geom_violin(aes(x=WithPredictor, y=Dev.expl))  +  facet_wrap(.~Predictor) 


frml_bin1_bestAIC <- allAIC_filt[allAIC_filt$AIC==min(allAIC_filt$AIC),"formula"]
frml_bin1_bestDevexpl <- allAIC_filt[allAIC_filt$Dev.expl==max(allAIC_filt$Dev.expl),"formula"]

# Look at the model with the best AIC value
mod.gambin_bestAIC <- gam(as.formula(frml_bin1_bestAIC)
                      # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                      , data=tck_borrelia_adj
                      , method="REML"
                      , family=binomial)
summary(mod.gambin_bestAIC)
gam.check(mod.gambin_bestAIC)
plot(mod.gambin_bestAIC, pages=1)

# Look at the model with the best deviance explained
mod.gambin_bestDevExpl <- gam(as.formula(frml_bin1_bestDevexpl)
                          # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                          , data=tck_borrelia_adj
                          , method="REML"
                          , family=binomial)
summary(mod.gambin_bestDevExpl)
gam.check(mod.gambin_bestDevExpl)
plot(mod.gambin_bestDevExpl, pages=1)

# Compare accuracy rate
mean((predict(mod.gambin_bestAIC, type = "response")>0.5) == tck_borrelia_adj$borrPresent)
mean((predict(mod.gambin_bestDevExpl, type = "response")>0.5) == tck_borrelia_adj$borrPresent)

# Honestly, the "deviance explained" looks better and makes more sense, in my opinion.

## FINAL MODEL:
frml.bin1 <- frml_bin1_bestDevexpl
mod.gambin<- gam(as.formula(frml.bin1)
                 # , offset=log(numberTested) # I include the offset in the formula instead, so that it is included in predictions. Including it here does NOT incorporate numberTested in predictiosn.
                 , data=tck_borrelia_adj
                 , method="REML"
                 , family=binomial)
gam.check(mod.gambin)
summary(mod.gambin)
plot(mod.gambin, scale=0, pages=1)
vis.gam(mod.gambin, view = c("logNLtckDensity","dayOfYear"), theta=45)
vis.gam(mod.gambin, view = c("year","dayOfYear"), theta=45)
vis.gam(mod.gambin, view = c("plotID","dayOfYear"), theta=45)

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
17/(189+17) # False positive rate
29/(73+29) # False negative rate
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

# make year a factor
tck_borrelia_filtBin <- tck_borrelia_filtBin %>%
  mutate(year=factor(year))
#### Part II: the second GAM model for abundance ####

### Try the whole model section process again with a binomial
# I'm increasing spline smoothness here, because the AIC values change drastically when I don't include them
allPred2 <- c("s(dayOfYear)","nlcdClass","s(elevation)"
             ,"s(logNLtckDensity, sp=1)","s(logadultDensity, sp=1)","s(lognymphDensity, sp=1)"
             ,"s(plotID, bs='re')", "s(plotID, dayOfYear, bs='re')"
             ,"s(year, bs='re')", "s(year, dayOfYear, bs='re')"
             , "domainID", "s(domainID, bs='re')"
)
if ( FALSE ) {
  allAIC_2 <- stepGAM(dat=tck_borrelia_filtBin, predictors = allPred2, response = "cbind(numberPositive,numberTested)", family = binomial, ignore.combos=list(c("domainID","s(domainID, bs='re')")
                                                                                                                                                              , c("s(logNLtckDensity, sp=1)","s(logadultDensity, sp=1)","s(lognymphDensity, sp=1)")
                                                                                                                                                              , c("s(logNLtckDensity, sp=1)","s(lognymphDensity, sp=1)")
                                                                                                                                                              , c("s(logNLtckDensity, sp=1)","s(logadultDensity, sp=1)")
                                                                                                                                                              ))
  allAIC_2_filt <- allAIC_2 %>% filter(!is.na(AIC)) %>% arrange(AIC) %>%mutate(rank=seq(1:length(AIC))) %>% 
    mutate(AIC=as.numeric(AIC)
           , Dev.expl = as.numeric(Dev.expl)
           , REML = as.numeric(REML)) %>%
    filter(AIC<2000)
  save(allAIC_2, file="allAIC_2.RData")
  save(allAIC_2_filt, file="allAIC_2_filt.RData")
  } else {
  load("allAIC_2.RData")
  load("allAIC_2_filt.RData")
}
summary(allAIC_2)

allAIC_2%>%
  arrange(AIC) %>%View()

# Let's see distribution of AIC values
allAIC_2_filt %>% 
  ggplot() +geom_line(aes(x=rank, y=AIC))+
  geom_line(aes(x=rank, y=Dev.expl*(max(allAIC_2_filt$AIC))), col="red", alpha=0.2) + scale_y_continuous(sec.axis=sec_axis(~./(max(allAIC_2_filt$AIC)), name="Deviance Explained")) +
  theme(axis.title.y.right = element_text(colour = "red"))

# Compare AIC for models with and without each predictor
allAIC_2_filt %>% gather(-c(AIC, REML, Dev.expl, formula, rank), key=Predictor, value=WithPredictor) %>%
  select(AIC, REML, Dev.expl, Predictor, WithPredictor) %>%
  mutate(WithPredictor = ifelse(WithPredictor==0, FALSE, TRUE)) %>%
  ggplot() + geom_violin(aes(x=WithPredictor, y=AIC))  +  facet_wrap(.~Predictor)

# Compare deviance explained for models with and without each predictor
allAIC_2_filt %>% gather(-c(AIC, REML, Dev.expl, formula,rank), key=Predictor, value=WithPredictor) %>%
  select(AIC, REML, Dev.expl, Predictor, WithPredictor) %>%
  mutate(WithPredictor = ifelse(WithPredictor==0, FALSE, TRUE)) %>%
  ggplot() + geom_violin(aes(x=WithPredictor, y=Dev.expl))  +  facet_wrap(.~Predictor)


frml_bin2_bestAIC <- allAIC_2_filt[allAIC_2_filt$AIC==min(allAIC_2_filt$AIC),"formula"]
frml_bin2_bestDevexpl <- allAIC_2_filt[allAIC_2_filt$Dev.expl==max(allAIC_2_filt$Dev.expl),"formula"]
frml_bin2_bestAIC==frml_bin2_bestDevexpl
# frml_bin2_bestAIC_adjdayOfYear <- "cbind(numberPositive,numberTested) ~ s(dayOfYear, sp=1) + nlcdClass + s(logNLtckDensity, sp=1) + s(plotID, bs='re') + s(plotID, dayOfYear, bs='re') + s(year, bs='re') + s(year, dayOfYear, bs='re') + domainID"
# they're the same


## Best AIC/Dev expl model
mod.gambin2_bestAIC <- gam(as.formula(frml_bin2_bestAIC)
                           , data=tck_borrelia_filtBin
                           , method="REML"
                           , family=binomial)
summary(mod.gambin2_bestAIC)
plot(mod.gambin2_bestAIC, pages=1)
gam.check(mod.gambin2_bestAIC)

## Best adjusted model to limit dayOfYear
mod.gambin2_bestAIC_adj <- gam(as.formula(frml_bin2_bestAIC_adjdayOfYear)
                       , data=tck_borrelia_filtBin
                       , method="REML"
                       , family=binomial)
summary(mod.gambin2_bestAIC_adj)
plot(mod.gambin2_bestAIC_adj, pages=1)
gam.check(mod.gambin2_bestAIC_adj)

### BEST MODEL:
mod.gambin_2 <- gam(as.formula(frml_bin2_bestDevexpl)
                    , data=tck_borrelia_filtBin
                    , method="REML"
                    , family=binomial)
summary(mod.gambin_2)
plot(mod.gambin_2, pages=1)
gam.check(mod.gambin_2)
vis.gam(mod.gambin_2, view = c("dayOfYear","logNLtckDensity"), theta=45)
vis.gam(mod.gambin_2, view = c("year","dayOfYear"), theta=45)

library("brms")

#### BRMS: Part I #####
if (FALSE ) {
  brm_bin1 <- brm(bf(frml.bin1)
                 , seed=48
                 , data=tck_borrelia_adj
                 , family=bernoulli
                 # , prior = gambin_priors
                 , control=list(adapt_delta=0.99, max_treedepth=15)
                 )
  save(brm_bin1, file="brm_bin1.RData")
  # save(brm_bin_noelev, file="brm_bin_noelev.RData")
  
} else {
  load("brm_bin1.RData")
  # load("brm_bin_noelev.RData")
}
summary(brm_bin1)
plot(brm_bin1)
conditional_effects(brm_bin1, effects = "dayOfYear")
conditional_effects(brm_bin1, effects = "lognymphDensity")


# Look at how estimated probability maps to actual data for borrelia
fitted(brm_bin1) %>%
  cbind(brm_bin1$data) %>%
  ggplot() +geom_point(aes(x=Estimate, y=borrPresent))

# Inspect error rate (false positive and false negative)
predict(brm_bin1) %>%
  cbind(brm_bin1$data) %>%
  mutate(Pred = ifelse(Estimate>0.5, 1, 0)) %>%
  select(borrPresent, Pred) %>% table()
12/(12+194) # false positive rate
29/(29+73) # false negative rate

# See if there is correlation between probability and how many positives there actually were
fitted(brm_bin1) %>%
  cbind(tck_borrelia_adj[,c("proportionPositive","numberPositive", "numberTested")]) %>%
  ggplot() + geom_point(aes(x=numberPositive, y=Estimate, cex=numberTested)) 

## Look at effect of various predictors
# Effect of day of year
fitted(brm_bin1) %>%
  cbind(brm_bin1$data) %>%
  left_join(tck_borrelia_adj) %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=Estimate, col=nlcdClass)) + geom_point(aes(x=dayOfYear, y=borrPresent), col="red",alpha=0.2) + geom_smooth(aes(x=dayOfYear, y=borrPresent))
conditional_smooths(brm_bin1, smooths="s(dayOfYear)")

# Effect of tck density
fitted(brm_bin1) %>%
  cbind(brm_bin1$data) %>%
  left_join(tck_borrelia_adj) %>%
  ggplot() + geom_point(aes(x=lognymphDensity, y=Estimate, col=domainID)) + geom_point(aes(x=lognymphDensity, y=borrPresent), col="red",alpha=0.2) + geom_smooth(aes(x=lognymphDensity, y=borrPresent)) 
conditional_smooths(brm_bin1, smooths="s(lognymphDensity)")

# Since each sample is IID, we can include all positive results in poisson component of model, and use fitted probabilities to determine
# whether each zero is a "binomial" zero or a "poisson" zero.
tck_borrelia_filtbin_brm <- tck_borrelia_adj %>%
  cbind(predict(brm_bin1)) %>%
  filter(borrPresent>0 | Estimate>0.5)

## Now try to fit a second distribution?
tck_borrelia_filtbin_brm %>%
  ggplot() + geom_histogram(aes(x=proportionPositive), bins=20)
tck_borrelia_filtbin_brm %>%
  ggplot() + geom_histogram(aes(x=numberPositive), bins=20)

##### BRMS: Part II ####
# Re-format frml for brm
frml_bin2_adj  <- gsub(pattern="cbind(numberPositive,numberTested)",replacement="numberPositive | trials(numberTested)",frml_bin2_bestDevexpl, fixed = TRUE)

if (FALSE) {
  brm_bin2 <- brm(bf(as.formula(frml_bin2_adj))
                  , data=tck_borrelia_filtbin_brm
                  , family=binomial()
                  , seed=24
                  , control = list(adapt_delta=0.95, max_treedepth=15)
  )
  save(brm_bin2, file = "brm_bin2.RData")
} else {
  load("brm_bin2.RData")
}
summary(brm_bin2)

plot(brm_bin2)
conditional_smooths(brm_bin2, smooths = "s(dayOfYear)")
conditional_smooths(brm_bin2, smooths = "s(lognymphDensity, sp = 1)")
conditional_smooths(brm_bin2, smooths = "s(logadultDensity, sp = 1)")

conditional_effects(brm_bin2, effects = "dayOfYear")
conditional_effects(brm_bin2, effects = "logadultDensity")
conditional_effects(brm_bin2, effects = "lognymphDensity")


predict(brm_bin2) %>%
  cbind(tck_borrelia_filtbin_brm[,c("numberPositive","dayOfYear","numberTested","logNLtckDensity","year","plotID")]) %>%
  ggplot() + geom_point(aes(x=log(numberPositive+1), y=log(Estimate+1))) +
  geom_segment(aes(x=log(numberPositive+1), xend=log(numberPositive+1), y=log(Q2.5+1), yend=log(Q97.5+1)), col="red")








  