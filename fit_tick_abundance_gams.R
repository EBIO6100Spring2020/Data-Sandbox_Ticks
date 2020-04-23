### Fitting Tick Abundance Data using GAMs
### Wynne Moss
### April 21 2020
# This script uses GAMs in the package mgcv to predict tick abundances at NEON sites
# Data were previously cleaned and a test/train dataset was made
# The training dataset was additionally cleaned up by removing empty domains and scaling variables
# This script will fit a variety of zero-inflated poisson models and evaluate model fit
library(mgcv)
library(dplyr)
library(broom)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(caret)
# data used to fit model
abun <- readRDS("data_derived/nymph_abun_all.Rdata") # data used to fit models (optional)

####  1: READ IN MODELS/FIT MODELS ####
# check if models have already been run
models.exist <- file.exists("Model_objects/zipmodels_tick_abun.Rdata")
# if so, read them in and can skip the next part
if(models.exist==TRUE){
  readRDS("Model_objects/zipmodels_tick_abun.Rdata") -> zipmodlist
  list2env(zipmodlist, .GlobalEnv)
  rm(zipmodlist)
}

# if the models haven't been pulled/run this will run them
# careful as this may take a while
if(models.exist == FALSE){
  abun <- readRDS("data_derived/nymph_abun_all.Rdata")
  head(abun)
  # fit zero-inflated Poissons
  # first part of the model: poisson process (abundance)
  # second part of the model: zero inflation process (prob of 0)
  # a null model
  zip0 <- gam(list(estimatedCount ~ offset(logSampledArea),
                   ~1),
              data = abun, family = ziplss())
  
  # random effects only
  zip1 <- gam(list(estimatedCount ~ s(plotID, bs = "re") + offset(logSampledArea) ,
                   ~s(plotID, bs = "re")),
              data = abun, family = ziplss())
  
  # add in season as a predictor of presence
  zip2 <- gam(list(estimatedCount ~ s(plotID, bs = "re") + offset(logSampledArea) ,
                   ~ season + s(plotID, bs = "re")),
              data = abun, family = ziplss())
  
  # add in nlcd as a predictor of presence
  zip3 <- gam(list(estimatedCount ~s(plotID, bs = "re") + offset(logSampledArea) ,
                   ~ season + nlcdClass + s(plotID, bs = "re")),
              data = abun, family = ziplss())
  
  # add in day of year as a predictor of abundance
  zip4 <- gam(list(estimatedCount ~ s(sDayofYear) + s(plotID, bs = "re") + offset(logSampledArea),
                   ~season + nlcdClass + s(plotID, bs = "re")),
              data = abun, family = ziplss())
  
  # add year as a random effect 
  zip5 <- gam(list(estimatedCount ~ s(sDayofYear) + s(fYear, bs = "re") + s(plotID, bs = "re") + offset(logSampledArea) ,
                   ~season + nlcdClass + s(plotID, bs = "re")),
              data = abun, family = ziplss())
  
  # add year as a random effect NESTED witin plot 
  zip6 <- gam(list(estimatedCount ~ s(sDayofYear) +offset(logSampledArea) + s(yearPlot, bs = "re"),
                   ~season + nlcdClass + s(plotID, bs = "re")),
              data = abun, family = ziplss())
  
  pattern <- ls(pattern="zip") # get names of models
  # pattern <- grep("zip", names(.GlobalEnv), value = TRUE)
  models <- do.call("list", mget(pattern))
  
  # save models
  saveRDS(models, "Model_objects/zipmodels_tick_abun.Rdata" )
}
####  2: COMPARE MODELS WITH AIC ####
# compare models with AIC
pattern <- ls(pattern="zip") # get names of models
zipmodlist <- do.call("list", mget(pattern))
do.call(rbind, lapply(zipmodlist, glance)) %>% data.frame() %>%
  mutate(model = pattern) %>% arrange(AIC)

# take a look at "best" model via AIC
summary(zip6)
####  3: FIT "BEST" MODEL WITH BRMS ####
# now that I have the best model, fit in brms so I can better extract uncertainty
# warning that this takes quite a while! 
# if you want to run this, manually change the below statement to yes!
run.mod <- "no"
if(run.mod == "yes"){
  require(brms)
  abun <- readRDS("data_derived/nymph_abun_all.Rdata")
  b_zip6 <- brm(bf(
    estimatedCount ~  s(sDayofYear)+ (1|yearPlot) + offset(logSampledArea), # observation-level RE
    zi ~ nlcdClass + season + (1|plotID)),
    data = abun, family = zero_inflated_poisson(), iter =1000, chains = 4)
  b_zip7 <- update(b_zip6, iter = 2000) # run chains longer
  saveRDS(b_zip7, "Model_objects/b_zip7.Rdata")
}

# if you don't want to run it can just load it
b_zip7 <- readRDS("Model_objects/b_zip7.Rdata")
summary(b_zip7)
plot(b_zip7)

# I don't think this model worked well.


####  4: EVALUATE MODEL FIT USING TEST DATASET ####
nymphtest <- readRDS("data_derived/nymph_abun_test.Rdata")
head(nymphtest)
nymph_pred <- nymphtest %>% select(sDayofYear, yearPlot, logSampledArea, plotID, season, nlcdClass)
# are there any ranef NOT in the training set?
levels(nymph_pred$plotID)[which(!levels(nymph_pred$plotID) %in% levels(abun$plotID))]
# nope
levels(nymph_pred$yearPlot)[which(!levels(nymph_pred$yearPlot) %in% levels(abun$yearPlot))]
# nope

# first I need to figure out how the predict function works ............ 
predict.test <- data.frame(predict(zip6, newdata = nymph_pred, type = "link", se.fit = TRUE))
predict.test <- cbind(predict.test, predict(zip6, newdata = nymph_pred, type = "response", se.fit = TRUE))
colnames(predict.test) <- c("log_count", "logit_prob", "se_log_count", "se_logit_prob", "pred_count", "se_pred_count")
predict.test$pred_count_poisson <- exp(predict.test$log_count)
predict.test$pred_prob <- plogis(predict.test$logit_prob)
predict.test$pred_mean_count <- predict.test$pred_count*predict.test$pred_prob
head(predict.test)
# predicts log-count if present (which we exponentiate), as well as logit prob present
# type = response predicts expcount*prob which is the OVERALL expected count 

# now try this but with a different offset -- see what changes
nymph_pred2 <- nymph_pred %>% mutate(logSampledArea = 1) # sampling a smaller area
predicts2 <- data.frame(predict(zip6, newdata = nymph_pred2, type = "link", se.fit = TRUE))
predicts2 <- cbind(predicts2, predict(zip6, newdata = nymph_pred2, type = "response", se.fit = TRUE))
colnames(predicts2) <- c("log_count", "logit_prob", "se_log_count", "se_logit_prob", "pred_count", "se_pred_count")
predicts2$pred_count_poisson <- exp(predicts2$log_count)
predicts2$pred_prob <- plogis(predicts2$logit_prob)
predicts2$pred_mean_count <- predicts2$pred_count*predicts2$pred_prob
head(predicts2)
head(predict.test)
# when I lower the logsampled area, the predicted count is lower
# the predicted_mean_count is also lower.
# so, it is using the log sampled area in making predictions
# when we give it an offset of 0 it should be predicting DENSITY (count per UNIT area)
rm(nymph_pred2, predicts2, predict.test)

# I think we want the type response which is MEAN count
nymph_pred %>% mutate(logSampledArea=0) -> nymph_pred0
predict.count.nooff <- data.frame(predict(zip6, newdata = nymph_pred0, type = "response", se.fit = TRUE)) %>%
  select(predicted_count_nooff = fit, se.predicted_count_noof= se.fit) 

predict.count <- data.frame(predict(zip6, newdata = nymph_pred, type = "response", se.fit = TRUE)) %>%
  select(predicted_count = fit, se.predicted_count = se.fit) 
predict.pois <- data.frame(predict(zip6, newdata = nymph_pred, se.fit = TRUE)) %>% 
  select(predicted_logpois_count = fit.1, predicted_logodds = fit.2, 
         se_predicted_pois_count = se.fit.1, se_predicted_logodds = se.fit.2) %>%
  mutate(predicted_pois_count = exp(predicted_logpois_count), predicted_prob = plogis(predicted_logodds))

# add in prob prediction
nymphtest$predicted_density <- predict.count.nooff$predicted_count_nooff
nymphtest$predicted_density_se <- predict.count.nooff$se.predicted_count_noof
nymphtest$predicted_prob <- predict.pois$predicted_prob
nymphtest$predicted_logodds <- predict.pois$predicted_logodds
nymphtest$preidcted_logodds_se <- predict.pois$se_predicted_logodds
nymphtest$predicted_count <- predict.pois$predicted_pois_count

# add in observed values
nymphtest$observed_presence <- nymphtest$nymph_presence
nymphtest$observed_count <- nymphtest$estimatedCount
nymphtest$observed_dens <- nymphtest$density

rm(predict.count, predict.pois, predict.count.nooff, nymph_pred0)

nymphtest %>% ggplot() + geom_histogram(aes(predicted_prob, fill = as.factor(observed_presence))) + facet_wrap(~observed_presence) + theme_classic() +
  xlab("Model predicted probability of presence") +
  labs(fill = "Observed Tick Presence")

# plots where model predicts no ticks but they were observed?
nymphtest %>% filter(observed_presence==1, predicted_prob  < 0.05)
# only 3! pretty good.
nymphtest %>% filter(observed_presence==0, predicted_prob  > 0.6)

# plots where model predicts lots of ticks, how many actually had them?
# when model predicts > 50% chance of ticks, how many of these records have them?
nymphtest %>% filter(predicted_prob > 0.5) %>% pull(observed_presence) %>% table()
nymphtest %>% filter(predicted_prob > 0.5) %>% pull(predicted_prob) %>% mean()

# when model predicts < 50% chance of ticks, how many records have them?
nymphtest %>% filter(predicted_prob < 0.5) %>% pull(observed_presence) %>% table()

nymphtest %>% ggplot(aes(x = predicted_prob, y = observed_presence)) + geom_jitter(width = .001, height = 0.01, alpha = .3, size = 3, col = "forestgreen") +
  xlab("Model predicted probability of presence")+
  ylab("Observed tick presence") +
  theme_classic()

head(nymphtest)
ggplot(nymphtest) +
  geom_line(aes(x=sDayofYear, y = predicted_density, group = yearPlot, color= yearPlot))+
  facet_wrap(~nlcdClass)+
  theme(legend.position = "none") 

# for plots where nymphs were present, how close was the count?
nymphtest %>% filter(observed_presence == 1) %>%
  ggplot(aes(y=log(observed_dens), x = log(predicted_density)))+
  geom_point(alpha = .5, size = 3, col = "magenta")+
  geom_abline(slope = 1, intercept = 0, col ="magenta")+
  xlim(-2,2)+
  ylim(-4,2)+
  theme_classic()+
  xlab("Predicted density from model")+
  ylab("Observed density")

nymphtest %>% filter(observed_presence == 1) %>%
  ggplot(aes(y=log(observed_count), x = log(predicted_count)))+
  geom_point(alpha = .5, size = 3, col = "magenta")+
  geom_abline(slope = 1, intercept = 0, col ="magenta")+
  theme_classic()+
  xlab("Predicted count from model")+
  ylab("Observed count")

nymphtest %>% filter(observed_presence == 1) %>%
  ggplot(aes(y=log(observed_count), x = log(predicted_count*predicted_prob)))+
  geom_point(alpha = .5, size = 3, col = "magenta")+
  geom_abline(slope = 1, intercept = 0, col ="magenta")+
  theme_classic()+
  xlab("Predicted count from model")+
  ylab("Observed count")

nymphtest %>%
  ggplot(aes(y=observed_count, x = predicted_count*predicted_prob))+
  geom_point(alpha = .5, size = 3, col = "magenta")+
  geom_abline(slope = 1, intercept = 0, col ="magenta")+
  theme_classic()+
  xlab("Predicted count from model")+
  ylab("Observed count")

nymphtest %>%
  ggplot(aes(x=observed_count, y= predicted_count*predicted_prob))+
  geom_point(alpha = .5, size = 3, col = "magenta")+
  geom_abline(slope = 1, intercept = 0, col ="magenta")+
  theme_classic()+
  ylab("Predicted count from model")+
  xlab("Observed count")

# what are the sites with really bad predictions?
nymphtest$resid <- nymphtest$observed_count - nymphtest$predicted_prob*nymphtest$predicted_count
hist(nymphtest$resid)
nymphtest %>% filter(resid < -25) -> bad_preds

abun %>% filter(yearPlot %in% bad_preds$yearPlot) %>%
  ggplot(aes(x=sDayofYear, y = estimatedCount)) +
  geom_point()+
  facet_wrap(~plotID) +
  geom_point(data = bad_preds, aes(x=sDayofYear, y = observed_count), col = "red")+
  geom_point(data = bad_preds, aes(x=sDayofYear, y = predicted_count*predicted_prob), col = "blue")
# why are these so bad??

nymphtest %>% filter(resid>-25) -> good_preds
abun %>% filter(yearPlot %in% good_preds$yearPlot) %>%
  ggplot(aes(x=sDayofYear, y = estimatedCount)) +
  geom_point()+
  facet_wrap(~plotID) +
  geom_point(data = good_preds, aes(x=sDayofYear, y = observed_count), col = "red")+
  geom_point(data = good_preds, aes(x=sDayofYear, y = predicted_count*predicted_prob), col = "blue")


# our predictions are often way off... 
# is 2018 weird? 
nymph_train <- readRDS("data_derived/nymph_abun_all.Rdata")
nymph_test <- readRDS("data_derived/nymph_abun_test.Rdata")
nymph_all <- rbind(nymph_train, nymph_test)
nymph_all %>% ggplot(aes(x=sDayofYear, y = density, group = plotID)) +
  geom_point(aes(color = plotID), show.legend=FALSE)+
  facet_wrap(~year)
rm(nymph_train, nymph_test)
# yes! ha. we picked a bad dataset to train/test on.

# let's look at how we did for plots that aren't outliers!
nymph_all %>% group_by(plotID) %>% mutate(meanTicks = mean(density)) %>% ungroup()-> nymph_all
nymph_all %>% group_by(plotID) %>% summarise(meanTicks = mean(density)) -> plotMeans
hist(plotMeans$meanTicks)
plotMeans$plotID[which(plotMeans$meanTicks> mean(plotMeans$meanTicks)+2*sd(plotMeans$meanTicks))] -> outliers


nymphtest %>% filter(resid>-25) -> good_preds

abun %>% filter(!plotID %in% outliers) %>%
  ggplot(aes(x=sDayofYear, y = estimatedCount)) +
  geom_point()+
  facet_wrap(~plotID) +
  geom_point(data = nymphtest, aes(x=sDayofYear, y = observed_count), col = "red")+
  geom_point(data = nymphtest, aes(x=sDayofYear, y = predicted_count*predicted_prob), col = "blue")

nymphtest %>% filter(!plotID %in% outliers) %>%
  ggplot(aes(x=predicted_count*predicted_prob, y = observed_count, color = plotID), show.legend = FALSE)+
  geom_point() 

rm(bad_preds, good_preds, outliers, plotMeans)
####  5: EVALUATE MODEL PERFORMS USING FULL DATASET #########

# our model doesn't do a great job of predicting the test dataset, but how well does it fit the actual data?
pred_all <- data.frame(predict(zip6, type = "response", se.fit=TRUE)) %>% select(pred_count = fit, pred_count_se = se.fit)
pred_all_p <- data.frame(predict(zip6)) %>% select(pred_logcount = X1, pred_logprob = X2) %>% 
  mutate(pred_count_p = exp(pred_logcount)*plogis(pred_logprob))

pred_all <- cbind(abun, pred_all, pred_all_p)

# just do the first 25
pred_all %>% filter(plotID %in% levels(plotID)[1:25]) %>%
ggplot(aes(x=sDayofYear, group = yearPlot)) +
  geom_line(aes(y = pred_count, x = sDayofYear, group = yearPlot, color = as.factor(year)))+
  facet_wrap(~plotID, scales = "free_y") +
  geom_point(aes(y = estimatedCount, x = sDayofYear, group = yearPlot, color = as.factor(year)))

pred_all %>% filter(plotID %in% levels(plotID)[26:50]) %>%
  ggplot(aes(x=sDayofYear, group = yearPlot)) +
  geom_line(aes(y = pred_count, x = sDayofYear, group = yearPlot, color = as.factor(year)))+
  facet_wrap(~plotID, scales = "free_y") +
  geom_point(aes(y = estimatedCount, x = sDayofYear, group = yearPlot, color = as.factor(year)))

pred_all %>% filter(plotID %in% levels(plotID)[51:75]) %>%
  ggplot(aes(x=sDayofYear, group = yearPlot)) +
  geom_line(aes(y = pred_count, x = sDayofYear, group = yearPlot, color = as.factor(year)))+
  facet_wrap(~plotID, scales = "free_y") +
  geom_point(aes(y = estimatedCount, x = sDayofYear, group = yearPlot, color = as.factor(year)))

pred_all %>% filter(plotID %in% levels(plotID)[76:100]) %>%
  ggplot(aes(x=sDayofYear, group = yearPlot)) +
  geom_line(aes(y = pred_count, x = sDayofYear, group = yearPlot, color = as.factor(year)))+
  facet_wrap(~plotID, scales = "free_y") +
  geom_point(aes(y = estimatedCount, x = sDayofYear, group = yearPlot, color = as.factor(year)))

pred_all %>% filter(plotID %in% levels(plotID)[76:100]) %>%
  ggplot(aes(x=sDayofYear, group = yearPlot)) +
  geom_line(aes(y = pred_count_p, x = sDayofYear, group = yearPlot, color = as.factor(year)))+ # switch to other prediction
  facet_wrap(~plotID, scales = "free_y") +
  geom_point(aes(y = estimatedCount, x = sDayofYear, group = yearPlot, color = as.factor(year)))

pred_all$pred_logprob

confusionMatrix(data = as.factor(as.numeric(plogis(pred_all$pred_logprob)>0.5)),
                      reference = as.factor(pred_all$nymph_presence))
# model does pretty well with predicting pres/abs
confusionMatrix(data= as.factor(as.numeric(nymphtest$predicted_prob>0.5)),
                reference = as.factor(nymphtest$nymph_presence))                
# this is true even for test dataset


####  6: EVALUATE HOW A NULL MODEL PERFORMS ####
if(file.exists("Model_objects/zip_null_re.Rds")==FALSE){
  zipre <- gam(list(estimatedCount ~  s(yearPlot, bs = "re"),
                    ~ s(plotID, bs = "re")),
               offset = logSampledArea,
               data = abun, family = ziplss())
  saveRDS(zipre, "Model_objects/zip_null_re.Rds")
} else {zipre <- readRDS("Model_objects/zip_null_re.Rds")}

predicts_null_re <- data.frame("pred_count" = predict(zipre, type = "response")) 
predicts_null_re <- cbind(abun, predicts_null_re)
predicts_null_p <- data.frame(predict(zipre)) %>% mutate(predicted_prob = plogis(X2), predicted_poiscount = exp(X1)) %>%
  select(predicted_prob, predicted_poiscount)
predicts_null_re <- cbind(predicts_null_re, predicts_null_p)
rm(predicts_null_p)

# of the plots where there were no nymphs, what was the predicted prob?

predicts_null_re %>% ggplot() + geom_histogram(aes(predicted_prob, fill = as.factor(nymph_presence))) + facet_wrap(~nymph_presence) + theme_classic() +
  xlab("Model predicted probability of presence") +
  labs(fill = "Observed Tick Presence")

# plots where model predicts no ticks but they were observed?
predicts_null_re %>% filter(nymph_presence==1, predicted_prob  < 0.05) # many from unusual nlcds
predicts_null_re %>% filter(nymph_presence==0, predicted_prob  > 0.50) # all from TREE! 


# a lot more now

predicts_null_re %>% ggplot(aes(x = predicted_prob, y = nymph_presence)) + geom_jitter(width = .001, height = 0.01, alpha = .1, size = 3, col = "forestgreen") +
  xlab("Model predicted probability of presence")+
  ylab("Observed tick presence") +
  theme_classic()


# for plots where nymphs were present, how close was the count?
predicts_null_re %>% filter(nymph_presence == 1) %>%
  ggplot(aes(y=log(estimatedCount), x = log(predicted_poiscount*predicted_prob)))+
  geom_point(alpha = .5, size = 3, col = "magenta")+
  geom_abline(slope = 1, intercept = 0, col ="magenta")+
  xlim(-3,6)+
  ylim(-3,6)+
  theme_classic()+
  xlab("log-Predicted count from model")+
  ylab("log-Observed count")

predicts_null_re %>% filter(nymph_presence == 1) %>%
  ggplot(aes(y=estimatedCount, x = predicted_poiscount*predicted_prob))+
  geom_point(alpha = .5, size = 3, col = "magenta")+
  theme_classic()+
  xlab("log-Predicted count from model")+
  ylab("log-Observed count")


predicts_null_re %>% filter(nymph_presence == 1) %>%
  ggplot(aes(y=log(estimatedCount), x = log(predicted_poiscount)))+
  geom_point(alpha = .5, size = 3, col = "magenta")+
  geom_abline(slope = 1, intercept = 0, col ="magenta")+
  xlim(-3,6)+
  ylim(-3,6)+
  theme_classic()+
  xlab("Predicted log-count from model")+
  ylab("Observed log-count")



####  7: INTERPRET MODEL COEFFICIENTS (IGNORE ALL THIS FOR NOW) #####
# let's only predict in plot-years where we have data
plotyearclasscombos <- unique(paste(abun$plotID, abun$yearPlot, abun$nlcdClass, sep ="_"))

newdata <- expand.grid("plotyearclasscombos" = plotyearclasscombos,
                       "sDayofYear" = seq(-2,2, by = .1),
                       "season" = c(0,1),
                       "logSampledArea" = 0) %>% 
  separate(plotyearclasscombos, sep = c(8,9,22,23), remove = FALSE,
                           into = c("plotID", "toss", "yearPlot", "toss2", "nlcdClass")) %>%
  select(-toss, -toss2, -plotyearclasscombos)  %>%
  separate(yearPlot, sep= 9, into = c("toss", "year"), remove = FALSE) %>% select(-toss) %>%
  mutate(year = as.factor(year)) %>%
  filter(!str_detect(as.character(yearPlot), "2019")) -> newdata
# get rid of season-day of year combos that don't exist
abun$year <- as.factor(abun$year)
abun %>% filter(season ==1) %>% summarise(minDay = min(sDayofYear), maxDay = max(sDayofYear)) -> tickSeason
# tick season goes from - 1.13 to 1.09
# remove days > maxDay AND season 1
newdata %>% filter(!(season==1 & sDayofYear>tickSeason$maxDay))->newdata
newdata %>% filter(!(season ==1 & sDayofYear<tickSeason$minDay)) ->newdata
newdata %>% filter(!(season ==0 & sDayofYear>tickSeason$minDay & sDayofYear < tickSeason$maxDay)) ->newdata

# get a common color palette
myColors <- brewer.pal(6,"Set1")
names(myColors) <- levels(newdata$year)
colScale <- scale_colour_manual(name = "year",values = myColors)
fillScale <- scale_fill_manual(name = "year",values = myColors)

# predict for new data
summary(zip6o)
colnames(newdata)

preds <- data.frame(predict(zip6o, newdata=newdata, se.fit = TRUE))
colnames(preds) <- c("log_count", "logit_prob", "log_count_se", "logit_prob_se")
preds <- cbind(preds, newdata)
preds %>% mutate(
  conf.logit.low =logit_prob - 2*logit_prob_se,
  conf.logit.hi = logit_prob + 2*logit_prob_se,
  conf.logc.low = log_count -2*log_count_se,
  conf.logc.hi = log_count + 2*log_count_se
) %>% mutate(
  prob.low =plogis(conf.logit.low),
  prob.hi = plogis(conf.logit.hi),
  prob = plogis(logit_prob),
  count.low = exp(conf.logc.low),
  count.hi = exp(conf.logc.hi),
  count = exp(log_count)
) -> preds
head(preds)



# pick a few plots to illustrate
preds %>% ggplot(aes(x=sDayofYear, y = count)) +
  geom_point(aes(color = year))+
  facet_wrap(~nlcdClass)+
  colScale
abun %>% ggplot(aes(x=sDayofYear, y = estimatedCount))+
  geom_point(aes(color = year))+
  facet_wrap(~nlcdClass)+
  colScale

randplot <- sample(unique(newdata$plotID), size = 10, replace = FALSE)
abun %>% filter(plotID %in% randplot) -> abun_rand

preds %>% filter(plotID %in% randplot) %>% droplevels() %>%
ggplot(aes(x=sDayofYear, y = count, group = year))+
  geom_line(aes(color = year))+
  geom_ribbon(aes(ymin=count.low, ymax=count.hi, color = year, fill = year), alpha = .2)+
  facet_wrap(~plotID) +
  colScale+
  fillScale

# pick just one site
preds %>% filter(plotID == "BLAN_005") %>% ggplot(aes(x=sDayofYear, y = count, group = year)) +
  geom_line(aes(color = year))+
  geom_ribbon(aes(color = year, ymin = count.low, ymax = count.hi, fill = year),
                  alpha = .2)+
  geom_point(data = filter(abun, plotID == "BLAN_005"), aes(y=estimatedCount/totalSampledArea, group = year, color = year))+
  colScale+
  fillScale +
  xlab("Scaled day of year")+
  ylab("Tick Abundance")


head(preds)
newdata <- cbind(newdata, preds)
head(newdata)
ggplot(newdata,aes(x = sDayofYear, group = yearPlot)) +
  geom_point(aes(y = fit))+
  facet_wrap(~nlcdClass)

# random effects 
# ranef values is the uncertainty on a new location 
# prediction interval should be huge 
# simulate in the mgcv package
# AUC confusing matrix
# future directions

