######## using a hurdle model ##########
library(tidyverse)
library(pscl) # hurdle models
library(lme4) # mixed effect models
library(boot) # logit/inv logit
library(car) # Anova (type III)
library(mgcv) # for GAMs
library(randomForest) # for random forest

tck <- read.csv("data_derived/MASTER_all_tck_data_merged.csv")

#### Filter tick dataset and adjust ####
tck %>%
  mutate(tested=ifelse(is.na(testingID),"no","yes")) %>%
  select(plotID, tested) %>% table(useNA="ifany")


# Plot over time, by year, by plot, to see how zero-inflated it is.
tck_allsamples_borr <- tck %>%
  mutate(tested=ifelse(is.na(Borrelia_sp.),0,1)
         ,isPositive=ifelse(Borrelia_sp.=="Positive", 1,0)
  ) %>%
  group_by(domainID, siteID, nlcdClass, plotID, elevation, collectDate, dayOfYear, year, month) %>%
  summarize(numberTested=sum(tested, na.rm=TRUE)
            ,n=n()
            , numberPositive=sum(isPositive,na.rm = TRUE)
            , nAdult=sum(lifeStage=="Adult")
            , nNymph=sum(lifeStage=="Nymph")
            , nLarva=sum(lifeStage=="Larva")) %>%
  ungroup() %>%
  mutate(proportionTested=numberTested/nNymph
         , proportionPositive=numberPositive/numberTested
         , tested = ifelse(numberTested > 0 , TRUE, FALSE)
         , testingStatus = ifelse(numberTested > 0, "Tested", ifelse(nNymph>0, "Nymphs present, not tested", "No nymphs"))
         ) %>%
  mutate(nlcdClass=factor(nlcdClass, levels=c("emergentHerbaceousWetlands","cultivatedCrops","pastureHay","grasslandHerbaceous"
                                              ,"dwarfScrub","shrubScrub","sedgeHerbaceous"
                                              ,"woodyWetlands","deciduousForest","evergreenForest","mixedForest"))) %>%
  mutate(borrPresent = ifelse(numberPositive>0,1,0)
         , domainID = factor(as.character(domainID))
         , siteID = factor(as.character(siteID))
         , plotID = factor(as.character(plotID))) 

tck_allsamples_borr %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=plotID, fill=proportionPositive, col=testingStatus), pch=21) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_manual(values=c(Tested="blue", `Nymphs present, not tested`="grey", `No nymphs`="white")) +
  facet_grid(.~year)


# Let's also filter by plots that NEVER had any positive borrelia
infectedPlots <- tck %>%
  filter(Borrelia_sp.=="Positive") %>%
  select(plotID) %>%pull() %>% unique()

tck_borrelia_fullzeroes <- tck %>%
  filter(plotID %in% as.character(infectedPlots)) %>%
  mutate(tested=ifelse(is.na(Borrelia_sp.),0,1)
         ,isPositive=ifelse(Borrelia_sp.=="Positive", 1,0)
  ) %>%
  group_by(domainID, siteID, nlcdClass, plotID, elevation, collectDate, dayOfYear, year, month) %>%
  summarize(numberTested=sum(tested, na.rm=TRUE)
            ,n=n()
            , numberPositive=sum(isPositive,na.rm = TRUE)
            , nAdult=sum(lifeStage=="Adult")
            , nNymph=sum(lifeStage=="Nymph")
            , nLarva=sum(lifeStage=="Larva")) %>%
  ungroup() %>%
  mutate(proportionTested=numberTested/nNymph
         , proportionPositive=numberPositive/numberTested
         , tested = ifelse(numberTested > 0 , TRUE, FALSE)
         , testingStatus = ifelse(numberTested > 0, "Tested", ifelse(nNymph>0, "Nymphs present, not tested", "No nymphs"))) %>%
  mutate(nlcdClass=factor(nlcdClass, levels=c("emergentHerbaceousWetlands","cultivatedCrops","pastureHay","grasslandHerbaceous"
                                              ,"dwarfScrub","shrubScrub","sedgeHerbaceous"
                                              ,"woodyWetlands","deciduousForest","evergreenForest","mixedForest"))) %>%
  mutate(borrPresent = ifelse(numberPositive>0,1,0)
         , domainID = factor(as.character(domainID))
         , siteID = factor(as.character(siteID))
         , plotID = factor(as.character(plotID))) 

tck_borrelia_fullzeroes %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=plotID, fill=proportionPositive, col=testingStatus), pch=21) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_manual(values=c(Tested="blue", `Nymphs present, not tested`="grey", `No nymphs`="white")) +
  facet_grid(.~year)


tck_borrelia <- tck %>%
  filter(plotID %in% as.character(infectedPlots)) %>%
  mutate(tested=ifelse(is.na(Borrelia_sp.),0,1)
         ,isPositive=ifelse(Borrelia_sp.=="Positive", 1,0)
  ) %>%
  group_by(domainID, siteID, nlcdClass, plotID, elevation, collectDate, dayOfYear, year, month) %>%
  summarize(numberTested=sum(tested, na.rm=TRUE)
            , tested = ifelse(numberTested > 0, TRUE, FALSE)
            ,n=n()
            , numberPositive=sum(isPositive,na.rm = TRUE)
            , nAdult=sum(lifeStage=="Adult")
            , nNymph=sum(lifeStage=="Nymph")
            , nLarva=sum(lifeStage=="Larva")) %>%
  ungroup() %>%
  mutate(proportionTested=numberTested/nNymph
         , proportionPositive=numberPositive/numberTested
         , tested = ifelse(numberTested > 0 , TRUE, FALSE)
         , testingStatus = ifelse(numberTested > 0, "Tested", ifelse(nNymph>0, "Nymphs present, not tested", "No nymphs"))) %>%
  mutate(nlcdClass=factor(nlcdClass, levels=c("emergentHerbaceousWetlands","cultivatedCrops","pastureHay","grasslandHerbaceous"
                                                ,"dwarfScrub","shrubScrub","sedgeHerbaceous"
                                                ,"woodyWetlands","deciduousForest","evergreenForest","mixedForest"))) %>%
  mutate(borrPresent = ifelse(numberPositive>0,1,0)
         , domainID = factor(as.character(domainID))
         , siteID = factor(as.character(siteID))
         , plotID = factor(as.character(plotID))) %>%
  mutate(numberPositive=ifelse(numberTested==0,NA,numberPositive)) %>%
  filter(numberTested!=0)

tck_borrelia %>%
  ggplot() + geom_point(aes(x=dayOfYear, y=plotID, fill=proportionPositive, col=testingStatus), pch=21) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_manual(values=c(Tested="blue", `Nymphs present, not tested`="grey", `No nymphs`="white")) +
  facet_grid(.~year)

#### Preliminary Plotting ####

# First, what is the distribution of Borrelia prevalence?
tck_borrelia %>%
  ggplot() +geom_histogram(aes(x=log(numberPositive+1)), bins=100)
tck_borrelia %>%
  filter(numberPositive>0) %>%
  ggplot() +geom_histogram(aes(x=(numberPositive)), bins=100)
tck_borrelia %>%
  filter(numberPositive>0) %>%
  ggplot() +geom_histogram(aes(x=log(numberPositive)))

# What is relationship between proportion positive and number positive
tck_borrelia %>%
  ggplot() +geom_point(aes(x=log(numberPositive), y=log(proportionPositive)))
# What in the world...? Is that straight line? It exists across year, domainID, month, nlcdClass....
tck_borrelia %>%
  filter(nlcdClass=="deciduousForest") %>% # Switch out for year, domain, month, etc
  ggplot() +geom_point(aes(x=log(numberPositive), y=log(proportionPositive), cex=log(nNymph)))

#### Using basic poisson model #####
mod.pois1 <- glm(numberPositive ~ factor(month) + nNymph + nlcdClass
                , family = poisson (link="log")
                ,  data=tck_borrelia
                )
summary(mod.pois1)
mu <- exp(predict(mod.pois1, type="response")) # estimated mean
exp <- sum(dpois(x=0, lambda=mu)) # get probability of zero, then add those up to get total zeros expected?
round(exp)
sum(tck_borrelia$numberPositive==0)
# There's definitely underfitting of zeros.

#### Basic hurdle model #####
# Now, let's try a simple hurdle model with just month and nNymph.
library("pscl") # for hurdle
mod.hurdle <- hurdle(numberPositive ~ 1
       , data=tck_borrelia
       # , offset=numberTested
       , dist = "poisson"
       )
summary(mod.hurdle)
predict(mod.hurdle)
# Why does it predict a constant..? Beause no predictors?

mod.hurdle1 <- hurdle(numberPositive ~ nNymph
                     , data=tck_borrelia
                     # , offset=numberTested
                     , dist = "poisson"
)
summary(mod.hurdle1)
plot(tck_borrelia$numberPositive, predict(mod.hurdle1))
plot((tck_borrelia$numberPositive-predict(mod.hurdle1))~ predict(mod.hurdle1))



mod.hurdle2 <- hurdle(numberPositive ~ nNymph + nlcdClass
                      , data=tck_borrelia
                      # , offset=numberTested
                      , dist = "poisson"
)
summary(mod.hurdle2)
plot(tck_borrelia$numberPositive, predict(mod.hurdle2))
plot((tck_borrelia$numberPositive-predict(mod.hurdle2))~ predict(mod.hurdle2))

# The models keep failing when I try to use "complicated" models. Not sure what it is-- not enough data?

#### GAM ####

### First thing: when alone, does the number of nymphs change the number of positive
mod.gam <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) 
    , data=tck_borrelia
    , method="REML")
summary(mod.gam)
plot(mod.gam, pages=1)

mod.gam1 <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + nNymph
               , data=tck_borrelia
               , method="REML"
)
summary(mod.gam1)
plot(mod.gam1, pages=1)

mod.gam1 <- gam(proportionPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + nNymph
                , data=tck_borrelia
                , method="REML"
)
summary(mod.gam1)
plot(mod.gam1, pages=1)

tck_borrelia %>% ggplot() + geom_point(aes(x=log(numberTested), y=log(nNymph)))

mod.gam2 <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + nNymph 
                , data=tck_borrelia
                , method="REML")
summary(mod.gam2)
plot(mod.gam2, pages=1)
plot(mod.gam2$residuals ~ mod.gam2$fitted.values)
plot(mod.gam2$fitted.values ~ tck_borrelia$numberPositive, ylim=c(-200,200))
# There's one outlier

mod.gam3 <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + nNymph + s(plotID, bs="re")
                , data=tck_borrelia
                , method="REML")
summary(mod.gam3)
plot(mod.gam3, pages=1)
plot(mod.gam3$residuals ~ mod.gam3$fitted.values)
plot(mod.gam3$fitted.values ~ tck_borrelia$numberPositive, ylim=c(-200,200))
# There's one outlier

mod.gam3a <- gam(numberPositive ~ offset(numberTested) + s(dayOfYear, sp=2) + s(plotID, bs="re")
                , data=tck_borrelia
                , method="REML")
summary(mod.gam3a)
plot(mod.gam3a, pages=1)
plot(mod.gam3a$residuals ~ mod.gam3a$fitted.values)
plot(mod.gam3a$fitted.values ~ tck_borrelia$numberPositive, ylim=c(-200,200))
# There's one outlier

mod.gam4 <- gam(log(numberPositive+1) ~ offset(log(numberTested)) + s(dayOfYear, sp=2) + nNymph + s(plotID, bs="re")+nlcdClass
                , data=tck_borrelia
                , method="REML")
mod.gam4 <- gam(log(numberPositive+1) ~ offset(log(numberTested)) + s(dayOfYear, sp=6) + nNymph 
                , data=tck_borrelia
                , method="REML")
summary(mod.gam4)
plot(mod.gam4, pages=1)
gam.check(mod.gam4)


plot(mod.gam4$residuals ~ mod.gam4$fitted.values)
plot(mod.gam4$fitted.values ~ log(tck_borrelia$numberPositive+1)) 
abline(a = 0, b=1)
# There's one outlier

#### Binomial ######
### Hmmmm not working. Let's try a poor-man's hurdle model
bin.intercept <- glm(borrPresent ~ 1
    , data=tck_borrelia
    , family="binomial")
inv.logit(bin.intercept$coefficients)

bin.mod1 <- glm(borrPresent ~ -1 +factor(month) + nlcdClass + nNymph
    , data=tck_borrelia
    # , offset=numberTested
    , family="binomial"
    )
summary(bin.mod1)
Anova(bin.mod1, type = 3)


# Now, take the non-zeros and predict
tck_borrelia_nozeros <- tck_borrelia %>%
  mutate(month=factor(month)) %>%
  filter(numberPositive>0)
# qpois.mod2 <- glm(numberPositive ~ -1 + factor(month) + nlcdClass + nNymph
#                 , data=tck_borrelia_nozeros
#                 , family="quasipoisson"
#                 , offset=numberTested
#                 ) #### Doesn't work-- dispersion is like 1e14!!
qpois.mod2 <- glm(numberPositive ~ -1 + factor(month) + nlcdClass + nNymph
                  , data=tck_borrelia_nozeros
                  , family="quasipoisson"
                  # , offset=numberTested
)
summary(qpois.mod2)
exp(coef(qpois.mod2))
Anova(qpois.mod2, type=3)


tck_borrelia_nozeros %>%
  ggplot() +geom_histogram(aes(x=log(numberPositive)))
tck_borrelia_nozeros %>%
  ggplot() +geom_histogram(aes(x=proportionPositive))

# See if the distribution looks right

cbind(tck_borrelia_nozeros, pred=predict(qpois.mod2)) %>%
  ggplot() + geom_point(aes(x=(numberPositive), y=exp(pred), col=numberTested)) +
  geom_abline(aes(intercept=0, slope=1))
cbind(tck_borrelia_nozeros, pred=predict(qpois.mod2)) %>%
  ggplot() + geom_histogram(aes(x=numberPositive))+geom_histogram(aes(x=exp(pred)), col="red", alpha=0.5)


### Adding random effects ####

ggnorm.mod3 <- lmer(log(numberPositive) ~ -1 + factor(month) + nlcdClass + (1|domainID)
                , data=tck_borrelia_nozeros
                # , family="poisson"
                # , offset=numberTested
)
summary(ggnorm.mod3)
Anova(ggnorm.mod3)
ranef(ggnorm.mod3)


#### Trying GAM ####

gam1 <- gam(numberPositive ~ s(dayOfYear)+ s(nNymph)+ nlcdClass
    ,data=tck_borrelia
    , method = "REML"
    , family = "quasipoisson")
gam.check(gam1)
summary(gam1)
plot(gam1, shift=coef(gam1)[1], pages=1, seWithMean = TRUE, scheme = 2)
# 
# gam2 <- gam(numberPositive ~ s(dayOfYear) + nNymph+ ti(dayOfYear,nNymph)+nlcdClass
#             ,data=tck_borrelia_fullzeroes
#             , method = "REML"
#             )
# gam.check(gam2)
# summary(gam2)
# plot(gam2, shift=coef(gam2)[1], seWithMean = TRUE, scheme=1, pages=1)
# vis.gam(gam2, view=c("dayOfYear","nNymph"), too.far = 0.05, theta=45)

gam3 <- gam(numberPositive ~ s(dayOfYear) + s(nNymph)+ ti(dayOfYear,nNymph) + nlcdClass
            ,data=tck_borrelia_fullzeroes
            , method = "REML"
)
gam.check(gam3)
summary(gam3)
plot(gam3, shift=coef(gam3)[1], seWithMean = TRUE, scheme=1, pages=1)
vis.gam(gam3, view=c("dayOfYear","nNymph"))


gam4 <- gam(proportionPositive ~ ti(dayOfYear,nNymph) + nlcdClass + s(dayOfYear) + nNymph
            ,data=tck_borrelia_fullzeroes
            , method = "REML"
)
summary(gam4)
plot(gam4, shift=coef(gam4)[1], seWithMean = TRUE, scheme=1, pages=1)
vis.gam(gam4, view=c("dayOfYear","nNymph"))

gam5 <- gam(proportionPositive ~ ti(dayOfYear,nNymph) + nlcdClass + s(dayOfYear) + s(nNymph)
            ,data=tck_borrelia_fullzeroes
            , method = "REML"
)
summary(gam5)
plot(gam5, shift=coef(gam5)[1], seWithMean = TRUE, scheme=1, pages=1)
vis.gam(gam5, view=c("dayOfYear","nNymph"))


tck_borrelia_nozeros %>%
  ggplot(aes(x=dayOfYear)) +geom_point(aes(y=numberPositive)) + geom_point(aes(y=proportionPositive*400), col="red") +
  scale_y_continuous(name="Number positive for Borrelia", sec.axis=sec_axis(~./400, name="Proportion Positive")) +
  geom_point(aes(y=log(nNymph)*50), col="blue") +
  xlim(0,365)
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
  