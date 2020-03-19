######## using a hurdle model ##########
library(tidyverse)
library(pscl) # hurdle models
library(lme4) # mixed effect models
library(boot) # logit/inv logit
library(car) # Anova (type III)

tck <- read.csv("data_derived/MASTER_all_tck_data_merged.csv")

#### Filter tick dataset and adjust ####
tck %>%
  mutate(tested=ifelse(is.na(testingID),"no","yes")) %>%
  select(plotID, tested) %>% table(useNA="ifany")
# There are only a few plots that are actually tested for any tick pathogen.
# Let's filter those out.
testedPlots <- tck %>%
  mutate(tested=ifelse(is.na(testingID),0,1)) %>%
  group_by(plotID) %>%
  summarize(anyTested=sum(tested)) %>% filter(anyTested>0) %>%pull(plotID) %>%as.character()

# Okay, let's also filter by plots that NEVER had any positive borrelia
infectedPlots <- tck %>%
  filter(Borrelia_sp.=="Positive") %>%
  select(plotID) %>%pull() %>% unique()

tck_borrelia <- tck %>%
  filter(plotID %in% as.character(infectedPlots)) %>%
  mutate(tested=ifelse(is.na(Borrelia_sp.),0,1)
         ,isPositive=ifelse(Borrelia_sp.=="Positive", 1,0)
  ) %>%
  group_by(domainID, siteID, nlcdClass,plotID, elevation, collectDate, dayOfYear, year, month) %>%
  summarize(numberTested=sum(tested)
            ,n=n()
            , numberPositive=sum(isPositive,na.rm = TRUE)
            , nAdult=sum(lifeStage=="Adult")
            , nNymph=sum(lifeStage=="Nymph")
            , nLarva=sum(lifeStage=="Larva")) %>%
  ungroup() %>%
  mutate(proportionTested=numberTested/n
         , proportionPositive=numberPositive/numberTested) %>%
  mutate(proportionTested=ifelse(proportionTested==0,NA,proportionTested)
         , proportionPositive=ifelse(proportionPositive==0,NA,proportionPositive)
         , nlcdClass=factor(nlcdClass, levels=c("emergentHerbaceousWetlands","cultivatedCrops","pastureHay","grasslandHerbaceous"
                                                ,"dwarfScrub","shrubScrub","sedgeHerbaceous"
                                                ,"woodyWetlands","deciduousForest","evergreenForest","mixedForest"))) %>%
  mutate(borrPresent = ifelse(numberPositive>0,1,0)
         , domainID = factor(as.character(domainID))
         , siteID = factor(as.character(siteID))
         , plotID = factor(as.character(plotID)))


#### Preliminary Plotting ####

# First, what is the distribution of Borrelia prevalence?
tck_borrelia %>%
  ggplot() +geom_histogram(aes(x=log(numberPositive+1)))
tck_borrelia %>%
  filter(numberPositive>0) %>%
  ggplot() +geom_histogram(aes(x=(numberPositive)))
tck_borrelia %>%
  filter(numberPositive>0) %>%
  ggplot() +geom_histogram(aes(x=log(numberPositive+1)))


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
mod.hurdle <- hurdle(numberPositive ~ factor(month)
       , data=tck_borrelia
       # , offset=numberTested
       )

mod.hurdle <- hurdle(numberPositive ~ 1
                     , data=tck_borrelia
                     # , offset = numberTested
                     # , dist = "negbin"
                     )

summary(mod.hurdle)
inv.logit(-2.3345)
#### Binomial ######
### Hmmmm not working. Let's try a poor-man's hurdle model
bin.intercept <- glm(borrPresent ~ 1
    , data=tck_borrelia
    , family="binomial")
inv.logit(bin.intercept$coefficients)

bin.mod1 <- glm(borrPresent ~ -1 +factor(month) + nlcdClass + nNymph
    , data=tck_borrelia
    , family="binomial"
    )
summary(bin.mod1)
Anova(bin.mod1, type = 3)


# Now, take the non-zeros and predict
tck_borrelia_nozeros <- tck_borrelia %>%
  mutate(month=factor(month)) %>%
  filter(numberPositive>0)
qpois.mod2 <- glm(numberPositive ~ -1 + factor(month) + nlcdClass + nNymph
                , data=tck_borrelia_nozeros
                , family="quasipoisson"
                # , offset=numberTested
                )
qpois.mod2 <- glm(numberPositive ~ 1
                  , data=tck_borrelia_nozeros
                  , family="quasipoisson"
                  # , offset=numberTested
)
summary(qpois.mod2)
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






# Let's try to simulate real data?
rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

nsim <- 1000
bin_coeff <- bin.mod1$coefficients
bin_pred <- predict(bin.mod1, newdata=tck_borrelia_nozeros,type="response")
pois_coeff <- c(0,bin.mod2$coefficients, 0, 0, 0)

var_overdisp <- 37.73
# x.pred <- tck_borrelia %>%
#   select(month)
x.pred <- predict(bin.mod2, type="response")
pred.mat <- matrix(ncol=length(x.pred), nrow=nsim, dimnames = list(seq(1:nsim), x.pred))
for ( n in 1:nsim ) {
  # Binomial
  # temp.bin <- rbinom(n=length(x.pred), size=1, prob= exp(bin_coeff))
  # temp.X <- pois_coeff[as.character(x.pred)]
  # temp.pois <- rnorm(n=length(temp.X), mean=temp.X, sd=temp.X*sqrt(var_overdisp))
  # temp.Y <- temp.bin*exp(temp.pois)
  # pred.mat[n,] <- temp.Y
  # Try neg binomial?
  pred.mat[n,] <- rqpois(n=length(x.pred), mu=x.pred, theta=var_overdisp)
}

pred.mat.long <- pred.mat %>% 
  as_tibble() %>%
  gather(key="Predicted", value="Simulated") %>%
  mutate(Predicted=as.numeric(Predicted), Observed= rep(tck_borrelia_nozeros$numberPositive, each=nsim))

pred.mat.long %>%
  ggplot() + geom_point(aes(x=log(Observed), y=log(Simulated)), alpha=0.5)


cbind(tck_borrelia_nozeros, pred = inv.logit(predict(bin.mod2))) %>%
  ggplot() +
  geom_line(aes(x=month, y=predPropPos)) +
  geom_point(aes(x=month, y=propPositive))
  