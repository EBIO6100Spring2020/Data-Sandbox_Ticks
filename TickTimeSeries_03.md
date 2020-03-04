About
-----

The goal of this analysis is to use time series analysis to predict tick densities at NEON plots. Previously, data were cleaned by M. Chen, see: `TickMerged_02_DataFiltering`. I am working with data generated from this script.

In this script, I will use a random walk state-space model derived from Dietze's ecological forecasting book (see exercise 6).

Load libraries and Data
-----------------------

Currently, just working with the training dataset.

JAGS is used to fit the random walk models; this code may give you trouble if you do not have JAGS installed. See: <https://www.r-bloggers.com/getting-started-with-jags-rjags-and-bayesian-modelling/>

``` r
library(dplyr)
library(ggplot2)
library(tidyr)
library(padr)
library(rjags)
# read in training set
tick_abun <- read.csv("data_derived/subset_data/ticks_IXOSCA/time_testset_IXOSCA.csv")
```

Filter and clean data
---------------------

Currently all plots are in the dataset, but for now we only want to work with plots that have never had tick nymphs.

``` r
tick_abun %>% group_by(plotID) %>% filter(lifeStage == "Nymph") %>% mutate(totalTicks = sum(estimatedCount, na.rm=TRUE)) %>% filter(totalTicks > 0) %>% ungroup() %>% droplevels() -> nymph
```

We also need to get a measure of tick density which is the raw counts standardized by the tick drag length. In the future a Poisson with offset could probably be used instead of directly estimating density, but for simplicity we will do it first.

``` r
nymph <-nymph %>% mutate(density = estimatedCount/totalSampledArea) %>% filter(!is.na(density))
# fix date
nymph$collectDate <- as.Date(nymph$collectDate)
```

We do not have regular observations; ticks are observed roughly every 3 weeks in summer, but not in fall and winter. Some observations are also missing.

There are probably more elegant ways to deal with this by using continuous time models, but right now I only know how to do discrete time steps, so we need to fill in the data with the missing observations.

I will have each time step be a week, and fill in the missing weeks with NAs. I imagine this time period matters -- especially when we look at the time steps without observations (winter) and how much the uncertainty changes. Something to deal with later...

The package `padr` is amazing, and will fill in missing weeks with NAs. The way this code works is: (1) bin observations into week-long bins (2) if there are &gt;1 observations within that week, take average density (3) fill in missing weeks with NAs.

Each plot will be a separate list item.

``` r
nymph_plots <- list()
for(i in 1:length(levels(nymph$plotID))){
  nymph %>% filter(plotID == levels(plotID)[i]) %>% select(plotID, density, collectDate)%>%
    thicken("week") %>%
    group_by(plotID, collectDate_week) %>%
    summarise(week_density = mean(density, na.rm = TRUE)) %>%
    pad() -> nymph_plots[[i]]}
names(nymph_plots) <- levels(nymph$plotID)
```

Make sure it worked with a few different plots:

``` r
nymph_plots[[1]] %>% filter(!is.na(week_density)) %>%
  ggplot(aes(y=week_density, x = collectDate_week, group = 1))+
  geom_line()
```

![](TickTimeSeries_03_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
nymph_plots[[2]] %>% filter(!is.na(week_density)) %>%
  ggplot(aes(y=week_density, x = collectDate_week, group = 1))+
  geom_line()
```

![](TickTimeSeries_03_files/figure-markdown_github/unnamed-chunk-5-2.png)

``` r
nymph_plots[[3]] %>% filter(!is.na(week_density)) %>%
  ggplot(aes(y=week_density, x = collectDate_week, group = 1))+
  geom_line()
```

![](TickTimeSeries_03_files/figure-markdown_github/unnamed-chunk-5-3.png)

Pick the plot with the most data to run the practice model:

``` r
# longest time series
which(sapply(nymph_plots, nrow) == max(sapply(nymph_plots, nrow)))
```

    ## OSBS_001 
    ##       23

``` r
# most observations
nymph %>% group_by(plotID) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup() %>% slice(1)
```

    ## # A tibble: 1 x 2
    ##   plotID       n
    ##   <fct>    <int>
    ## 1 SCBI_002    44

Run state space model
---------------------

See Dietze's for details on how this model works but briefly, there are two components to the state space model: 1) a state process: this is the "true" unobserved tick density at a given time point. Also called a latent variable. This changes over time, with each time step a random deviation from the previous time step. 2) an observation process: this represents the observed data, and encompasses observation error. At some time points we observe tick density with error that is a random variable.

The "random walk"" part of this model describes how at each time step there is a change in density that is randomly drawn from a distribution.

First write the jags model:

``` r
randomWalk = "model{
  ### observations
  for(t in 1:n){
    y[t] ~ dnorm(x[t], tau_obs) # obs tick density at time t comes from a normal centered around true state, with obs error tao ()
  }
  
  ### process 
  for(t in 2:n){ # starting in the second time step
    x[t] ~ dnorm(x[t-1], tau_add) # true density at time t comes from a distribution centered on prior density, with random fluctuation
  }
  
  ### priors
  x[1] ~ dnorm(x_ic,tau_ic) # state in time 1 is randomly drawn from a distribution centered on mean
  tau_obs ~ dgamma(a_obs,r_obs) # observation error drawn from a gamma distribution (always positive)
  tau_add ~ dgamma(a_add,r_add) # time fluctuation drawn from a comes from gamma (always pos)
}
"
```

Tell jags what data should be used to fit the model:

1.  the actual observation data, y. We will use the log density to somewhat normalize it, but note that this is not a great idea and we will need to modify in the future...

2.  the number of time steps, n

3.  x\_ic which is the mean of the prior distribtuion for the density in time step 0 (we will use the global mean)

4.  tau\_ic which is the precision of the prior distribution for the density in time step 0 (Dietze uses sd which is kind of weird)

5.  a\_obs, the shape parameter for the prior for observation error (hyperparameter?)

6.  r\_obs, the rate parameter for the prior for observation error

7.  a\_add, the shape parameter for the prior for process error

8.  a\_obs, the rate parameter for the prior for process error

In this model, we are trying to estimate process and observation error (taus).

``` r
y = log(nymph_plots$SCBI_002$week_density+1)
data <- list(y = y, n = length(y), x_ic = mean(y, na.rm = TRUE),
             tau_ic = sd(y, na.rm=TRUE),
             a_obs = 1, r_obs = 1, a_add = 1, r_add = 1)
```

Set arguments for jags run:

1.  number of chains

2.  initial values for chains

``` r
nchain = 3
init <- list()
for(i in 1:nchain){ # for each chain initialize with diff random numbers
  y.samp = sample(y,length(y),replace=TRUE)
  # time diffs each random y sample take difference btw each value
  # overall obs error is taken from overall variation
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),tau_obs=5/var(log(y.samp)))
}
```

Fit jags model

``` r
# initialize
j.model   <- jags.model (file = textConnection(randomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 44
    ##    Unobserved stochastic nodes: 410
    ##    Total graph size: 461
    ## 
    ## Initializing model

``` r
# start model
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
# now get posterior distribution
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)
```

Visualize fit of state space model
----------------------------------

First summarise posterior distribution with confidence intervals:

``` r
out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
nymph_plots$SCBI_002$pred <- ci[2,]
nymph_plots$SCBI_002$pred_lo <- ci[1,]
nymph_plots$SCBI_002$pred_hi <- ci[3,]
```

Plot data on the raw scale:

``` r
ggplot(data = nymph_plots$SCBI_002, aes(x=collectDate_week)) +
  geom_line(aes(y = pred, group = 1))+
  geom_point(aes(y = week_density), alpha = .8, col = "red")+
  theme_classic()+
  scale_y_continuous(limits = c(0, max(nymph_plots$SCBI_002$week_density, na.rm = TRUE)))+
  ylab("Tick Density")+
  xlab("Date")
```

    ## Warning: Removed 182 rows containing missing values (geom_point).

![](TickTimeSeries_03_files/figure-markdown_github/unnamed-chunk-12-1.png)

Try on log scale:

``` r
nymph_plots$SCBI_002$predlog <- apply(out[,x.cols],2,mean) ## model was fit on log scale
ggplot(data = nymph_plots$SCBI_002, aes(x=collectDate_week)) +
  geom_line(aes(y = predlog, group = 1))+
  geom_point(aes(x=collectDate_week, y = log(week_density+1)), alpha = .7, col = "red")+
  theme_classic()+
  theme(plot.title = element_text(hjust = .5))+
  scale_y_continuous(limits = c(-0.1, log(max(nymph_plots$SCBI_002$week_density, na.rm = TRUE))))+
  ggtitle("SCBI")+
  ylab("Tick Density")+
  xlab("Date")
```

    ## Warning: Removed 183 rows containing missing values (geom_point).

![](TickTimeSeries_03_files/figure-markdown_github/unnamed-chunk-13-1.png)

Try with another site:
----------------------

``` r
y = log(nymph_plots$SCBI_013$week_density+1)
data <- list(y = y, n = length(y), x_ic = mean(y, na.rm = TRUE),
             tau_ic = sd(y, na.rm=TRUE),
             a_obs = 1, r_obs = 1, a_add = 1, r_add = 1)
init <- list()
for(i in 1:nchain){ # for each chain initialize with diff random numbers
  y.samp = sample(y,length(y),replace=TRUE)
  # time diffs each random y sample take difference btw each value
  # overall obs error is taken from overall variation
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),tau_obs=5/var(log(y.samp)))
}
# initialize
j.model   <- jags.model (file = textConnection(randomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 41
    ##    Unobserved stochastic nodes: 413
    ##    Total graph size: 461
    ## 
    ## Initializing model

``` r
# start model
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
# now get posterior distribution
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)
```

``` r
out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
nymph_plots$SCBI_013$pred <- ci[2,]
nymph_plots$SCBI_013$pred_lo <- ci[1,]
nymph_plots$SCBI_013$pred_hi <- ci[3,]

ggplot(data = nymph_plots$SCBI_013, aes(x=collectDate_week)) +
  geom_line(aes(y = pred, group = 1))+
  geom_point(aes(y = week_density), alpha = .8, col = "red")+
  theme_classic()+
  scale_y_continuous(limits = c(0, max(nymph_plots$SCBI_013$week_density, na.rm = TRUE)))+
  ylab("Tick Density")+
  xlab("Date")
```

    ## Warning: Removed 185 rows containing missing values (geom_point).

![](TickTimeSeries_03_files/figure-markdown_github/unnamed-chunk-15-1.png)

Next steps/Questions
--------------------

-   Log transforming doesn't make tick density anywhere near normal, and the +1 creates some discrepencies because of all the 0s. Should at least use a link function

-   Need to loop through all the sites. Could do with one jags model, just need to index y list elements.

-   Right now the prior for log-density in year 0 is normally distributed. We probably want to draw it from some kind of overdispersed poisson.

-   The observation errors and process errors are also normal...does this make sense? Are some deviations more likely than others given the distribution of tick densities?

-   Lack of observation in the fall/winter means that the state variables start to diverge a lot. If we plot uncertainty, it balloons during the periods without observations. This makes sense mathematically, but not biologically... we need to somehow constrain tick densities. Would adding a covariate for month help?

-   Time is still discrete. How to fix this?

-   How to build in other covariates...?

-   Should this really be a multivariate time series with multiple plots? instead of a separate time series for each plot? Plots with lots of data should inform other time series, maybe...

-   If so, we need to model the hierarchy: plots nested within sites nested within domains.

-   Still need to forecast and evaluate accuracy. Is this as simple as going from 1:n to 1:n+1 for the state variable x?

-   What's the deal with Kalman filters and smoothing?