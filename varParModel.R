### Variance partioning model to partition variation in (1) total tick abundance (2) IXODES tick 
### abundance (3) AMBAME tick abundance

### Random effects to look at: site, plot, year, month, nlcd, nlcd nested within site

### 26 Feb 2020
### MEB

### Not using test/validation datasets for this analysis. Using complete datasets of all and each species.


###-----------------------------------------------------------------------------------------------

library(lme4)



### BOTH SPECIES

# read in complete dataset that includes both species
comdat <- read.csv("data_derived/subset_data/ticks_taxonomy_collapsed/complete_data_taxonomy_collapsed.csv")
head(comdat)

# BOTH SPECIES, all stages
vp.c <- glmer(Count ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = comdat)
summary(vp.c)
# variances are the same magnitude except for domain, which is 11 times smaller.

### BOTH SPECIES, adults only
comdat.a <- subset(comdat, lifeStage == "Adult")
vp.c.a <- glmer(Count ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = comdat.a)
summary(vp.c.a)
# month has the most variance but everything besides nlcd class is on the same scale. nlcd is 12 times smaller.

### BOTH SPECIES, nymphs only
comdat.n <- subset(comdat, lifeStage == "Nymph")
vp.c.n <- glmer(Count ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = comdat.n)
summary(vp.c.n)
# domain shows the most variation for nymphs. duplicate siteID.1???

### BOTH SPECIES, larvae only
comdat.l <- subset(comdat, lifeStage == "Larva")
vp.c.l <- glmer(Count ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = comdat.l)
summary(vp.c.l)
# Model fails to converge. Duplicate siteID.1 is also listed again... SiteID.1, month, and nlcd have high variances.


###------------------------------------------------------------------------------------------


### AMBAME

# read in complete dataset that includes both species
ambdat <- read.csv("data_derived/subset_data/ticks_AMBAME/complete_data_AMBAME.csv")
head(ambdat)

# all stages
vp.a <- glmer(rawCount ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = ambdat)
summary(vp.a)
# model fails to converge. domain is way higher than everything else.

# adults only
ambdat.a <- subset(ambdat, lifeStage == "Adult")
vp.a.a <- glmer(rawCount ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = ambdat.a)
summary(vp.a.a)
# domain and month have the highest variance

# nymphs only
ambdat.n <- subset(ambdat, lifeStage == "Nymph")
vp.a.n <- glmer(rawCount ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = ambdat.n)
summary(vp.a.n)
# domain is way higher than everything else 

# larvae only
ambdat.l <- subset(ambdat, lifeStage == "Larva")
vp.a.l <- glmer(rawCount ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = ambdat.l)
summary(vp.a.l)
# most of variance in domain and year. month and ncld are 10 times smaller.



###------------------------------------------------------------------------------------------


### IXOSCA

# read in complete dataset that includes both species
ixodat <- read.csv("data_derived/subset_data/ticks_IXOSCA/complete_data_IXOSCA.csv")
head(ixodat)

# all stages
vp.i <- glmer(rawCount ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = ixodat)
summary(vp.i)
# domain is 5 times higher than everything else.

# adults only
ixodat.a <- subset(ixodat, lifeStage == "Adult")
vp.i.a <- glmer(rawCount ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = ixodat.a)
summary(vp.i.a)
# domain and month have the highest variance

# nymphs only
ixodat.n <- subset(ixodat, lifeStage == "Nymph")
vp.i.n <- glmer(rawCount ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = ixodat.n)
summary(vp.i.n)
# model failed to converge. domain is way higher than everything else 

# larvae only
ixodat.l <- subset(ixodat, lifeStage == "Larva")
vp.i.l <- glmer(rawCount ~ 1 + (1|domainID) + (1|siteID) + (1|plotID) + (1|nlcdClass) + (1|siteID/nlcdClass) + (1|year) + (1|month), family = poisson, data = ixodat.l)
summary(vp.i.l)
# looks weird. ncld nested within site and year are 192 and 174, respectively. everything else is 
# zero (or close to it).





### CONCLUSIONS
# AMBAME and IXOSCA show similar variances across random effects. Most is in domain and month.


### ISSUES
# siteID and siteID.1 are both listed. Is there a duplicate site or something causing this?
# Is poisson dist correct? Need to do zero inflated/neg binom?
# Some models fail to converge.
# Should probably plot these.
