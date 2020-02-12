### ABOUT
## Script to merge NEON Tick abundance and NEON Parasite prevalence data
## 30 Jan 2020
## WEM
library(dplyr)
library(tidyverse)
library(forcats)
library(vegan)
#### READ IN DATA ##########
# borrelia prevalence at the plot-date-taxonomy levels
# prevalence = # ticks infected with Borrelia sp. / # ticks caught
# prevalence is for all ticks tested of a given species from a given sample
TickPathAgg <- read.csv("data_derived/Tick_Borrelia_Prev_Aggregated.csv")


# tick individual infection data
# each row is a tick with a unique testingID
# separate column for each pathogen that individual was tested for
# NA = not tested
TickPathIndiv <- read.csv("data_derived/Tick_Pathogen_Individual.csv")


# tick abundance data from tick drags
# each row is from a given plot-date (1 drag per site-plot-date)
# each row has a count of ticks collected
# if > 0, tick sample ID was assigned
TickAbun <- read.csv("data_derived/Tick_Abundance_Counts.csv")
TickAbun$PlotID_Date <- paste(TickAbun$plotID, TickAbun$collectDate, sep = "_")
length(unique(as.factor(TickAbun$PlotID_Date))) # 7644

# tick taxonomy data from tick drags
# each sample ID is grouped into subsamples of unique taxonomy-age-sex
# to get the total counts from each sample we need to reformat data
TickTax <- read.csv("data_derived/Tick_Abundance_Taxonomy.csv")

# right now sex and age info is grouped into same column. clean this up:

TickTax %>% mutate(Age = case_when(sexOrAge == "Adult" ~ "Adult", sexOrAge == "Female" ~ "Adult",
                                   sexOrAge=="Larva"~"Larva", sexOrAge == "Male" ~ "Adult", 
                                   sexOrAge == "Nymph" ~ "Nymph"))%>%
          mutate(Age = replace_na(Age, "Unknown")) %>%
          mutate(Sex = case_when(sexOrAge == "Male" ~ "Male", sexOrAge == "Female" ~ "Female")) %>%
          mutate(Sex = replace_na(Sex, "Unknown")) %>%
          mutate(acceptedTaxonID = fct_recode(acceptedTaxonID, "UNKNOWN" = "")) %>%
          mutate(Tax_Age_Sex = paste(acceptedTaxonID, Age, Sex, sep = "_")) -> TickTax
levels(TickTax$acceptedTaxonID)[which(levels(TickTax$acceptedTaxonID)=="")] <- "UNKNOWN"
# widen first
TickTax %>% pivot_wider(id_cols = subsampleID, names_from = Tax_Age_Sex, values_from = individualCount, values_fill = list(individualCount = 0)) %>% right_join(TickTax) -> TickTax_wide
TickTax_wide %>% pivot_wider(id_cols = subsampleID, names_from = acceptedTaxonID, values_from=individualCount, values_fill = list(individualCount = 0)) %>% right_join(TickTax_wide) %>% select(subsampleID, namedLocation, plotID, collectDate, sampleID, individualCount, acceptedTaxonID, everything(.)) -> TickTax_wide

colSums(TickTax_wide[, 27:69], na.rm = TRUE)
TickTax_wide %>% group_by(sampleID) %>% summarise_at(vars(c(DERVAR:IXOANG, DERVAR_Adult_Female:IXOPAC_Nymph_Unknown)), sum, na.rm = TRUE)  -> TickTax_SampleID

##### JOIN TICK TAXONOMY AND COUNTS #####
TickAbunTax <- left_join(TickAbun, TickTax_SampleID, by = "sampleID") 
sp.cols <- which(colnames(TickAbunTax) == "DERVAR"):which(colnames(TickAbunTax)=="IXOANG")
TickAbunTax[, sp.cols][is.na(TickAbunTax[, sp.cols])] <- 0
TickAbunTax$TaxonRichness <- rowSums(decostand(TickAbunTax[, sp.cols], "pa", na.rm = TRUE), na.rm = TRUE)
TickAbunTax[, 21:23][is.na(TickAbunTax[, 21:23])] <- 0
TickAbunTax$totalAbun = TickAbunTax$adultCount + TickAbunTax$nymphCount + TickAbunTax$larvaCount
hist(TickAbunTax$TaxonRichness)
hist(log(TickAbunTax$totalAbun+1))
TickAbunTax$Date <- as.Date(TickAbunTax$collectDate)

TickAbunTax$Month <- month(TickAbunTax$Date)
ggplot(aes(y = totalAbun, x = Month,  group = domainID), data = TickAbunTax)+
        geom_jitter(aes(color = domainID)) 
##### GET TICK DISEASE PREVALENCE BY SITE #####
# fist merge the tick prevalence with tick taxonomy
CalcPrev <- function(x){sum(x, na.rm = TRUE)/sum(!is.na(x))}
TickPrev <- TickPathIndiv %>% group_by(subsampleID) %>% summarise_at(vars(Borrelia.burgdorferi:Borrelia.lonestari),CalcPrev)
TickPrevTax <- left_join(TickTax, TickPrev)
TickPrevTax_sampleID <- TickPrevTax %>% group_by(sampleID) %>% summarise_at(vars(Borrelia.burgdorferi:Borrelia.lonestari), mean, na.rm = TRUE)
TickAbunPrev <- left_join(TickAbunTax, TickPrevTax_sampleID)
colnames(TickAbunPrev)
TickAbunPrev$Risk <- NA
cc <- complete.cases(TickAbunPrev[, c("totalAbun", "Borrelia.sp.")])
TickAbunPrev$Risk[cc] <- TickAbunPrev$totalAbun[cc]*TickAbunPrev$Borrelia.sp.[cc]
ggplot(aes(y = log(Risk+1), x = Month,  group = domainID), data = TickAbunPrev)+
  geom_jitter(aes(color = domainID)) 
plot(log(Risk+1)~log(totalAbun), TickAbunPrev)
