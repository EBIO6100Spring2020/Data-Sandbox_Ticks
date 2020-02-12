### ABOUT
## Script to merge NEON Tick abundance and NEON Parasite prevalence data
## 30 Jan 2020
## WEM
library(dplyr)
library(tidyverse)
library(forcats)
library(vegan)
library(lubridate)
library(ggthemes)  # for a mapping theme
library(ggalt)  # for custom map projections
library(viridis)
library(ggrepel)  # for annotations

#### PLOT THEME #############
theme_tick <- function(){
  theme_bw() +
    theme(text = element_text(family = "Arial"),
          axis.text = element_text(size = 11), 
          axis.title = element_text(size = 14),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 16, vjust = 1, hjust = 0),
          legend.text = element_text(size = 9),          
          legend.title = element_blank(),                              
          legend.position = "bottom", 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

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
ggplot(aes(y = log(totalAbun+1), x = Month,  group = domainID), data = TickAbunTax)+
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


### get some summary statistics ############
TickAbunPrev %>% group_by(domainID) %>% summarise(
  max.abun = max(totalAbun),
  mean.abun = mean(totalAbun, na.rm = TRUE),
  mean.Borrelia.prev = mean(Borrelia.sp., na.rm=TRUE),
  min.year = min(year(as.Date(collectDate))),
  max.year = max(year(as.Date(collectDate))),
  n.years = length(unique(year(as.Date(collectDate)))),
  latitude = first(decimalLatitude),
  longitude = first(decimalLongitude)) %>% filter(mean.abun > 0) -> Domain.Aves

### Plot ####
north_america <- map_data("world", region = c("USA", "Canada"))
north_america <- north_america[!(north_america$subregion %in% "Hawaii"),]
(tick_map1 <- ggplot() +
    geom_map(map = north_america, data = north_america,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj(paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96",
                      " +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"),  ylim = c(25, 70), xlim = c(-150, -60)) +
  
    # Add points for the site locations
    geom_point(data = Domain.Aves, 
               aes(x = longitude, y = latitude, fill = mean.abun),
               alpha = 0.8, size = 4, colour = "grey30",
               shape = 21)+
  scale_fill_viridis(option = "magma", direction = -1, begin = 0.2)+
  geom_label_repel(data = Domain.Aves,
                   aes(x = longitude, y = latitude,
                       label = domainID),
                   # Setting the positions of the labels
                   box.padding = .5, size = 3, nudge_x = .5, nudge_y = .5)+
  theme_tick()+
  theme(plot.title = element_text(hjust =0.5), axis.ticks = element_blank(), axis.text = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Mean Tick Abundance"))
ggsave(tick_map1, filename = "output/figs/TickAbundance.png", height = 6, width = 6)

(tick_map2 <- ggplot() +
    geom_map(map = north_america, data = north_america,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj(paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96",
                      " +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"),  
               ylim = c(25, 70), xlim = c(-150, -60)) +
    # Add points for the site locations
    geom_point(data = Domain.Aves, 
               aes(x = longitude, y = latitude, fill = mean.Borrelia.prev),
               alpha = 0.8, size = 4, colour = "grey30",
               shape = 21)+
    scale_fill_viridis(option = "magma", direction = -1, begin = 0.2)+
  geom_label_repel(data = Domain.Aves,aes(x = longitude, y = latitude,label = domainID),
                   box.padding = .5, size = 3, nudge_x = .5, nudge_y = .5)+
  theme_tick()+
  theme(plot.title = element_text(hjust =0.5), axis.ticks = element_blank(), axis.text = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Mean Borrelia Prevalence"))
ggsave(tick_map2, filename = "output/figs/Borrelia_prev.png",
       height = 6, width = 6)
