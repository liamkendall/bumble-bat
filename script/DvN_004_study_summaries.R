#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "20/11/2024"

###script: 004 - Summaries

#libraries (CRAN)
library(plyr)
library(dplyr)
library(tidyr)
library(countrycode)

###load analysis dataframes (DvN_003)
load(file="output/DvN_analysis_dataframe final.RData")
load(file="output/DvN_polldep_analysis_dataframe final.RData")

#####load data (from DvN_002)
load("output/DvN_effect_size_dataframe final.rData")

####Data summaries

#summarise number of effect (rows) per treatment
effect.tots=diel.all.diffs %>% 
  group_by(treatment) %>% 
  summarise(n=n())

#number of studies
study.tots=diel.all.diffs %>% 
  distinct(study_ID) %>% 
  summarise(n=n())

#number of species
species.tots=diel.all.diffs %>% 
  distinct(phylo) %>% 
  summarise(n=n())

#species by treatment
species.tots.trt=diel.all.diffs %>% 
  group_by(treatment) %>% 
  distinct(phylo) %>% 
  summarise(n=n())

#by metric
measure.tots=diel.all.diffs %>% 
  group_by(treatment_effectiveness_metric) %>% 
  summarise(n=n())

#get counts of genera and family
genus.tots=diel.all.diffs %>% 
  distinct(genus) %>% 
  summarise(n=n())

family.tots=diel.all.diffs %>% 
  distinct(family) %>% 
  summarise(n=n())

family.spp.tots=diel.all.diffs %>% 
  distinct(family,phylo) %>% 
  group_by(family) %>% 
  summarise(n=n())

##summarise average daylength and DTR
daylength.tots=es.list$environment %>% 
  distinct(study_ID,Daylength) %>%
  summarise(mean=mean(Daylength),
            sd=sd(Daylength),
            se=sd/sqrt(n()))

DTR.tots=es.list$environment %>%
  distinct(study_ID,DTR) %>%
  summarise(mean=mean(DTR),
            sd=sd(DTR),
            se=sd/sqrt(n()))

#elevation
elevation.tots=diel.all.diffs %>% 
  distinct(study_ID,elevation) %>%
  summarise(mean=mean(elevation),
            sd=sd(elevation),
            se=sd/sqrt(n()))

####get continent information
diel.all.diffs$continent=countrycode(diel.all.diffs$country, origin = "country.name", destination = "continent")

#summarise continent by study
continent.tots=diel.all.diffs %>% 
  #if na change to Africa (Canary Islands)
  mutate(continent=ifelse(is.na(continent),"Africa",continent)) %>% #Canary Islands
  distinct(study_ID,continent) %>% 
  group_by(continent) %>% 
  #summarise as percentage
  summarise(perc=n()/135,
            n=n())
 
#summarise number of studies by treatment
treatment.tots=diel.all.diffs %>% 
  distinct(study_ID,treatment) %>% 
  group_by(treatment) %>% 
  summarise(n=n())

#summarise number of studies that studied crops
crop.tots=diel.all.diffs %>% 
  distinct(study_ID,plant_crop) %>% 
  filter(plant_crop%in%"y") %>% 
  summarise(n=n())

#summarise number of plant species that are crops
crop.spp.tots=diel.all.diffs %>% 
  distinct(accepted_name,plant_crop) %>% 
  filter(plant_crop%in%"y") %>% 
  summarise(n=n())

####pollination dependency summaries
poll.dep.summaries=poll.dep.es
#species
unique(poll.dep.es$accepted_name) %>% length
#studies
unique(poll.dep.es$study_ID) %>% length
#effect sizes
dim(poll.dep.es)
