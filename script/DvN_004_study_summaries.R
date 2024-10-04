#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "05/06/2024"

###script: 004 - Summaries

#additional libaries
library(countrycode)

###load from DvN_003_Analysis_dataframe_setup.R
load(file="output/DvN_analysis_dataframe Sep16.RData")
diel.all.diffs
load(file="output/DvN_plant_phylogeny Sep16.RData")
diel.tree.out


####Data summaries

#summarise number of effect (rows) per treatment
effect.tots=diel.all.diffs %>% 
  group_by(treatment) %>% 
  summarise(n=n())

study.tots=diel.all.diffs %>% 
  distinct(study_ID) %>% 
  summarise(n=n())

species.tots=diel.all.diffs %>% 
  distinct(phylo) %>% 
  summarise(n=n())

species.tots.trt=diel.all.diffs %>% 
  group_by(treatment) %>% 
  distinct(phylo) %>% 
  summarise(n=n())

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
  #if na change to Africa
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
poll.dep.summaries=es.list$dvn_effects %>% filter(treatment%in%"bag_open") %>% 
  filter(accepted_name%in%diel.all.diffs$accepted_name)

#species
unique(poll.dep.es$accepted_name) %>% length
#studies
unique(poll.dep.es$study_ID) %>% length
#effect sizes
dim(poll.dep.es)
