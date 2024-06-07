#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "05/06/2024"

###script: 003 - final set-up for analysis and phylogeny generation

####
#additional functions
source('script/DvN_additional_functions.R')

####Wrangle data to wide format
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

#phylo/taxonomy packages
library(WorldFlora)
library(phytools)
library("V.PhyloMaker2")

#meta-analysis
library(metafor)
library(orchaRd)

#plotting
library(ggplot2)

#####data input
#this dataframe is a product of the set-up scripts
load("output/DvN_effect_size_dataframe June5.rData")

#### PREPARE DATA ####

###calculate mean pollination dependency for each species that had bag open comparison
#then calculate species level
#and genus level

#dependency dataframe
diel.dependency=es.list$dvn_effects %>% 
  filter(treatment == "bag_open") %>% 
  #filter yi between -10 and 10
  filter(yi > -10 & yi < 10) #removes one extreme yi

#metric level dependency
diel.spp.dependency=diel.dependency %>% 
  group_by(accepted_name,treatment_effectiveness_metric) %>%
  summarise(spp.poll.dep=mean(yi),
            spp.poll.dep.var=mean(vi),
            n=n())

#species level
diel.spp.nm.dependency=diel.dependency %>% 
  group_by(accepted_name) %>%
  summarise(spp.nm.pd=mean(yi),
            spp.nm.pd.var=mean(vi),
            n=n())

diel.genus.nm.dependency=diel.dependency %>% 
  mutate(genus=word(accepted_name,1)) %>% 
  group_by(genus) %>%
  summarise(spp.g.pd=mean(yi),
            spp.g.pd.var=mean(vi),
            n=n())

#####assemble final dataframe
#remove extreme effect sizes
#scale variables
#add pollination dependency values
#remove one species without photosynthetic pathway

diel.all.diffs=es.list$dvn_effects %>% 
  filter(treatment%in%c("day_night","open_day","open_night"))%>%
  #filter yi between -10 and 10
  filter(yi > -10 & yi < 10) %>% #removes five extreme yis
  mutate(genus=word(accepted_name,1)) %>% 
  mutate(sDaylength=as.numeric(scale(Daylength)),
         sDTR=as.numeric(scale(DTR)),
         sElevation=as.numeric(scale(elevation))) %>%  
  left_join(diel.spp.dependency,
            by = c("accepted_name", "treatment_effectiveness_metric")) %>%  
  left_join(diel.spp.nm.dependency,
            by = c("accepted_name"))%>%  
  left_join(diel.genus.nm.dependency,
            by = c("genus"))%>%  
  mutate(pd.imp=ifelse(is.na(spp.poll.dep),spp.nm.pd,spp.poll.dep)) %>% 
  mutate(pd.imp.all=ifelse(is.na(pd.imp),spp.g.pd,pd.imp)) %>% 
  ungroup() %>% 
  mutate(sFlower_length=as.numeric(scale((flower_length_midpoint_mm))),
         sStyle_length=as.numeric(scale((style_length_midpoint_mm))),
         sFlower_width=as.numeric(scale((flower_width_midpoint_mm))),
         sPlant_height=as.numeric(scale((plant_height_midpoint_m)))) %>% 
  filter(!is.na(ps_pathway)) 

###changes to traits
#only 6 blue so merged with purple
#only 3 brown so merged with yelloe

diel.all.diffs$color=revalue(diel.all.diffs$color,c("blue"="blue_purple","purple"="blue_purple",
                                                    "brown"="brown_yellow","yellow"="brown_yellow"))

#only 7 capitulum so merged with open
diel.all.diffs$flower_shape=revalue(diel.all.diffs$flower_shape,c("capitulum"="open"))

#BLoom period: if NA, recoded to variable ("both")
diel.all.diffs$bloom_period_simple=ifelse(is.na(diel.all.diffs$bloom_period_simple),
                                          "both",as.character(diel.all.diffs$bloom_period_simple))

###Generate phylogeny with V.phylomaker2
diel.phylo=diel.all.diffs %>% select(phylo,genus,family) %>% 
  distinct(phylo,genus,family) %>% 
  mutate(species=gsub("_"," ",phylo),
         genus=word(gsub("_"," ",phylo),1)) %>% 
  select(species,genus,family) %>% as.data.frame()

# generate a phylogeny for the species list
tree <- phylo.maker(sp.list = diel.phylo, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, scenarios="S3")

#take scenario 3 tree
diel.tree.out=tree$scenario.3

###save dataframe and phylogeny
save(diel.all.diffs,file="output/DvN_analysis_dataframe June5.RData")
save(diel.tree.out,file="output/DvN_plant_phylogeny June5.RData")


