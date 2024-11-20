#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "20/11/2024"

###script: 003 - final set-up for analysis and phylogeny generation

####
#additional functions
source('script/DvN_additional_functions.R')

####libraries (CRAN)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(phytools)
library(metafor)
library(ggplot2)

#non-cran libraries
#to install non-cran libraries
#library("devtools")
#devtools::install_github("jinyizju/V.PhyloMaker2") v.0.1.0
#devtools::install_github("daniel1noble/orchaRd") v.2.0

#load non-cran libraries
library("V.PhyloMaker2")
library(orchaRd) 

#####load data (from DvN_002)
load("output/DvN_effect_size_dataframe final.rData")

#### PREPARE DATA ####

#subset bag vs. open contrasts (i.e. pollination dependency) for joining to main dataframe
diel.dependency=es.list$dvn_effects %>% 
  filter(treatment == "bag_open")%>% 
#  #rename yi and vi
  rename(poll.dep=yi,
         poll.dep.var=vi) %>%
  select(study_ID,poll.dep,poll.dep.var,
         effectiveness_value_complete.exclusion,
         effectiveness_value_open.pollination,
         effectiveness_value_day.pollination,
         effectiveness_value_night.pollination)

#####assemble final dataframe
#remove extreme effect sizes
#scale variables
#add pollination dependency values
#remove one species without photosynthetic pathway

diel.all.diffs=es.list$dvn_effects %>% 
  filter(treatment%in%c("day_night","open_day","open_night"))%>% 
  left_join(diel.dependency)%>%
  #filter yi between -10 and 10
  filter(yi > -10 & yi < 10) %>% #removes five extreme yis
  mutate(genus=word(accepted_name,1)) %>% 
  mutate(sDaylength=as.numeric(scale(Daylength)),
         sDTR=as.numeric(scale(DTR)),
         sElevation=as.numeric(scale(elevation))) %>%  
  mutate(hr.ratio=day_pollination_hrs/night_pollination_hrs) %>%
  mutate(hr.ratio.missing=ifelse(is.na(hr.ratio),"Absent","Present")) %>%
  mutate(hr.ratio.out=ifelse(is.na(hr.ratio),Daylength/(24-Daylength),hr.ratio)) %>%
  mutate(day_pollination_hrs2=ifelse(is.na(day_pollination_hrs),Daylength,day_pollination_hrs)) %>%
  mutate(night_pollination_hrs2=ifelse(is.na(night_pollination_hrs),24-Daylength,night_pollination_hrs)) %>%
  mutate(exposure=ifelse(treatment%in%"open_day",day_pollination_hrs2/24,hr.ratio.out)) %>% 
  mutate(exposure=ifelse(treatment%in%"open_night",night_pollination_hrs2/24,exposure)) %>% 

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

###pollination dependency dataframe
poll.dep.es=es.list$dvn_effects %>% filter(treatment%in%"bag_open") %>% 
  filter(accepted_name%in%diel.all.diffs$accepted_name)%>%
  #filter yi between -10 and 10
  filter(yi > -10 & yi < 10)

###save dataframe and phylogeny
save(diel.all.diffs,file="output/DvN_analysis_dataframe final.RData")
save(poll.dep.es,file="output/DvN_polldep_analysis_dataframe final.RData")
save(diel.tree.out,file="output/DvN_plant_phylogeny final.RData")


