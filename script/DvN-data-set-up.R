####Wrangle data to wide format
library(plyr)
library(dplyr)
library(tidyr)
library(httr)
library(stringr)

#phylo/taxonomy packages
library(WorldFlora)
library(phytools)
library(rtrees)

#meta-analysis
library(metafor)
library(orchaRd)
#plotting
library(ggplot2)

source('script/DvN-functions.R')

#raw dataframe
diel.raw=read.csv("data/DvN_LiteratureDataClean_2023-12-07.csv",#sep=";",dec=".",
                  row.names = 1,
                  encoding = "UTF-8")

#take relevant columns
diel=diel.raw %>% 
  select(study_ID,country,Site_period,coder_initials,
         coordinates_lat,coordinates_lon,
         year_start:night_pollination_hrs,
         plant_species,plant_cultivar,
         treatment_condition,
         additional_treatment,
         treatment_effectiveness_metric,
         sample_size:effect_error_units) %>% 
  filter(!sample_size<=1) %>% 
  filter(!is.na(sample_size)==T) %>% 
  filter(!is.na(effectiveness_value)==T) #%>% 
  #filter(treatment_effectiveness_metric%in%c("fruit set",
  #                                           "seed set"))

####chcek stuff
table(diel$treatment_condition)
table(diel$treatment_effectiveness_metric)
table(diel$treatment_condition,diel$treatment_effectiveness_metric)

#####
sum(is.na(diel$effect_error))

########
#add column of required imputation / transformation

diel=diel %>% mutate(impute.type=case_when(treatment_effectiveness_metric%in%"fruit set" &
                          is.na(effect_error) == T & 
                          is.na(sample_size) == F &
                          effectiveness_value<=1& 
                          !effectiveness_value==0  ~ "binary",
                          
                          is.na(effect_error) == F & 
                            effect_error > 0 &
                            is.na(sample_size) == F & 
                            effect_error_units == "SE" ~  "SED",
                          
                            treatment_effectiveness_metric%in%"fruit set" &
                            effectiveness_value ==0|
                            treatment_effectiveness_metric%in%"fruit set" &
                            effectiveness_value ==1 ~ "zero", 
                          
                            effect_error_units == "CI95" ~
                              "CI95",
                          
                            is.na(effect_error)==T ~ "CV",
                            .default = "none")) 

#34 missings to be imputed later
#12 fruit set
#22 seed set

diel %>% group_by(
  impute.type,treatment_effectiveness_metric) %>% 
  summarise(sum_na = sum(is.na(effect_error)))

######FIX SDs 
#those conformable to binary arrays ("binary")
#from SEs ("SED")
diel.sd=diel %>% 
  rowwise() %>% 
  #fruit sets
  mutate(SDc = ifelse(impute.type%in%"binary", #binary function call
                      binary.summary(effectiveness_value,sample_size)[3], #returns just SD
                      
                      ifelse(impute.type%in%"SED", #else SE to SD transformation
                             effect_error*sqrt(sample_size),
                      effect_error)))%>% 
  #confidence intervals to SD
  mutate(SDc=ifelse(impute.type%in%"CI95",
                    sqrt(sample_size)*((effectiveness_value+effect_error)-       
                                         (effectiveness_value-effect_error))/3.92,
                    SDc)) %>% 
  ###make sure 1s and 0s have 0 variance fruit set
  mutate(SDc=ifelse(impute.type%in%"zero",0,SDc))

###missingness summary
#18 missing from seed set dataframe
#6 missing from pollen deposition
#5 missing from seed mass
diel.sd %>% group_by(impute.type,treatment_effectiveness_metric) %>% 
  summarise(sum_na = sum(is.na(SDc)))

###use Bishop & Nakagawa (2021) method for imputing missing SDs

#When a measure of variance was not obtainable, 
#we imputed standard deviation based on the fitted 
#relationship between mean, SD and n in the available data 
#(91 data rows imputed; adjusted-R2 of these models were between 0.29 and 0.46).

no.missing.all=diel.sd %>% 
  filter(!is.na(SDc))

missing.seed.set=diel.sd%>% 
  filter(treatment_effectiveness_metric%in%"seed set") %>% 
  filter(is.na(SDc))

no.missing.seed.set=diel.sd%>% 
  filter(treatment_effectiveness_metric%in%"seed set") %>% 
  filter(!is.na(SDc))

missing.seed.mass=diel.sd%>% 
  filter(treatment_effectiveness_metric%in%"seed mass") %>% 
  filter(is.na(SDc))


no.missing.seed.mass=diel.sd%>% 
  filter(treatment_effectiveness_metric%in%"seed mass") %>% 
  filter(!is.na(SDc))

missing.pollen.dep=diel.sd%>% 
  filter(treatment_effectiveness_metric%in%"pollen deposition") %>% 
  filter(is.na(SDc))


no.missing.pollen.dep=diel.sd%>% 
  filter(treatment_effectiveness_metric%in%"pollen deposition") %>% 
  filter(!is.na(SDc))


###fit models
seed.set.sd.lm=lm(SDc~effectiveness_value*sample_size*treatment_condition,
         data=no.missing.seed.set)#

summary(seed.set.sd.lm) #0.71

seed.mass.sd.lm=lm(SDc~effectiveness_value+sample_size+treatment_condition,
                  data=no.missing.seed.mass)#

summary(seed.mass.sd.lm) #0.79

pollen.dep.sd.lm=lm(SDc~effectiveness_value+sample_size+treatment_condition,
                   data=no.missing.pollen.dep)#

summary(pollen.dep.sd.lm) #0.81

#impute missing SDs
missing.seed.set$SDc=predict(seed.set.sd.lm,missing.seed.set)
missing.seed.mass$SDc=predict(seed.mass.sd.lm,missing.seed.mass)
missing.pollen.dep$SDc=predict(pollen.dep.sd.lm,missing.pollen.dep)

##combine dataframe
diel.out=rbind.fill(no.missing.all,
                    missing.seed.set,
                    missing.seed.mass,
                    missing.pollen.dep) 

#check imputation types visually
ggplot(data=diel.out,
       aes(x=effectiveness_value,
           y=SDc,fill=impute.type))+
  geom_point(pch=21,size=4,col="black")+
  theme_bw()+
  theme(aspect.ratio=1)+
  ylab("SD")+
  xlab("Mean")+
  facet_wrap(~treatment_effectiveness_metric,scales="free")

diel.out %>% filter(treatment_effectiveness_metric%in%"fruit mass" & SDc >200)

diel.out %>% filter(treatment_effectiveness_metric%in%"fruit set" & SDc >15)

diel.out %>% filter(treatment_effectiveness_metric%in%"seed mass" & SDc >5)

#Widen dataframe
diel.final=diel.out %>%   
  select(study_ID, 
         country,
         Site_period, 
         year_start,
         month_start,
         month_end,
         coordinates_lat,
         coordinates_lon,
         plant_species,
         plant_cultivar, 
         additional_treatment, 
         treatment_effectiveness_metric,
         treatment_condition,
         sample_size,
         effectiveness_value,
         SDc) 

#Imputed values make sense (but definitely some different scales going on)
source('script/DvN-phylogeny-set-up.R')

diel.final.out=diel.final%>% 
  left_join(plant.taxonomy.df.out,by=c("plant_species" = "old_name")) %>% 
  mutate(phylo=gsub(" ","_",accepted_name)) %>% 
  relocate(accepted_name:phylo,.after=plant_species) #%>% glimpse

###manual fix of ipomoea aff.
sum(is.na(diel.final.out$accepted_name)) #7

#select plant_species and accepted_name and filter only NAs in accepted name
diel.final.out %>% 
  select(plant_species,accepted_name) %>% 
  filter(is.na(accepted_name)) %>% 
  distinct(plant_species)

#####add environmental co-variates
list.files("data/")
temperature.range=read.csv("data/DvN_Site_DailyTemperatureRange.csv",row.names = 1)
day.length=read.csv("data/DvN_Site_Daylength.csv",row.names = 1)%>% 
  mutate(year_start=as.integer(year_start)) 

#View(temperature.range)

environment=temperature.range %>% 
  left_join(day.length,
            by=c("study_ID","coordinates_lat","coordinates_lon","country",
                 "year_start","month_start","month_end")) %>% 
  mutate(year_start=as.character(year_start)) 

diel.env.final=diel.final.out %>% 
  left_join(environment,
            by=c("study_ID","coordinates_lat","coordinates_lon","country",
                 "year_start","month_start","month_end")) %>% #glimpse
  relocate(midDate:Daylength,.after=coordinates_lon)

###add traits dataframe to diel.env.final
#source("script/DvN-trait-set-up.R")

traits=read.csv("data/SppTraitMatrix.csv",row.names=1,stringsAsFactors = T)
rownames(traits)=ifelse(rownames(traits)%in%"Ipomoea aff. Marcellia","Ipomoea aff-marcellia",
                        rownames(traits))

rownames(traits)=ifelse(rownames(traits)%in%"Silene latifolia  subsp. alba",
                        "Silene latifolia-alba",
                        rownames(traits))

rownames(traits)=ifelse(rownames(traits)%in%"Castilleja purpurea  var. citrina",
                        "Castilleja purpurea-citrina",
                        rownames(traits))

rownames(traits)=ifelse(rownames(traits)%in%"Castilleja purpurea  var. lindheimeri",
                        "Castilleja purpurea-lindheimeri",
                        rownames(traits))

#check difference between traits rownames and plant species from diel
setdiff(rownames(traits),
        unique(diel.env.final$accepted_name))

setdiff(unique(diel.env.final$accepted_name),
        rownames(traits))
#BRA!

diel.env.final.out=diel.env.final %>% 
  left_join(traits %>% 
              mutate(accepted_name=rownames(traits)),
            by=c("accepted_name")) 

