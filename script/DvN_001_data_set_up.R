####Kendall & Nicholson 2024. Pollination across the diel cycle
####Script 1. Data input and effect size calculation

####Libraries
###Wrangle data to wide format
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

###phylo/taxonomy packages
library(WorldFlora)
library(phytools)
#library(rtrees)
#library(V.PhyloMaker2)

###meta-analysis
library(metafor)

###plotting
library(orchaRd)
library(ggplot2)

#bespoke functions
#binomial summary stats, simplified orchard plots
source('script/DvN_additional_functions.R')

#load raw data
diel.raw=read.csv("data/DvN_LiteratureDataClean_2024-09-13.csv",#sep=";",dec=".",
                  row.names = 1,
                  encoding = "UTF-8")

#take relevant columns
diel=diel.raw %>% 
  dplyr::select(study_ID,country,Site_period,coder_initials,
         coordinates_lat,coordinates_lon,
         year_start:night_pollination_hrs,
         plant_species,plant_cultivar,plant_crop,
         treatment_condition,
         additional_treatment,
         treatment_effectiveness_metric,
         sample_size:effect_error_units) %>% 
  filter(!sample_size<=1) %>% #remove sample sizes of 1
  filter(!is.na(sample_size)==T) %>% #remove NAs
  filter(!is.na(effectiveness_value)==T)

table(diel$plant_crop) #1356,229
  
####check stuff
table(diel$treatment_condition)
table(diel$treatment_effectiveness_metric)
table(diel$treatment_condition,diel$treatment_effectiveness_metric)

#####Summarising missing errors
sum(is.na(diel$effect_error))
table(diel$effect_error_units)

########
#set up required imputation / transformation columns
#binary = use binomial
#SED = SE to SD transformation
#zero = zero error as complete fruit set (probability of 0 or 1)
#CI91t and CI95 = SD from confidence interval - if sample size less than 60, use t-distribution otherwise
#CV = use linear regression (See below)
#none = no imputation required

diel=diel %>% mutate(impute.type=case_when(treatment_effectiveness_metric%in%"fruit set" &
                          is.na(effect_error) == T & 
                          is.na(sample_size) == F &
                          effectiveness_value<=1& 
                          !effectiveness_value==0  &
                          !effect_error_units%in%"CI95"~ "binary",
                          
                          is.na(effect_error) == F & 
                            effect_error > 0 &
                            is.na(sample_size) == F & 
                            effect_error_units == "SE" ~  "SED",
                          
                            treatment_effectiveness_metric%in%"fruit set" &
                            effectiveness_value ==0|
                            treatment_effectiveness_metric%in%"fruit set" &
                            effectiveness_value ==1 ~ "zero", 
                          
                            effect_error_units == "CI95" &
                            sample_size < 60 ~
                              "CI95t",
                          effect_error_units == "CI95" &
                            sample_size > 60 ~
                            "CI95",
                            is.na(effect_error)==T ~ "CV",
                          effect_error_units=="min-max" ~ "CV",
                            .default = "none")) 

table(diel$impute.type)

#summary of imputation type
table(diel$impute.type,diel$coder_initials) #CI95, 40, 20 

#by treatment_effectiveness_metric
diel %>% group_by(
  impute.type,treatment_effectiveness_metric) %>% 
  summarise(sum_na = sum(is.na(effect_error)))

######Imputation
#those conformable to binary arrays ("binary")
#from SEs ("SED")

diel.sd=diel %>% 
  mutate(t.dist=qt(0.975,sample_size-1)) %>% #get t-distribution
  rowwise() %>% 
  mutate(SDc = ifelse(impute.type%in%"binary", #binary function call
                      binary.summary(effectiveness_value,sample_size)[3], #returns just SD
                      
                      ifelse(impute.type%in%"SED", #else SE to SD transformation
                             effect_error*sqrt(sample_size),
                      effect_error)))%>% 
  #confidence intervals to SD
  mutate(SDc=ifelse(impute.type%in%"CI95",
                    sqrt(sample_size)*(effect_error_U-       
                                         effect_error_L)/3.92,
                    SDc)) %>% 
  mutate(SDc=ifelse(impute.type%in%"CI95t",
                    sqrt(sample_size)*(effect_error_U-       
                                         effect_error_L)/t.dist,
                    SDc)) %>% 
  ###make sure 1s and 0s have 0 variance fruit set
  mutate(SDc=ifelse(impute.type%in%"zero",0,SDc)) %>% 
  #colwise() %>% 
  dplyr::select(-t.dist)

###remaining missingness summary
diel.sd %>% group_by(impute.type,treatment_effectiveness_metric) %>% 
  summarise(sum_na = sum(is.na(SDc)))

#23 missing from seed set 
#6 missing from pollen deposition
#5 missing from seed mass

###use Bishop & Nakagawa (2021) method for imputing missing SDs

#When a measure of variance was not obtainable, 
#we imputed standard deviation based on the fitted 
#relationship between mean, SD and n in the available data 

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
           y=SDc,
           fill=impute.type,
           size=(sample_size)))+
  geom_point(pch=21,col="black")+
  theme_bw()+
  theme(aspect.ratio=1)+
  ylab("SD")+
  xlab("Mean")+
  facet_wrap(~treatment_effectiveness_metric,scales="free")

#Widen dataframe
diel.final=diel.out %>%   
  select(study_ID, 
         country,
         Site_period, 
         year_start,
         month_start,
         month_end,
         day_pollination_hrs,
         night_pollination_hrs,
         coordinates_lat,
         coordinates_lon,
         plant_species,
         plant_cultivar,
         plant_crop,
         additional_treatment, 
         treatment_effectiveness_metric,
         treatment_condition,
         sample_size,
         effectiveness_value,
         SDc) 


#run script to check and validate plant names
source('script/DvN_001A_taxonomy_standardisation.R')

#join to main dataframe
diel.final.out=diel.final%>% 
  left_join(plant.taxonomy.df,by=c("plant_species" = "old_name")) %>% 
  mutate(phylo=gsub(" ","_",accepted_name)) %>% 
  relocate(accepted_name:phylo,.after=plant_species) #%>% glimpse

###mcheck for missing names
sum(is.na(diel.final.out$accepted_name)) #0

#select plant_species and accepted_name and filter only NAs in accepted name
diel.final.out %>% 
  select(plant_species,accepted_name) %>% 
  filter(is.na(accepted_name)) %>% 
  distinct(plant_species)

#####add environmental co-variates
list.files("data/")
temperature.range=read.csv("data/DvN_Site_DailyTemperatureRange.csv",row.names = 1)

temperature.range$study_ID %>% sort

###fix bad study names
temperature.range$study_ID=ifelse(temperature.range$study_ID=="JÃ¼rgens_2014","Jürgens_2014",temperature.range$study_ID)
temperature.range$study_ID=ifelse(temperature.range$study_ID=="Larrea-AlcÃ¡zar_2011","Larrea-Alcázar_2011",temperature.range$study_ID)


day.length=read.csv("data/DvN_Site_Daylength.csv",row.names = 1)%>% 
  mutate(year_start=as.integer(year_start)) 

###fix bad study names
day.length$study_ID=ifelse(day.length$study_ID=="J\xfcrgens_2014",
                           "Jürgens_2014",day.length$study_ID)

day.length$study_ID=ifelse(day.length$study_ID=="Larrea-Alc\xe1zar_2011",
                           "Larrea-Alcázar_2011",day.length$study_ID)

#fix one study with bad year
day.length$year_start=ifelse(day.length$study_ID=="Otero-Arnaiz_2003",2003,day.length$year_start)


elevation=read.csv("data/DvN_Sites_elevation.csv",row.names = 1)

###fix bad study names
elevation$study_ID=ifelse(elevation$study_ID=="J\xfcrgens_2014",
                           "Jürgens_2014",elevation$study_ID)
elevation$study_ID=ifelse(elevation$study_ID=="Larrea-Alc\xe1zar_2011",
                           "Larrea-Alcázar_2011",elevation$study_ID)

#combine all env into one df
environment=temperature.range %>% 
  full_join(day.length,
            by=c("study_ID","coordinates_lat","coordinates_lon","country",
                 "year_start","month_start","month_end")) %>% 
  mutate(year_start=as.integer(year_start))  %>% 
  full_join(elevation,
            by=c("study_ID","coordinates_lat","coordinates_lon","country"))

#check merge - Nas for each column
environment %>% 
  summarise_all(funs(sum(is.na(.))))

#merge with diel dataframe
diel.env.final=diel.final.out %>% 
  left_join(environment,
            by=c("study_ID","coordinates_lat","coordinates_lon","country",
                 "year_start","month_start","month_end")) %>% #glimpse
  relocate(midDate:elevation,.after=coordinates_lon)

diel.env.final %>% 
  summarise_all(funs(sum(is.na(.))))

####load traits dataframe
traits=read.csv("data/SppTraitData_clean.csv",row.names=1,stringsAsFactors = T)#[,-1]
rownames(traits)=traits$plant_species

#check names
setdiff(rownames(traits),
        unique(diel.env.final$accepted_name))

setdiff(unique(diel.env.final$accepted_name),
        rownames(traits))
#BRA!

#combine 
diel.env.final.out=diel.env.final %>% 
  left_join(traits %>% 
              mutate(plant_height_midpoint_m=ifelse(is.na(plant_height_midpoint_m),0,plant_height_midpoint_m)) %>% 
              mutate(accepted_name=plant_species) %>% 
              select(-plant_species),
            by=c("accepted_name")) 


#check for missing values
sum(is.na(diel.env.final.out$plant_species)) #0

#count sum NAs in each column
diel.env.final.out %>% 
  summarise_all(~sum(is.na(.)))

save(diel.env.final.out,file="output/DvN_raw_dataframe_pre_effect_size_calculation Sep16.RData")
