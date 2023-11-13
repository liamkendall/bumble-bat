####Wrangel data to wide format
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

#plotting
library(ggplot2)

source('script/DvN-functions.R')

#raw dataframe
diel.raw=read.csv("data/DvN_LiteratureDataClean_2023-11-13.csv",sep=",",dec=".",
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
  filter(!is.na(effectiveness_value)==T) %>% 
  filter(treatment_effectiveness_metric%in%c("fruit set",
                                             "seed set"))

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
diel.sd %>% group_by(impute.type,treatment_effectiveness_metric) %>% 
  summarise(sum_na = sum(is.na(SDc)))

###use Bishop & Nakagawa (2021) method for imputing missing SDs

#When a measure of variance was not obtainable, 
#we imputed standard deviation based on the fitted 
#relationship between mean, SD and n in the available data 
#(91 data rows imputed; adjusted-R2 of these models were between 0.29 and 0.46).

missing.seed=diel.sd %>% 
  filter(is.na(SDc))
no.missing.seed=diel.sd %>% 
  filter(!is.na(SDc))

sd.lm=lm(SDc~effectiveness_value*sample_size*treatment_condition,
         data=no.missing.seed)#

summary(sd.lm)

#With treatment condition: Adjusted R-squared:  0.7275

###simpler model

#sd.lm2=lm(SDc~effectiveness_value*sample_size,
#         data=no.missing.seed)#
#
#summary(sd.lm2)

#Without treatment condition: Adjusted R-squared:  0.6353

#impute missing SDs
missing.seed$SDc=predict(sd.lm,missing.seed)

##combine dataframe
diel.out=rbind.fill(missing.seed,
                    no.missing.seed) 

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

#Imputed values make sense (but definitely some different scales going on)

#Widen dataframe
diel.final=diel.sd %>%   
  select(study_ID, 
         Site_period, 
         year_start,
         month_start,
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




