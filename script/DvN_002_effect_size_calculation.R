#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "20/11/2024"

###script: 002 - Effect size (SMD) calculation

##############################
###Effect size calculations###
##############################

#libraries (CRAN)
library(plyr)
library(dplyr)
library(tidyr)
library(metafor)

#load dataframe from DvN_001
load(file="output/DvN_raw_dataframe_pre_effect_size_calculation final.RData")

#loaded file
diel.env.final.out

#rename file
diel.es=diel.env.final.out

#widen dataframe
diel.esc.1=diel.es %>% group_by(study_ID, 
           Site_period, 
           year_start,
           month_start,
           coordinates_lat,
           coordinates_lon,
           plant_species,
           accepted_name,
           plant_cultivar, 
           additional_treatment) %>% 
  pivot_wider(names_from = treatment_condition,
              values_from=c(sample_size,SDc,effectiveness_value)) %>% 
  ungroup() #%>% 

####calculate effect sizes
# day vs. night

day.night.es=escalc(data=diel.esc.1,
                  measure="SMD",
                  #day pollination (control)
                  m2i=`effectiveness_value_day pollination`, 
                  sd2i = `SDc_day pollination`,
                  n2i=`sample_size_day pollination`, 
                  
                  #night pollination (treatment)
                  m1i=`effectiveness_value_night pollination`, 
                  sd1i = `SDc_night pollination`,
                  n1i=`sample_size_night pollination`,
                  append=T) %>% 
  filter(!is.na(yi)) %>% 
  mutate(treatment="day_night")

#open vs. night
open.night.es=escalc(data=diel.esc.1,
             measure="SMD",
             #open pollination (control)
             m2i=`effectiveness_value_open pollination`, 
             sd2i = `SDc_open pollination`,
             n2i=`sample_size_open pollination`, 
             
             #night pollination (treatment)
             m1i=`effectiveness_value_night pollination`, 
             sd1i = `SDc_night pollination`,
             n1i=`sample_size_night pollination`,
             append=T) %>% 
  filter(!is.na(yi)) %>% 
  mutate(treatment="open_night")

#open vs. day 
open.day.es=escalc(data=diel.esc.1,
                     measure="SMD",
                     #open pollination (control)
                     m2i=`effectiveness_value_open pollination`, 
                     sd2i = `SDc_open pollination`,
                     n2i=`sample_size_open pollination`, 
                     
                     #day pollination (treatment)
                     m1i=`effectiveness_value_day pollination`, 
                     sd1i = `SDc_day pollination`,
                     n1i=`sample_size_day pollination`,
                     append=T) %>% 
  filter(!is.na(yi)) %>% 
  mutate(treatment="open_day")

###bag vs. open
bag.open.es=escalc(data=diel.esc.1,
                    measure="SMD",
                    #bag (control)
                    m2i=`effectiveness_value_complete exclusion`, 
                    sd2i = `SDc_complete exclusion`,
                    n2i=`sample_size_complete exclusion`, 
                    
                    #open pollination (treatment)
                    m1i=`effectiveness_value_open pollination`, 
                    sd1i = `SDc_open pollination`,
                    n1i=`sample_size_open pollination`,
                    append=T) %>% 
  filter(!is.na(yi)) %>% 
  mutate(treatment="bag_open")

###cbind the effect size dataframes
diel.es.out=rbind.fill(day.night.es,
                  open.night.es,
                  open.day.es,
                  bag.open.es)%>% 
  mutate(effect_ID=1:n()) %>% 
  relocate(effect_ID,.after = study_ID)

#reload environment and trait dataframes
load(file="output/DvN_trait data final.RData")
load(file="output/DvN_environmental data final.RData")

###place effect sizes, traits and environments into a list
es.list=list(diel.es.out,
             environment,
             traits)

#name list elements
names(es.list)=c("dvn_effects","environment","traits")

#save list
save(es.list,file="output/DvN_effect_size_dataframe final.rData")

