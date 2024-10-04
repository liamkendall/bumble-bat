#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "05/06/2024"

###script: 002 - Effect size (SMD) calculation

##############################
###Effect size calculations###
##############################

#load dataframe from DvN_001
load(file="output/DvN_raw_dataframe_pre_effect_size_calculation Sep16.RData")
diel.env.final.out

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

####day-night effect sizes
#effect sizes
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

###open vs. night 354
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

###open vs. day  346
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

###cbind the g's
diel.es.out=rbind.fill(day.night.es,
                  open.night.es,
                  open.day.es,
                  bag.open.es)%>% 
  mutate(effect_ID=1:n()) %>% 
  relocate(effect_ID,.after = study_ID)

###effect size dataframes is in long format and then saved as an rData file in the data folder
es.list=list(diel.es.out,
             environment,
             traits)

#name list elements
names(es.list)=c("dvn_effects","environment","traits")
save(es.list,file="output/DvN_effect_size_dataframe Sep16.rData")

