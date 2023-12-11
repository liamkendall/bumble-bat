##############################
###Effect size calculations###
##############################

diel.es=diel.env.final.out %>%  
  group_by(study_ID, 
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
  ungroup() 

##non-nestedness of D v. N comes from site_periods/other treatments 
##where they did something else there (hand-pollination) but not D v. N

####day-night effect sizes
#effect sizes
day.night.es=escalc(data=diel.es,
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

###open vs. night 
open.night.es=escalc(data=diel.es,
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

###open vs. day 
open.day.es=escalc(data=diel.es,
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

###bag vs. day 
bag.day.es=escalc(data=diel.es,
                   measure="SMD",
                   #open pollination (control)
                   m2i=`effectiveness_value_complete exclusion`, 
                   sd2i = `SDc_complete exclusion`,
                   n2i=`sample_size_complete exclusion`, 
                   
                   #day pollination (treatment)
                   m1i=`effectiveness_value_day pollination`, 
                   sd1i = `SDc_day pollination`,
                   n1i=`sample_size_day pollination`,
                   append=T) %>% 
  filter(!is.na(yi)) %>% 
  mutate(treatment="bag_day")


###bag vs. night
bag.night.es=escalc(data=diel.es,
                  measure="SMD",
                  #open pollination (control)
                  m2i=`effectiveness_value_complete exclusion`, 
                  sd2i = `SDc_complete exclusion`,
                  n2i=`sample_size_complete exclusion`, 
                  
                  #day pollination (treatment)
                  m1i=`effectiveness_value_night pollination`, 
                  sd1i = `SDc_night pollination`,
                  n1i=`sample_size_night pollination`,
                  append=T) %>% 
  filter(!is.na(yi)) %>% 
  mutate(treatment="bag_night")

###bag vs. open
bag.open.es=escalc(data=diel.es,
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

###cross vs. open
hand.open.es=escalc(data=diel.es,
                   measure="SMD",
                   #open pollination (control)
                   m2i=`effectiveness_value_open pollination`, 
                   sd2i = `SDc_open pollination`,
                   n2i=`sample_size_open pollination`, 
                   
                   #hand pollination (treatment)
                   m1i=`effectiveness_value_hand pollination`, 
                   sd1i = `SDc_hand pollination`,
                   n1i=`sample_size_hand pollination`,
                   append=T) %>% 
  filter(!is.na(yi)) %>% 
  mutate(treatment="open_hand")

###cbind the g's
diel.es.out=rbind(day.night.es,
                  open.night.es,
                  open.day.es,
                  bag.night.es,
                  bag.day.es,
                  bag.open.es,
                  hand.open.es)%>% 
  mutate(effect_ID=1:n()) %>% 
  relocate(effect_ID,.after = study_ID)

diel.es.out %>% glimpse

###effect size dataframes is in long format and then saved as an rData file in the data folder
es.list=list(diel.es.out,
             dn.tree,
             traits.out)

#name list elements
names(es.list)=c("dvn_effects","phylo","traits")
save(es.list,file="data/effect_size_dataframes.rData")
