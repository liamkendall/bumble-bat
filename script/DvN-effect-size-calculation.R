##############################
###Effect size calculations###
##############################

diel.es=diel.final %>%  
  group_by(study_ID, 
           Site_period, 
           year_start,
           month_start,
           coordinates_lat,
           coordinates_lon,
           plant_species,
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
                  append=F)

colnames(day.night.es)=c("SMD.dn","VI.dn")

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
             append=F)

colnames(open.night.es)=c("SMD.on","VI.on")

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
                     append=F)

colnames(open.day.es)=c("SMD.od","VI.od")

###open vs. bag (TBC)

###cbind the g's
diel.es.out=cbind(diel.es,
                  day.night.es,
                  open.night.es,
                  open.day.es)

##then subset at your leisure
day.night.df=diel.es.out %>% 
  filter(!is.na(SMD.dn))

open.day.df=diel.es.out %>% 
  filter(!is.na(SMD.od))

open.night.df=diel.es.out %>% 
  filter(!is.na(SMD.on))
