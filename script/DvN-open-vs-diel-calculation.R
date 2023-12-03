################################################
###Effect size calculations - open as control###
################################################

diel.open.treatments=diel.env.final %>%  
  filter(treatment_condition%in%c("day pollination",
                                  "night pollination")) %>% 
  rename("N_poll"="sample_size",
         "M_poll"="effectiveness_value",
         "SD_poll"="SDc") %>% 
  mutate(effect_ID=1:n())%>% #dim
  group_by(study_ID,
           Site_period, 
           year_start,
           month_start,
           coordinates_lat,
           coordinates_lon,
           accepted_name,
           plant_cultivar, 
           additional_treatment,
           treatment_effectiveness_metric) %>%
  mutate(label = cur_group_id())

###check studies have both day and night
diel.study.check = diel.open.treatments %>% 
  group_by(study_ID,
           label,
           Site_period, 
           year_start,
           month_start,
           coordinates_lat,
           coordinates_lon,
           accepted_name,
           plant_cultivar, 
           additional_treatment,
           treatment_effectiveness_metric) %>% 
  summarise(n=n())

bad.studies=diel.study.check %>% 
  filter(!n%in%c(2))

diel.open.control=diel.env.final %>%  
filter(treatment_condition%in%c("open pollination")) %>% 
  rename("N_open"="sample_size",
         "M_open"="effectiveness_value",
         "SD_open"="SDc")

diel.open.out=diel.open.treatments %>% 
  filter(!label%in%bad.studies$label) %>% 
  left_join(diel.open.control %>% select(-treatment_condition),
            by = c("study_ID", "country", "Site_period", "year_start", "month_start",
                   "month_end", "coordinates_lat", "coordinates_lon", "midDate", 
                   "JDay_midDate", "Lat_int", "Daylength", "plant_species", 
                   "accepted_name", "family", "phylo", "plant_cultivar", "additional_treatment",
                   "treatment_effectiveness_metric", "DTR")) %>% 
  filter(!is.na(M_open)) %>% 
  filter(!is.na(N_open))%>% 
  filter(!is.na(SD_open))

table(diel.open.out$treatment_condition)

####day-night effect sizes
#effect sizes
diel.open.es=escalc(data=diel.open.out,
                    measure="SMD",
                    #open pollination (control)
                    m2i=`M_open`, 
                    sd2i = `SD_open`,
                    n2i=`N_open`, 
                    
                    #night pollination (treatment)
                    m1i=`M_poll`, 
                    sd1i = `SD_poll`,
                    n1i=`N_poll`,
                    append=F)


###cbind the g's
diel.open.es.out=cbind(diel.open.out,
                       diel.open.es) 

