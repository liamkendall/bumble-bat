################################################
###Effect size calculations - bag as control###
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

diel.bag.control=diel.env.final %>%  
  filter(treatment_condition%in%c("complete exclusion")) %>% 
  rename("N_bag"="sample_size",
         "M_bag"="effectiveness_value",
         "SD_bag"="SDc")

diel.bag.out=diel.open.treatments %>% 
  filter(!label%in%bad.studies$label) %>% 
  left_join(diel.bag.control %>% select(-treatment_condition),
            by = c("study_ID", "country", "Site_period", "year_start", "month_start",
                   "month_end", "coordinates_lat", "coordinates_lon", "midDate", 
                   "JDay_midDate", "Lat_int", "Daylength", "plant_species", 
                   "accepted_name", "family", "phylo", "plant_cultivar", "additional_treatment",
                   "treatment_effectiveness_metric", "DTR")) %>% 
  filter(!is.na(M_bag)) %>% 
  filter(!is.na(N_bag))%>% 
  filter(!is.na(SD_bag))

table(diel.bag.out$M_bag)

####day-night effect sizes
#effect sizes
diel.bag.es=escalc(data=diel.bag.out,
                    measure="SMD",
                    #complete exclusion (control)
                    m2i=M_bag, 
                    sd2i = `SD_bag`,
                    n2i=`N_bag`, 
                    
                    #pollination (treatment)
                    m1i=`M_poll`, 
                    sd1i = `SD_poll`,
                    n1i=`N_poll`,
                    append=F)

###cbind the g's
diel.bag.es.out=cbind(diel.bag.out,
                       diel.bag.es) 
