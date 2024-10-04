####fix taxonomy for plant species

#run during set-up

WFO.download(timeout=100^9)

diel.tax.df=diel.final

plants=data.frame(plant_species=unique(diel.tax.df$plant_species))

plants.wfo=WFO.match(spec.data=plants, 
                     spec.name="plant_species",
                     WFO.data=WFO.data, 
                     counter=1, 
                     verbose=TRUE)

###picks one record (See details)
plants.wfo.one=WFO.one(plants.wfo,priority = "Accepted") %>% 
  select(plant_species,scientificName,family)#

#20 changed
length(setdiff(plants.wfo.one$plant_species,plants.wfo.one$scientificName))

####accepted species names
plants.wfo.one

#####
plant.taxonomy.df=plants.wfo.one%>% 
  rename(old_name=plant_species,
         accepted_name=scientificName,
         family=family)

#Check for differences between the two dataframes (due to Worldflora formatting for database gleaning)
setdiff(plant.taxonomy.df$old_name,diel.final$plant_species)
setdiff(diel.final$plant_species,plant.taxonomy.df$old_name)

#Fix mismatched names
plant.taxonomy.df$old_name=ifelse(plant.taxonomy.df$old_name%in%"Ipomoea Marcellia",
                                  "Ipomoea aff. Marcellia",
                                  plant.taxonomy.df$old_name)

plant.taxonomy.df$old_name=ifelse(plant.taxonomy.df$old_name%in%"Fragaria ×ananassa",
                                      "Fragaria x ananassa",
                                  plant.taxonomy.df$old_name)

plant.taxonomy.df$accepted_name=ifelse(plant.taxonomy.df$old_name%in% "Ipomoea aff. Marcellia",
                                           "Ipomoea aff-marcellia",
                                       plant.taxonomy.df$accepted_name)

plant.taxonomy.df$accepted_name=ifelse(plant.taxonomy.df$accepted_name%in%"Fragaria × ananassa",
                                           "Fragaria ananassa",
                                       plant.taxonomy.df$accepted_name)

plant.taxonomy.df$accepted_name=ifelse(plant.taxonomy.df$old_name%in%"Silene alba",
                                       "Silene latifolia-alba",
                                       plant.taxonomy.df$accepted_name)

#done


