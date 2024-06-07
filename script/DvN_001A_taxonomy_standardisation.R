####fix taxonomy for plant species

#run during set-up

#run the following line to find world plant database or WFO.download() to get WFO.data the first
WFO.download(timeout=100^9)
#WorldFlora::WFO.remember(WFO.file = "WFO_Backbone.zip", WFO.data = "WFO.data", WFO.pos = 1)

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

##trait database
#plant.traits.df=diel.raw %>% 
#  distinct(study_ID,journal,DOI,plant_species) %>% 
#  filter(plant_species%in%diel.final$plant_species) %>% #11 species filtered here 
#  left_join(plants.wfo.one)%>% 
#  rename(old_name=plant_species,
#         accepted_name=scientificName,
#         family=family) 

#write.csv(plant.traits.df,"DvN_Plant_Species_Study_List.csv")

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


