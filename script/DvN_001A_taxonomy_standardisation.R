#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "20/11/2024"

###script: 001A - Standardisation of plant taxonomy
#runs within script: DvN_001

#Download plant taxonomy database
##This function can take a very long time to download. If you have issues downloading it using this function, we suggest that you downloaded independently from https://www.worldfloraonline.org/downloadData and use the WFO.remember function. 
WFO.download(timeout=100^9)

diel.tax.df=diel.final

plants=data.frame(plant_species=unique(diel.tax.df$plant_species))

#match plant species in df to WFO database
plants.wfo=WFO.match(spec.data=plants, 
                     spec.name="plant_species",
                     WFO.data=WFO.data, 
                     counter=1, 
                     verbose=TRUE)

###picks one record (See details)
plants.wfo.one=WFO.one(plants.wfo,priority = "Accepted") %>% 
  select(plant_species,scientificName,family)#

#20 species changed
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
