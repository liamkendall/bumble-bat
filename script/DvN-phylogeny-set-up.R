####get phylo for plant species

#run set-up and effect size calculation first

#run the following line to find world plant database or WFO.download() to get WFO.data the first
#WorldFlora::WFO.remember(WFO.file = NULL, WFO.data = "WFO.data", WFO.pos = 1)

diel.phy.df=diel.final

plants=data.frame(plant_species=unique(diel.phy.df$plant_species))

plants.wfo=WFO.match(spec.data=plants, 
                     spec.name="plant_species",
                     WFO.data=WFO.data, 
                     counter=1, 
                     verbose=TRUE)

###picks one record (See details)
plants.wfo.one=WFO.one(plants.wfo,priority = "Accepted") %>% 
  select(plant_species,scientificName,family)#

#21 changed
length(setdiff(plants.wfo.one$plant_species,plants.wfo.one$scientificName))

####accepted species names
plants.wfo.one

##trait database
plant.traits.df=diel.raw %>% 
  distinct(study_ID,journal,DOI,plant_species) %>% 
  filter(plant_species%in%diel.final$plant_species) %>% #11 species filtered here 
  left_join(plants.wfo.one)%>% 
  rename(old_name=plant_species,
         accepted_name=scientificName,
         family=family) 

#write.csv(plant.traits.df,"DvN_Plant_Species_Study_List.csv")

#####
plant.taxonomy.df=plants.wfo.one%>% 
  rename(old_name=plant_species,
         accepted_name=scientificName,
         family=family)# %>% 
 # mutate(accepted_name=word(accepted_name,1,2,sep=" "))

###manual fix for aff. marcellia
#plant.taxonomy.df$accepted_name=ifelse(duplicated(plant.taxonomy.df$accepted_name),
#       "Ipomoea aff. marcellia",
#       plant.taxonomy.df$accepted_name)

plant.taxonomy.df$accepted_name=ifelse(plant.taxonomy.df$old_name%in%"Silene alba",
                                       "Silene latifolia  subsp. alba",
                                       plant.taxonomy.df$accepted_name)

plant.taxonomy.df$accepted_name=ifelse(plant.taxonomy.df$old_name%in%"Castilleja citrina",
       "Castilleja purpurea  var. citrina",
       plant.taxonomy.df$accepted_name)

plant.taxonomy.df$accepted_name=ifelse(plant.taxonomy.df$old_name%in%"Castilleja lindheimeri",
       "Castilleja purpurea  var. lindheimeri",
       plant.taxonomy.df$accepted_name)

plant.taxonomy.df$accepted_name=ifelse(plant.taxonomy.df$old_name%in%"Ipomoea Marcellia",
                                       "Ipomoea aff. marcellia",
                                       plant.taxonomy.df$accepted_name)


#####phylogeny

#tree
dn.tree = get_tree(sp_list = plant.taxonomy.df%>% 
                     filter(!duplicated(accepted_name)) %>% 
                     rename(species=accepted_name),
                   taxon = "plant",
                   scenario = "at_basal_node",
                   show_grafted = F)

unique(dn.tree$tip.label)

dn.tree$tip.label=ifelse(dn.tree$tip.label%in%"Silene_latifolia_subsp._alba",
                                       "Silene_latifolia-alba",
                         dn.tree$tip.label)

dn.tree$tip.label=ifelse(dn.tree$tip.label%in%"Castilleja_purpurea_var._citrina",
                                       "Castilleja_purpurea-citrina",
                         dn.tree$tip.label)

dn.tree$tip.label=ifelse(dn.tree$tip.label%in%"Castilleja_purpurea_var._lindheimeri",
                                       "Castilleja_purpurea-lindheimeri",
                         dn.tree$tip.label)

dn.tree$tip.label=ifelse(dn.tree$tip.label%in% "Ipomoea_aff._marcellia",
                         "Ipomoea_aff-marcellia",
                         dn.tree$tip.label)


###configure plant.taxonomy.df
plant.taxonomy.df.out=plant.taxonomy.df
plant.taxonomy.df.out$accepted_name
plant.taxonomy.df.out$accepted_name=ifelse(plant.taxonomy.df.out$accepted_name%in%"Silene latifolia  subsp. alba",
                         "Silene latifolia-alba",
                         plant.taxonomy.df.out$accepted_name)

plant.taxonomy.df.out$accepted_name=ifelse(plant.taxonomy.df.out$accepted_name%in%"Castilleja purpurea  var. citrina",
                         "Castilleja purpurea-citrina",
                         plant.taxonomy.df.out$accepted_name)

plant.taxonomy.df.out$accepted_name=ifelse(plant.taxonomy.df.out$accepted_name%in%"Castilleja purpurea  var. lindheimeri",
                         "Castilleja purpurea-lindheimeri",
                         plant.taxonomy.df.out$accepted_name)

plant.taxonomy.df.out$accepted_name=ifelse(plant.taxonomy.df.out$accepted_name%in% "Ipomoea aff. marcellia",
                         "Ipomoea aff-marcellia",
                         plant.taxonomy.df.out$accepted_name)

