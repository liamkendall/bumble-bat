####get phylo for plant species

#run set-up and effect size calculation first

#run the following line to find world plant database or WFO.download() to get WFO.data the first
#WFO.remember(WFO.file = NULL, WFO.data = "WFO.data", WFO.pos = 1)

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
         family=family) 

###manual fix for aff. marcellia
plant.taxonomy.df$accepted_name=ifelse(duplicated(plant.taxonomy.df$accepted_name),
       "Ipomoea aff. marcellia",
       plant.taxonomy.df$accepted_name)

####
plants.f=plant.taxonomy.df %>% 
  mutate(species=gsub(" ","_",accepted_name)) %>% 
  mutate(genus=word(accepted_name,1)) %>% 
  select(species,genus,family)

#get phylogeny
dn.tree = get_tree(sp_list = plants.f,
                     taxon = "plant",
                     scenario = "at_basal_node",
                     show_grafted = F)

#45 species added at genus level (*) 
#10 species added at family level (**) 

#bespoke trees for effect sizes

#day night tree
day.night.species=gsub(" ","_",unique(day.night.df$plant_species))
day.night.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                         day.night.species))

#open day tree
open.day.species=gsub(" ","_",unique(open.day.df$plant_species))
open.day.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                        open.day.species))

#open night tree
open.night.species=gsub(" ","_",unique(open.night.df$plant_species))
open.night.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                        open.night.species))

##may be some manual name matching to do with the Castilleja's

##convert to VCV matrices when model fitting


##
#random phytools visualisation stuff
#
#####
#plot(dn.tree)
#
##fruit set tree
#fs.tree=drop.tip(dn.tree,setdiff())
##seed set tree
#
#
### extract character of interest
#species.yi=diel.es %>% 
#  group_by(plant_species,treatment_effectiveness_metric) %>% 
#  summarise(yi=mean(yi)) %>% 
#  mutate(plant_species=gsub(" ","_",plant_species)) %>% 
#  as.data.frame()
#
##check names match
#drop.tip
#setdiff(species.yi$plant_species,dn.tree$tip.label)
#setdiff(dn.tree$tip.label,species.yi$plant_species)
#
##match to phylo tips
#species.tips.yi=species.yi[ order(match(species.yi$plant_species, dn.tree$tip.label)), ]
#
#species.tips.yi=setNames(species.tips.yi$yi,
#                         diel.trait$plant_species)
#
### estimate ancestral state under BM model
#fit.BM<-anc.ML(dn.tree,species.tips.yi)
#print(fit.BM)
#
#obj<-contMap(dn.tree,species.tips.yi,fsize=c(0.5,0.1),lwd=0.5,plot=T)
#obj2=setMap(obj,RColorBrewer::brewer.pal(n=3,name="RdBu"))
##ln(Day/Night)
#plot(obj2,fsize=c(0.7,1),leg.txt="Fruit set",lwd=3)

