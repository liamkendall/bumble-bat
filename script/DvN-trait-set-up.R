##trait matrix

#read trait file from data folder
traits=read.csv("data/SppTraitMatrix.csv",row.names=1,stringsAsFactors = T)
rownames(traits)=ifelse(rownames(traits)%in%"Ipomoea aff. Marcellia","Ipomoea aff-marcellia",
       rownames(traits))

rownames(traits)=ifelse(rownames(traits)%in%"Silene latifolia  subsp. alba",
                         "Silene latifolia-alba",
                         rownames(traits))

rownames(traits)=ifelse(rownames(traits)%in%"Castilleja purpurea  var. citrina",
                         "Castilleja purpurea-citrina",
                         rownames(traits))

rownames(traits)=ifelse(rownames(traits)%in%"Castilleja purpurea  var. lindheimeri",
                         "Castilleja purpurea-lindheimeri",
                         rownames(traits))

#check difference between traits rownames and plant species from diel
setdiff(rownames(traits),
        unique(diel.env.final$accepted_name))

setdiff(unique(diel.env.final$accepted_name),
        rownames(traits))
#BRA!

#####calculate gawdis for trait matrix
install.packages("gawdis")
library(gawdis)

#calculate gawdis for trait matrix
#trait.distGow <- gowdis(traits)
trait.distGaw <- gawdis(traits)

#pcoa(trait.distGaw) %>% biplot()

#Consider removing traits: lifespan breeding_system nectar odour color
table(traits$flower_symmetry)
table(traits$lifespan)
table(traits$life_form)
table(traits$flower_shape)
table(traits$breeding_system)
table(traits$bloom_period_simple)
table(traits$nectar)
table(traits$odour)
table(traits$color)

#cluster distance matrix, plot the clusters and then cut into groups
library(ape)
trait.distGaw.clust=hclust(as.dist(trait.distGaw),method="ward.D2")

plot(trait.distGaw.clust)

#cut into seven groups (for example), add group ID to traits dataframe
trait.groups=cutree(trait.distGaw.clust,k=4)
traits$group=trait.groups

#select numeric columns and summarise by group
num.trait.groups=traits %>% 
  select_if(is.numeric) %>% 
  group_by(group) %>% 
  summarise_all(mean)

#select categorical columns and summarise by group
cat.trait.groups=traits %>% 
  select(flower_symmetry:color,group) %>% 
  group_by(group) %>% 
  summarise_all(function(x) names(which.max(table(x))))

#combine numeric and categorical trait summary dfs
trait.groups=cat.trait.groups %>% 
  left_join(num.trait.groups,by="group")
