######
#Calculate species pollination dependency
######

diel.pd=diel.env.final[-c(690,696,702,708),]  %>%  
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
  filter(treatment_condition%in%c("open pollination",
                                  "complete exclusion")) %>% 
  pivot_wider(names_from = treatment_condition,
              values_from=c(sample_size,SDc,effectiveness_value)) %>% 
  ungroup() %>% 
  mutate(effect_ID=1:n()) %>% 
  relocate(effect_ID,.after = study_ID) 

####dependency effect sizes
diel.pd.es=escalc(data=diel.pd,
                    measure="SMD",
                    #day pollination (control)
                    m2i=`effectiveness_value_complete exclusion`, 
                    sd2i = `SDc_complete exclusion`,
                    n2i=`sample_size_complete exclusion`, 
                    
                    #night pollination (treatment)
                    m1i=`effectiveness_value_open pollination`, 
                    sd1i = `SDc_open pollination`,
                    n1i=`sample_size_open pollination`,
                    append=F)

colnames(diel.pd.es)=c("Yi.pd","VI.pd")


diel.pd.out=cbind(diel.pd,diel.pd.es) %>% 
  filter(!is.na(Yi.pd)) %>% 
  group_by(accepted_name,treatment_effectiveness_metric) %>% 
  summarise(pd=mean(Yi.pd))

pd.values=diel.pd.out %>% 
  ungroup() %>% 
  #select(2:4) %>% 
  pivot_wider(names_from = 2,
              values_from = 3)%>% 
  rename("seed_set_pd"=`seed set`,
         "fruit_set_pd"=`fruit set`)

plot(pd.values$seed_set_pd,pd.values$fruit_set_pd)

pd.coverage=diel.open.es.out %>% 
  ungroup() %>% 
  distinct(accepted_name,treatment_effectiveness_metric)%>% 
  select(accepted_name,treatment_effectiveness_metric) %>% 
  mutate(sum=1)%>% 
  pivot_wider(names_from = 2,
              values_from = 3) %>% 
  rename("seed_set_incl"=`seed set`,
         "fruit_set_incl"=`fruit set`) %>% 
  left_join(pd.values)

write.csv(pd.coverage,"data/DvN_Missing_dependency.csv")

##hand pollination imputations

diel.scale.phd=diel.phd.out %>% 
  group_by(treatment_effectiveness_metric) %>%
  mutate(sPhd=as.numeric(scale(p.hd)))%>% 
  ungroup()

diel.fs.pd.out=cbind(diel.pd,diel.pd.es) %>% 
  filter(!is.na(Yi.pd)) %>% 
  filter(treatment_effectiveness_metric%in%"fruit set") %>% 
  group_by(accepted_name,phylo) %>% 
  summarise(pd=mean(Yi.pd))%>% 
  ungroup()

diel.fs.pd.out$sPd=as.numeric(scale(diel.fs.pd.out$pd))

diel.fs.pd.imp.out=diel.fs.pd.out %>% 
  full_join(diel.scale.phd%>% 
              filter(treatment_effectiveness_metric%in%"fruit set") %>% 
              select(-c(treatment_effectiveness_metric)))

diel.fs.pd.imp.out$iPd=ifelse(is.na(diel.fs.pd.imp.out$sPd),
                              diel.fs.pd.imp.out$sPhd,
                              diel.fs.pd.imp.out$sPd)
View(diel.fs.pd.imp.out)


diel.ss.pd.out=cbind(diel.pd,diel.pd.es) %>% 
  filter(!is.na(Yi.pd)) %>% 
  filter(treatment_effectiveness_metric%in%"seed set") %>% 
  group_by(accepted_name,phylo) %>% 
  summarise(pd=mean(Yi.pd))

####List of missing species
diel.open.fs.species
missing.fs.pd=setdiff(diel.open.fs.species,diel.fs.pd.out$phylo) #34
missing.ss.pd=setdiff(diel.open.ss.species,diel.ss.pd.out$phylo) #25

all.missing=unique(c(missing.fs.pd,missing.ss.pd))





##check phylogenetic signal of PD

diel.pd.signal=cbind(diel.pd,diel.pd.es) %>% 
  filter(!is.na(Yi.pd))

diel.pd.species=diel.pd.signal$phylo

diel.pd.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                       diel.pd.species))
diel_pd_vcv <- vcv(compute.brlen(diel.pd.tree), cor = T)

###
diel_pd_mod_1 <- rma.mv(yi = Yi.pd, V = VI.pd, #mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID,
                                            ~ 1 | effect_ID), 
                             data = diel.pd.signal)
summary(diel_pd_mod_1)

diel_pd_mod_2 <- rma.mv(yi = Yi.pd, V = VI.pd, 
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            ~ 1 | phylo),
                             R = list(phylo = diel_pd_vcv), 
                             data = diel.pd.signal)

AIC(diel_pd_mod_1,
    diel_pd_mod_2)


#pd tree
diel.pd.tree

## extract character of interest
pd.yi=diel.pd.out  %>% 
  mutate(phylo=gsub(" ","_",accepted_name)) %>% 
  as.data.frame() 

#match to phylo tips
pd.yi=pd.yi[ order(match(pd.yi$phylo, diel.pd.tree$tip.label)), ]

pd.yi=setNames(pd.yi$pd,
               pd.yi$phylo)

## estimate ancestral state under BM model
fit.BM<-anc.ML(diel.pd.tree,pd.yi)
print(fit.BM)

obj<-contMap(diel.pd.tree,pd.yi,fsize=c(0.5,0.1),lwd=0.5,plot=T)
obj2=setMap(obj,RColorBrewer::brewer.pal(n=3,name="RdBu"))
#ln(Day/Night)
plot(obj2,fsize=c(0.7,1),leg.txt="Pollination dependency",lwd=3)

fitContinous
