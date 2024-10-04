#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "05/06/2024"

###script: 005 - Model fitting with metafor

###load from DvN_003_analysis_dataframe_setup.R

load(file="output/DvN_analysis_dataframe Sep16.RData")
diel.all.diffs
load(file="output/DvN_plant_phylogeny Sep16.RData")
diel.tree.out



##################
##Diel pollination differences: overall meta-analytic model
#for loop across each treatment type: day vs. night, open vs. night, open vs. day, open vs. closed (*dependency)
##################

treatment.levels=unique(diel.all.diffs$treatment)
overall.list=list()

for(i in treatment.levels){
  
  out.df=diel.all.diffs %>% 
    filter(treatment %in% i)
  
  out.species=unique(out.df$phylo)
  out.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                           out.species))
  out.vcv <- vcv(out.tree, cor = T)
  
  overall.full <- rma.mv(yi = yi, 
                         V = vi,
                         mods=~exposure,
                         random = list(~1|treatment_effectiveness_metric,
                                       ~1|study_ID,
                                       ~1|effect_ID,
                                       ~1|phylo,
                                       ~1|accepted_name),
                         R = list(phylo = out.vcv),
                         data = out.df)
  
  overall.list[[i]]=overall.full
  
}

save(overall.list,file="output/DvN_meta_analytical_models Sep16.RData")

#######################
##Meta-analysis: Pollination dependency of plant species
######################

####pollination dependency dataframe
poll.dep.es

pd.species=unique(poll.dep.es$phylo)
pd.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                        pd.species))
pd.vcv <- vcv(pd.tree, cor = T)

pd.full <- rma.mv(yi = yi, 
                   V = vi,
                  # mods = ~ treatment_effectiveness_metric+exposure,
                   random = list(
                     ~1|study_ID,
                     ~1|effect_ID,
                     ~1|phylo,
                     ~1|accepted_name,
                     ~1|treatment_effectiveness_metric),
                   R = list(phylo = pd.vcv),
                   data = poll.dep.es)

save(pd.full,file="output/DvN_pollination_dependency_meta_model Sep16.RData")

######################
##Meta-regression: Diel pollination difference models as a function of pollination effectiveness metric
######################

treatment.levels=unique(diel.all.diffs$treatment)

eff.list=list()

for(i in treatment.levels){
  
  out.df=diel.all.diffs %>% 
    filter(treatment %in% i)
  
  out.species=unique(out.df$phylo)
  out.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                           out.species))
  out.vcv <- vcv(out.tree, cor = T)
  
  eff.full <- rma.mv(yi = yi, 
                     V = vi,
                     mods = ~ treatment_effectiveness_metric+exposure,
                     random = list(
                       ~1|study_ID,
                       ~1|effect_ID,
                       ~1|phylo,
                       ~1|accepted_name),
                     R = list(phylo = out.vcv),
                     data = out.df)
  
  eff.list[[i]]=eff.full
  
}

save(eff.list,file="output/DvN_meta_regression_metric_models Sep16.RData")


######################
##Meta-regression: Day vs. night pollination as a function of pollination dependency
######################

dn.pd.es=diel.all.diffs %>% 
  filter(treatment == "day_night" & !is.na(poll.dep))

dn.pd.species=unique(dn.pd.es$phylo)
dn.pd.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                           dn.pd.species))

dn.pd.vcv <- vcv(dn.pd.tree, cor = T)

DvN.PD.mod0 <- rma.mv(yi,vi,mods = ~ poll.dep+exposure,
                      random = list(~1|treatment_effectiveness_metric,
                                    ~1|study_ID,
                                    ~1|effect_ID,
                                    ~1|phylo,
                                    ~1|accepted_name),
                      R = list(phylo = dn.pd.vcv),
                      data = dn.pd.es)

DvN.PD.mod <- rma.mv(yi,vi,mods = ~ poll.dep + I(poll.dep^2)+exposure,
                     random = list(~1|treatment_effectiveness_metric,
                                   ~1|study_ID,
                                   ~1|effect_ID,
                                   ~1|phylo,
                                   ~1|accepted_name),
                     R = list(phylo = dn.pd.vcv),
                     data = dn.pd.es)

poll.dep.mods=list(`poll.dep day_night linear`=DvN.PD.mod0,
                   `poll.dep day_night quadratic`=DvN.PD.mod)

save(poll.dep.mods,file="output/DvN_meta_regression_pollination_dependency_models Sep16.RData")

#####
#Meta-regression: Day vs. night (vs. open) pollination as a function of categorical plant traits
# Categorical traits

#get column names
cat.traits=c(colnames(es.list$traits[,4:14]))[-7]# remove non simple bloom period

#treatment levels
treatment.levels=unique(diel.all.diffs$treatment)

#create list
cat.trait.mod.list=list()

#cat.trait.mod.list.0=list()
for(h in 1:length(cat.traits)){
  k=cat.traits[h]
  for(i in treatment.levels){
    
    out.df=diel.all.diffs %>% 
      filter(treatment %in% i) %>% 
      filter(!is.na(yi))
    
    out.df[,k]=as.factor(out.df[,k])
    
    out.species=unique(out.df$phylo)
    out.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                             out.species))
    out.vcv <- vcv(out.tree, cor = T)
    
    out=rma.mv(yi,vi,mods = as.formula(paste("~exposure+",k)),
               random = list(~1|treatment_effectiveness_metric,
                             ~1|study_ID,
                             ~1|effect_ID,
                             ~1|phylo,
                             ~1|accepted_name),
               R = list(phylo = out.vcv),
               data = out.df)
    
    cat.trait.mod.list[[paste(k,i)]]=out
    
    
  }
}

#save model list
save(cat.trait.mod.list,file="output/DvN_meta_regression_categorical_trait_models Sep16.RData")

#####
#Meta-regression: Day vs. night (vs. open) pollination as a function of continuous plant traits (linear) and environmental variables (linear and quadratic)

#get column names
cont.env.traits=c("sStyle_length","sPlant_height","sDTR","sDaylength","sElevation")

#treatment levels
treatment.levels=unique(diel.all.diffs$treatment)

#create list
cont.env.trait.mod.list=list()

for(h in 1:length(cont.env.traits)){
  
  k=cont.env.traits[h]
  
  for(i in treatment.levels){
    
    out.df=diel.all.diffs %>% 
      filter(treatment %in% i) %>% 
      filter(!is.na(yi))
    
    out.species=unique(out.df$phylo)
    out.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                             out.species))
    out.vcv <- vcv(out.tree, cor = T)
    
    out=rma.mv(yi,vi,mods = as.formula(paste("~exposure+",k)),
               random = list(~1|treatment_effectiveness_metric,
                             ~1|study_ID,
                             ~1|effect_ID,
                             ~1|phylo,
                             ~1|accepted_name),
               R = list(phylo = out.vcv),
               data = out.df)
    
    cont.env.trait.mod.list[[paste(k,i,"linear")]] = out
    
    if(k %in% c("sDTR","sDaylength","sElevation")){
      
      out2=rma.mv(yi,vi,mods = as.formula(paste("~exposure+",k,"+","I(",k,"^2)")),
                  random = list(~1|treatment_effectiveness_metric,
                                ~1|study_ID,
                                ~1|effect_ID,
                                ~1|phylo,
                                ~1|accepted_name),
                  R = list(phylo = out.vcv),
                  data = out.df)
      
      cont.env.trait.mod.list[[paste(k,i,"quadratic")]] =out2
    }
    
  }
}

save(cont.env.trait.mod.list,file="output/DvN_meta_regression_continuous_trait_environment_models Sep16.RData")


#####
#Publication bias: Day vs. night (vs. open) pollination differences - multilevel Egger tests and time-lag models
#Models of sqrt of vi with all random effects
#Models of study year
####

#treatment levels
treatment.levels=unique(diel.all.diffs$treatment)

#create lists
egger.list=list()
time.lag.list=list()

for(i in treatment.levels){
  
  out.df=diel.all.diffs %>% 
    filter(treatment %in% i)
  
  out.species=unique(out.df$phylo)
  out.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                           out.species))
  out.vcv <- vcv(out.tree, cor = T)
  
  egger.full <- rma.mv(yi = yi, 
                       V = vi,
                       mods = ~ sqrt(vi),
                       random = list(~1|treatment_effectiveness_metric,
                                     ~1|study_ID,
                                     ~1|effect_ID,
                                     ~1|phylo,
                                     ~1|accepted_name),
                       R = list(phylo = out.vcv),
                       data = out.df)
  
  #time lag 
  time.lag.full <- rma.mv(yi = yi, 
                          V = vi,
                          mods = ~ year,
                          random = list(~1|treatment_effectiveness_metric,
                                        ~1|study_ID,
                                        ~1|effect_ID,
                                        ~1|phylo,
                                        ~1|accepted_name),
                          R = list(phylo = out.vcv),
                          data = out.df %>% 
                            mutate(year=year_start-min(year_start)))
  egger.list[[i]]=egger.full
  time.lag.list[[i]]=time.lag.full
  
}

save(egger.list,file="output/DvN_publication_bias_egger_models June5.RData")
save(time.lag.list,file="output/DvN_publication_bias_time_lag_models June5.RData")
