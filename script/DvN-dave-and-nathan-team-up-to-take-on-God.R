#############################
###Dave and Nathan vs. God###
#############################

####
library(metafor)
library(orchaRd)

#Of note, an SMD of 0 means that there is no difference between groups, 
#and an SMD that is negative means that the experimental group has a lower 
#mean score than the control group (this is when the numerator for SMD is 
#calculated as experimental minus control and the negative sign is retained)

diel.open.fs=diel.open.es.out %>% 
  filter(treatment_effectiveness_metric%in%"fruit set") %>% #View
  filter(!effect_ID%in%c(91,93)) %>% 
  filter(!is.na(yi)) %>%
  filter(!is.na(DTR)) %>% #7
  filter(!is.na(Daylength))

diel.open.fs$sDTR=scale(diel.open.fs$DTR)
diel.open.fs$sDaylength=scale(diel.open.fs$Daylength)

forest(diel.open.fs$yi,diel.open.fs$vi)

##re-equate
bad.labels.fs=table(diel.open.fs$label) %>% as.data.frame() %>% 
  filter(Freq%in%1)

diel.open.fs=diel.open.fs %>% 
  filter(!label%in%bad.labels.fs$Var1)

diel.open.ss=diel.open.es.out %>% 
  filter(treatment_effectiveness_metric%in%"seed set")%>% 
  #filter(VI.dn<2)%>% 
  filter(!is.na(DTR)) %>% #4
  filter(!is.na(Daylength))

diel.open.ss$sDTR=scale(diel.open.ss$DTR)
diel.open.ss$sDaylength=scale(diel.open.ss$Daylength)
diel.open.ss$lDTR=log(diel.open.ss$DTR)
diel.open.ss$lDaylength=log(diel.open.ss$Daylength)

forest(diel.open.ss$yi,diel.open.ss$vi)

diel.open.fs.day=diel.open.fs %>% 
  filter(treatment_condition%in%"day pollination")

diel.open.fs.night=diel.open.fs %>% 
  filter(treatment_condition%in%"night pollination")

day.check=table(diel.open.fs.day$study_ID) %>% as.data.frame()
night.check=table(diel.open.fs.night$study_ID) %>% as.data.frame()

dn.check=day.check %>% 
  left_join(night.check,by="Var1") %>% 
  mutate(diff=Freq.x-Freq.y)

table(dn.check$diff)


########
#fruit set models
########

###set up phylo
diel.open.fs.species=diel.open.fs$phylo
diel.open.fs.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                            diel.open.fs.species))
diel_open_fs_vcv <- vcv(diel.open.fs.tree, cor = T)
###

diel_open_fs_mod_1 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID), 
                             data = diel.open.fs)
summary(diel_open_fs_mod_1)

diel_open_fs_mod_2 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID), 
                             data = diel.open.fs)
AIC(diel_open_fs_mod_1,
    diel_open_fs_mod_2)

summary(diel_open_fs_mod_2)

diel_open_fs_mod_3 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            ~ 1 | phylo), 
                             R = list(phylo = diel_open_fs_vcv),
                             data = diel.open.fs)
summary(diel_open_fs_mod_3)

AIC(diel_open_fs_mod_1,
    diel_open_fs_mod_2,
    diel_open_fs_mod_3)

#phylo does not improve it
diel_open_fs_mod_4 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            #~ 1 | phylo,
                                            ~ 1 | accepted_name),  # non-phylogenetic species effect
                           #  R = list(phylo = diel_open_fs_vcv),
                             data = diel.open.fs)

summary(diel_open_fs_mod_4)

AIC(diel_open_fs_mod_4,
    diel_open_fs_mod_3,
    diel_open_fs_mod_2,
    diel_open_fs_mod_1)

###no overall effect
###low effect of phylogeny

####add co-variates
ggplot(diel.open.fs,aes(x=sDTR,y=yi,col=treatment_condition))+
  geom_point()+
  geom_smooth(formula="y~x+I(x^2)",method="lm")+
  geom_smooth(formula="y~x",method="lm")

ggplot(diel.open.fs,aes(x=sDaylength,y=yi,col=treatment_condition))+
  geom_point()+
  geom_smooth(formula="y~x+I(x^2)",method="lm")+
  geom_smooth(formula="y~x",method="lm")

diel_open_fs_sEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition *
                                    sDTR*
                                    sDaylength,
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID, 
                                                 ~ 1 | phylo), 
                                  R = list(phylo = diel_open_fs_vcv),
                                  data = diel.open.fs)
summary(diel_open_fs_sEnv_mod_1)

AIC(diel_open_fs_sEnv_mod_1,
    diel_open_fs_mod_2)

diel.open.fs$lDTR=log(diel.open.fs$DTR)

diel.open.fs$lDaylength=log(diel.open.fs$Daylength)

diel_open_fs_logEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition *
                                    lDTR*
                                    lDaylength,
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID, 
                                                 ~ 1 | phylo), 
                                  R = list(phylo = diel_open_fs_vcv),
                                  data = diel.open.fs)
summary(diel_open_fs_logEnv_mod_1)

AIC(diel_open_fs_sEnv_mod_1,
    diel_open_fs_logEnv_mod_1,
    diel_open_fs_mod_2)

diel_open_fs_sEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition * sDTR, 
                                 random = list( ~ 1 | study_ID, 
                                                ~ 1 | effect_ID, 
                                                ~ 1 | phylo), 
                                   R = list(phylo = diel_open_fs_vcv),
                                 data = diel.open.fs)
summary(diel_open_fs_sEnv_mod_2)

AIC(diel_open_fs_sEnv_mod_1,
    diel_open_fs_sEnv_mod_2)


diel_open_fs_logEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition * lDTR, 
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID, 
                                                 ~ 1 | phylo), 
                                  R = list(phylo = diel_open_fs_vcv),
                                  data = diel.open.fs)
summary(diel_open_fs_logEnv_mod_2)

AIC(diel_open_fs_sEnv_mod_1,
    diel_open_fs_sEnv_mod_2,
    diel_open_fs_logEnv_mod_2)

fs_DTR_bubble <- mod_results(diel_open_fs_sEnv_mod_2, mod = "sDTR", 
                          group = "study_ID",
                          weights = "prop", by = "treatment_condition")

fs_DTR_bubble_plot=bubble_plot(fs_DTR_bubble, group = "study_ID",  mod = "sDTR", xlab = "Daily Temperature Range", legend.pos = "top.left")

fs_DTR_bubble_out=fs_DTR_bubble_plot+
  facet_wrap(~condition,scales="free")

ggsave(fs_DTR_bubble_out,file="figures/fs DTR bubble.jpg")#,width = 5,height=5,units="in")

####try with pollination dependency
diel.open.fs.pd = diel.open.fs %>% 
  left_join(diel.fs.pd.imp.out) %>% 
  filter(!is.na(iPd))

diel.open.fs.pd.species=diel.open.fs.pd$phylo
diel.open.fs.pd.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                               diel.open.fs.pd.species))
diel_open_fs_pd_vcv <- vcv(diel.open.fs.pd.tree, cor = T)

diel_open_fs_pd_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition *pd, 
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID, 
                                                 ~ 1 | phylo), 
                                    R = list(phylo = diel_open_fs_pd_vcv),
                                  data = diel.open.fs.pd)
summary(diel_open_fs_pd_mod_1)

diel_open_fs_pd_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition + pd, 
                              random = list( ~ 1 | study_ID, 
                                             ~ 1 | effect_ID, 
                                             ~ 1 | phylo), 
                                R = list(phylo = diel_open_fs_pd_vcv),
                              data = diel.open.fs.pd)
summary(diel_open_fs_pd_mod_2)

AIC(diel_open_fs_pd_mod_1,
    diel_open_fs_pd_mod_2)

pd_bubble <- mod_results(diel_open_fs_pd_mod_2, mod = "pd", 
                          group = "study_ID",
                          weights = "prop", by = "treatment_condition"
                         )

pd_bubble_plot=bubble_plot(pd_bubble, group = "study_ID",  
                           mod = "pd", xlab = "Pollination dependency", 
                           legend.pos = "top.left")

pd_bubble_plot_out=pd_bubble_plot+
  geom_hline(yintercept = 0)+
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face="bold",size=10))
pd_bubble_plot_out

diel_open_fs_pd_mod_3 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition * sDTR + pd * sDTR + pd * sDaylength, 
                                random = list( ~ 1 | study_ID, 
                                               ~ 1 | effect_ID, 
                                               ~ 1 | phylo), 
                                R = list(phylo = diel_open_fs_pd_vcv),
                                data = diel.open.fs.pd)
summary(diel_open_fs_pd_mod_3)

AIC(diel_open_fs_pd_mod_1,
    diel_open_fs_pd_mod_2,
    diel_open_fs_pd_mod_3)


###############
#seed set
#########

###set up phylo
diel.open.ss.species=diel.open.ss$phylo
diel.open.ss.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                            diel.open.ss.species))
diel_open_ss_vcv <- vcv(diel.open.ss.tree, cor = T)
###

diel_open_ss_mod_1 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID), 
                             data = diel.open.ss)
summary(diel_open_ss_mod_1)

diel_open_ss_mod_2 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID), 
                             data = diel.open.ss)
AIC(diel_open_ss_mod_1,
    diel_open_ss_mod_2)

summary(diel_open_ss_mod_2)

diel_open_ss_mod_3 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            ~ 1 | phylo), 
                             R = list(phylo = diel_open_ss_vcv),
                             data = diel.open.ss)
summary(diel_open_ss_mod_3)

AIC(diel_open_ss_mod_1,
    diel_open_ss_mod_2,
    diel_open_ss_mod_3)

#phylo does not improve it
diel_open_ss_mod_4 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            #~ 1 | phylo,
                                            ~ 1 | label),  # non-phylogenetic species effect
                             #  R = list(phylo = diel_open_ss_vcv),
                             data = diel.open.ss)

summary(diel_open_ss_mod_4)

AIC(diel_open_ss_mod_4,
    diel_open_ss_mod_3,
    diel_open_ss_mod_2,
    diel_open_ss_mod_1)

###no overall effect
###low effect of phylogeny
open_ss_model <- mod_results(diel_open_ss_mod_2, 
                          group = "study_ID", 
                          mod = "treatment_condition")

open_ss_orchard=orchard_plot(open_ss_model, xlab = "Standardised Mean Difference")

open_ss_orchard

####add co-variates
diel_open_ss_sEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition *  
                                   sDTR * 
                                   sDaylength,
                                 random = list( ~ 1 | study_ID, 
                                                ~ 1 | effect_ID,
                                                ~ 1 | phylo), 
                                 R = list(phylo = diel_open_ss_vcv),
                                 data = diel.open.ss)

summary(diel_open_ss_sEnv_mod_1)

diel_open_ss_logEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition *  
  lDTR * 
    lDaylength,
  random = list( ~ 1 | study_ID, 
                 ~ 1 | effect_ID,
                 ~ 1 | phylo), 
  R = list(phylo = diel_open_ss_vcv),
  data = diel.open.ss)

summary(diel_open_ss_logEnv_mod_1)

AIC(diel_open_ss_sEnv_mod_1,
    diel_open_ss_logEnv_mod_1,
    diel_open_ss_mod_2)

diel_open_ss_sEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~treatment_condition* sDTR,
  random = list( ~ 1 | study_ID, 
                 ~ 1 | effect_ID), 
  #  R = list(phylo = diel_open_ss_vcv),
  data = diel.open.ss)

AIC(diel_open_ss_sEnv_mod_1,
    diel_open_ss_sEnv_mod_2,
    diel_open_ss_mod_2)

dtr_bubble <- mod_results(diel_open_ss_sEnv_mod_2, mod = "sDTR", 
                                   group = "study_ID",
                                   weights = "prop", by = "treatment_condition")

dtr_bubble_plot=bubble_plot(dtr_bubble, group = "study_ID",  mod = "sDTR", xlab = "Diurnal Temperature Range", legend.pos = "top.left")

dtr_bubble_plot+
  facet_wrap(~moderator)
dtr_bubble_plot_out=dtr_bubble_plot+
  geom_hline(yintercept = 0)+
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face="bold",size=10))

ggsave(dtr_bubble_plot,file="figures/dtr bubble.jpg")#,width = 5,height=5,units="in")



open_model <- mod_results(diel_open_fs_mod_2, 
                          group = "study_ID", 
                          mod = "treatment_condition")

open_orchard=orchard_plot(open_model, xlab = "Standardised Mean Difference")

open_orchard

