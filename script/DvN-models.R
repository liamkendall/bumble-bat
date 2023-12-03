#############################
###Dave vs. Nathan###
#############################

####
library(metafor)
library(orchaRd)

#Of note, an SMD of 0 means that there is no difference between groups, 
#and an SMD that is negative means that the experimental group has a lower 
#mean score than the control group (this is when the numerator for SMD is 
#calculated as experimental minus control and the negative sign is retained)

diel.fs=day.night.df %>% 
  filter(treatment_effectiveness_metric%in%"fruit set") %>% #View
  #filter(!effect_ID%in%c(47,48)) %>% 
  filter(!is.na(yi)) %>%
  filter(!is.na(DTR)) %>% #7
  filter(!is.na(Daylength))%>% #dim
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
  mutate(label = cur_group_id()) %>% 
  ungroup()

diel.fs$sDTR=as.numeric(scale(diel.fs$DTR))
diel.fs$sDaylength=as.numeric(scale(diel.fs$Daylength))

forest(diel.fs$yi,diel.fs$vi)

##re-equate
#bad.labels.fs=table(diel.fs$label) %>% as.data.frame() %>% 
#  filter(Freq%in%1)
#
#diel.fs=diel.fs %>% 
#  filter(!label%in%bad.labels.fs$Var1)

diel.ss=day.night.df %>% 
  filter(treatment_effectiveness_metric%in%"seed set")%>% 
  #filter(!effect_ID%in%175)%>% 
  filter(!is.na(DTR)) %>% #4
  filter(!is.na(Daylength))%>% #dim
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
  mutate(label = cur_group_id()) %>% 
  ungroup()

diel.ss$sDTR=as.numeric(scale(diel.ss$DTR))
diel.ss$sDaylength=as.numeric(scale(diel.ss$Daylength))
diel.ss$lDTR=log(diel.ss$DTR)
diel.ss$lDaylength=log(diel.ss$Daylength)

forest(diel.ss$yi,diel.ss$vi)

forest(day.night.df$yi,
       day.night.df$vi)

day.night.df

########
#fruit set models
########

###set up phylo
diel.fs.species=diel.fs$phylo
diel.fs.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                            diel.fs.species))
diel_fs_vcv <- vcv(diel.fs.tree, cor = T)
###

diel_fs_mod_1 <- rma.mv(yi = yi, V = vi, #mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID), 
                             data = diel.fs)
summary(diel_fs_mod_1)

diel_fs_mod_2 <- rma.mv(yi = yi, V = vi, #mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID), 
                             data = diel.fs)
AIC(diel_fs_mod_1,
    diel_fs_mod_2)

summary(diel_fs_mod_2)

diel_fs_mod_3 <- rma.mv(yi = yi, V = vi, #mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            ~ 1 | phylo), 
                             R = list(phylo = diel_fs_vcv),
                             data = diel.fs)
summary(diel_fs_mod_3)

AIC(diel_fs_mod_1,
    diel_fs_mod_2,
    diel_fs_mod_3)

#phylo does not improve it
diel_fs_mod_4 <- rma.mv(yi = yi, V = vi, #mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            #~ 1 | phylo,
                                            ~ 1 | accepted_name),  # non-phylogenetic species effect
                             #  R = list(phylo = diel_fs_vcv),
                             data = diel.fs)

summary(diel_fs_mod_4)

AIC(diel_fs_mod_4,
    diel_fs_mod_3,
    diel_fs_mod_2,
    diel_fs_mod_1)

#phylo does not improve it
diel_fs_mod_5 <- rma.mv(yi = yi, V = vi, #mods = ~ treatment_condition,
                        random = list( ~ 1 | study_ID, 
                                       ~ 1 | effect_ID,
                                       ~ 1 | phylo,
                                       ~ 1 | accepted_name),  # non-phylogenetic species effect
                          R = list(phylo = diel_fs_vcv),
                        data = diel.fs)

summary(diel_fs_mod_5)

AIC(diel_fs_mod_5,
    diel_fs_mod_4,
    diel_fs_mod_3,
    diel_fs_mod_2,
    diel_fs_mod_1)


###no overall effect
###low effect of phylogeny

####add co-variates
ggplot(diel.fs,aes(x=sDTR,y=yi))+
  geom_point()+
  geom_smooth(formula="y~x+I(x^2)",method="lm")#+
  #geom_smooth(formula="y~x",method="lm")

ggplot(diel.fs,aes(x=sDaylength,y=yi))+
  geom_point()+
  geom_smooth(formula="y~x+I(x^2)",method="lm")#+
 # geom_smooth(formula="y~x",method="lm")

diel_fs_sEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ #treatment_condition *
                                    sDTR*
                                    sDaylength,
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID, 
                                                 ~ 1 | phylo), 
                                  R = list(phylo = diel_fs_vcv),
                                  data = diel.fs)
summary(diel_open_fs_sEnv_mod_1)

AIC(diel_fs_sEnv_mod_1,
    diel_fs_mod_2)

diel.fs$lDTR=log(diel.fs$DTR)
diel.fs$lDaylength=log(diel.fs$Daylength)

diel_fs_logEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ #treatment_condition *
                                      lDTR*
                                      lDaylength,
                                    random = list( ~ 1 | study_ID, 
                                                   ~ 1 | effect_ID, 
                                                   ~ 1 | phylo), 
                                    R = list(phylo = diel_fs_vcv),
                                    data = diel.fs)
summary(diel_fs_logEnv_mod_1)

AIC(diel_fs_sEnv_mod_1,
    diel_fs_logEnv_mod_1,
    diel_fs_mod_2)

diel_open_fs_sEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ sDTR, 
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID, 
                                                 ~ 1 | phylo), 
                                  R = list(phylo = diel_fs_vcv),
                                  data = diel.fs)
summary(diel_open_fs_sEnv_mod_2)

AIC(diel_open_fs_sEnv_mod_1,
    diel_open_fs_sEnv_mod_2)


diel_open_fs_logEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~  lDTR, 
                                    random = list( ~ 1 | study_ID, 
                                                   ~ 1 | effect_ID, 
                                                   ~ 1 | phylo), 
                                    R = list(phylo = diel_fs_vcv),
                                    data = diel.fs)
summary(diel_open_fs_logEnv_mod_2)

AIC(diel_open_fs_sEnv_mod_1,
    diel_open_fs_sEnv_mod_2,
    diel_open_fs_logEnv_mod_2)

fs_DTR_bubble <- mod_results(diel_open_fs_sEnv_mod_2, mod = "sDTR", 
                             group = "study_ID",
                             weights = "prop"#, by = "treatment_condition"
                             )

fs_DTR_bubble_plot=bubble_plot(fs_DTR_bubble, group = "study_ID",  mod = "sDTR", xlab = "Daily Temperature Range", legend.pos = "top.left")

fs_DTR_bubble_out=fs_DTR_bubble_plot+
  facet_wrap(~condition,scales="free")

ggsave(fs_DTR_bubble_out,file="figures/fs DTR bubble.jpg")#,width = 5,height=5,units="in")


diel.fs$abs.lat=scale(abs(diel.fs$Lat_int))
diel.fs$hemisphere=ifelse(diel.fs$Lat_int<0,"south","north")
diel.fs$sYear=scale(as.numeric(diel.fs$year_start))

#correalted with daylength but not DTR

diel_open_fs_sEnv_mod_3 <- rma.mv(yi = yi, V = vi,mods = ~  (sDTR+I(sDTR^2))*
                                    (abs.lat+I(abs.lat^2)), 
                                    random = list( ~ 1 | study_ID, 
                                                   ~ 1 | effect_ID, 
                                                   ~ 1 | phylo), 
                                    R = list(phylo = diel_fs_vcv),
                                    data = diel.fs)
AIC(diel_open_fs_sEnv_mod_3,
    diel_open_fs_sEnv_mod_2)

summary(diel_open_fs_sEnv_mod_3)

####try with pollination dependency
diel.fs.pd = diel.fs %>% 
  left_join(diel.fs.pd.imp.out) %>% 
  filter(!is.na(iPd))

diel.fs.pd.species=diel.fs.pd$phylo
diel.fs.pd.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                               diel.fs.pd.species))
diel_fs_pd_vcv <- vcv(diel.fs.pd.tree, cor = T)

diel_fs_pd_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ pd,#+I(pd^2), 
                                random = list( ~ 1 | study_ID, 
                                               ~ 1 | effect_ID, 
                                               ~ 1 | phylo), 
                                R = list(phylo = diel_fs_pd_vcv),
                                data = diel.fs.pd)
summary(diel_fs_pd_mod_1)



diel_fs_pd_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ pd+I(pd^2), 
                                random = list( ~ 1 | study_ID, 
                                               ~ 1 | effect_ID, 
                                               ~ 1 | phylo), 
                                R = list(phylo = diel_fs_pd_vcv),
                                data = diel.fs.pd)
summary(diel_fs_pd_mod_2)

AIC(diel_open_fs_pd_mod_1,
    diel_open_fs_pd_mod_2)

pd_bubble <- mod_results(diel_fs_pd_mod_2, mod = "pd", 
                         group = "study_ID",
                         weights = "prop" #, by = "treatment_condition"
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

diel_fs_pd_mod_3 <- rma.mv(yi = yi, V = vi,mods = ~ (pd + I(pd^2)) * (sDTR  * sDaylength), 
                                random = list( ~ 1 | study_ID, 
                                               ~ 1 | effect_ID, 
                                               ~ 1 | phylo), 
                                R = list(phylo = diel_fs_pd_vcv),
                                data = diel.fs.pd)

summary(diel_fs_pd_mod_3)

AIC(diel_open_fs_pd_mod_1,
    diel_open_fs_pd_mod_2,
    diel_open_fs_pd_mod_3)


###############
#seed set
#########

###set up phylo
diel.ss.species=diel.ss$phylo
diel.ss.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                            diel.ss.species))
diel_ss_vcv <- vcv(diel.ss.tree, cor = T)
###

diel_ss_mod_1 <- rma.mv(yi = yi, V = vi, #mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID), 
                             data = diel.ss)
summary(diel_ss_mod_1)

diel_ss_mod_2 <- rma.mv(yi = yi, V = vi, #mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID), 
                             data = diel.ss)
AIC(diel_ss_mod_1,
    diel_ss_mod_2)

summary(diel_ss_mod_2)

diel_ss_mod_3 <- rma.mv(yi = yi, V = vi,# mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            ~ 1 | phylo), 
                             R = list(phylo = diel_ss_vcv),
                             data = diel.ss)
summary(diel_ss_mod_3)

AIC(diel_ss_mod_1,
    diel_ss_mod_2,
    diel_ss_mod_3)

#phylo does not improve it
diel_ss_mod_4 <- rma.mv(yi = yi, V = vi, #mods = ~ treatment_condition,
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            #~ 1 | phylo,
                                            ~ 1 | family),  # non-phylogenetic species effect
                             #  R = list(phylo = diel_ss_vcv),
                             data = diel.ss)

summary(diel_ss_mod_4)

AIC(diel_ss_mod_4,
    #diel_ss_mod_3,
    diel_ss_mod_2,
    diel_ss_mod_1)

###no overall effect
###low effect of phylogeny
ss_model <- mod_results(diel_ss_mod_4, 
                             group = "study_ID")

ss_orchard=orchard_sans_pi(ss_model, xlab = "Standardised Mean Difference")

ss_orchard

####add co-variates
diel_open_ss_sEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition *  
                                    sDTR * 
                                    sDaylength,
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID,
                                                 ~ 1 | phylo), 
                                  R = list(phylo = diel_ss_vcv),
                                  data = diel.ss)

summary(diel_open_ss_sEnv_mod_1)

diel_open_ss_logEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_condition *  
                                      lDTR * 
                                      lDaylength,
                                    random = list( ~ 1 | study_ID, 
                                                   ~ 1 | effect_ID,
                                                   ~ 1 | phylo), 
                                    R = list(phylo = diel_ss_vcv),
                                    data = diel.ss)

summary(diel_open_ss_logEnv_mod_1)

AIC(diel_open_ss_sEnv_mod_1,
    diel_open_ss_logEnv_mod_1,
    diel_ss_mod_2)

diel_open_ss_sEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~treatment_condition* sDTR,
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID), 
                                  #  R = list(phylo = diel_ss_vcv),
                                  data = diel.ss)

AIC(diel_open_ss_sEnv_mod_1,
    diel_open_ss_sEnv_mod_2,
    diel_ss_mod_2)

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



open_model <- mod_results(diel_fs_mod_2, 
                          group = "study_ID", 
                          mod = "treatment_condition")

open_orchard=orchard_plot(open_model, xlab = "Standardised Mean Difference")

open_orchard

