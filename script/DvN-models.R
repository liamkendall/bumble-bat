#############################
###Dave vs. Nathan###
#############################

##load effect size dataframes rdata file

load("data/effect_size_dataframes.rData")

es.list
####
library(metafor)
library(orchaRd)

#Of note, an SMD of 0 means that there is no difference between groups, 
#and an SMD that is negative means that the experimental group has a lower 
#mean score than the control group (this is when the numerator for SMD is 
#calculated as experimental minus control and the negative sign is retained)

##day night effect size list
day.night.df=es.list$day.night.df

diel.fs=day.night.df %>% 
  filter(treatment_effectiveness_metric%in%"fruit set") %>% #View
  filter(!effect_ID%in%c(47,48)) %>% 
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

#left join trait.groups to diel.fs


diel.fs$sDTR=as.numeric(scale(diel.fs$DTR))
diel.fs$sDaylength=as.numeric(scale(diel.fs$Daylength))

forest(diel.fs$yi,diel.fs$vi)

#funnel plot
funnel(diel.fs$yi,diel.fs$vi)

diel.ss=day.night.df %>% 
  filter(treatment_effectiveness_metric%in%"seed set")%>% 
  filter(!effect_ID%in%175)%>% 
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

#forest plot
forest(diel.ss$yi,diel.ss$vi)

#funnel plot
funnel(diel.ss$yi,diel.ss$vi)

########
#fruit set models
########

###set up phylo
diel.fs.species=unique(diel.fs$phylo)
setdiff(diel.fs.species,es.list$phylo$tip.label)
diel.fs.tree=drop.tip(es.list$phylo, setdiff(es.list$phylo$tip.label,
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
                               R = list(phylo = diel_fs_vcv),
                             data = diel.fs)

summary(diel_fs_mod_4)

AIC(diel_fs_mod_4,
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
summary(diel_fs_sEnv_mod_1)

AIC(diel_fs_sEnv_mod_1,
    diel_fs_mod_2)

diel_fs_sEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ sDTR+I(sDTR^2)+sDaylength+I(sDaylength^2), 
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID, 
                                                 ~ 1 | phylo), 
                                  R = list(phylo = diel_fs_vcv),
                                  data = diel.fs)
summary(diel_fs_sEnv_mod_2)

AIC(diel_fs_sEnv_mod_1,
    diel_fs_sEnv_mod_2)


#fs_DTR_bubble <- mod_results(diel_open_fs_sEnv_mod_2, mod = "sDTR", 
#                             group = "study_ID",
#                             weights = "prop"#, by = "treatment_condition"
#                             )
#
#fs_DTR_bubble_plot=bubble_plot(fs_DTR_bubble, group = "study_ID",  mod = "sDTR", xlab = "Daily Temperature Range", legend.pos = "top.left")
#
#fs_DTR_bubble_out=fs_DTR_bubble_plot+
#  facet_wrap(~condition,scales="free")
#
#ggsave(fs_DTR_bubble_out,file="figures/fs DTR bubble.jpg")#,width = 5,height=5,units="in")

#check latitude
diel.fs$abs.lat=scale(abs(diel.fs$Lat_int))

#correalted with daylength but not DTR

diel_fs_sEnv_mod_4 <- rma.mv(yi = yi, V = vi,mods = ~  abs.lat+I(abs.lat^2), 
                                    random = list( ~ 1 | study_ID, 
                                                   ~ 1 | effect_ID, 
                                                   ~ 1 | phylo), 
                                    R = list(phylo = diel_fs_vcv),
                                    data = diel.fs)
AIC(diel_fs_sEnv_mod_4,
    diel_fs_sEnv_mod_3)

summary(diel_fs_sEnv_mod_4)

####try with pollination dependency
source('script/DvN-pollination-dependency.R')
diel.fs.pd = diel.fs %>% 
  left_join(diel.pd.out) %>% 
  filter(!is.na(pd))

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
                                               ~ 1 | accepted_name), 
                                #R = list(phylo = diel_fs_pd_vcv),
                                data = diel.fs.pd)
summary(diel_fs_pd_mod_2)

AIC(diel_fs_pd_mod_1,
    diel_fs_pd_mod_2)

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

ggsave(pd_bubble_plot_out,file="pd fruit set.jpg",
       width=6,height=6,dpi=600)

diel_fs_pd_mod_3 <- rma.mv(yi = yi, V = vi,mods = ~ pd + I(pd^2) + sDTR + I(sDTR^2) + sDaylength+ I(sDaylength^2), 
                                random = list( ~ 1 | study_ID, 
                                               ~ 1 | effect_ID, 
                                               ~ 1 | phylo), 
                                R = list(phylo = diel_fs_pd_vcv),
                                data = diel.fs.pd)

summary(diel_fs_pd_mod_3)

AIC(diel_fs_pd_mod_1,
    diel_fs_pd_mod_2,
    diel_fs_pd_mod_3)

###############
#TRAIT GROUPS
#############
diel.fs$group=factor(diel.fs$group)

diel_fs_trt_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ group, 
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID, 
                                            ~ 1 | phylo), 
                             R = list(phylo = diel_fs_vcv),
                             data = diel.fs)
summary(diel_fs_trt_mod_1)

AIC(diel_fs_mod_3,
    diel_fs_trt_mod_1)

#subset trait group 5
diel.fs.trt5=diel.fs %>% 
  filter(group==1)

trt_model <- mod_results(diel_fs_trt_mod_1, 
                          group = "study_ID", 
                          mod = "group")

trt_orchard=orchard_sans_pi(trt_model, xlab = "Standardised Mean Difference")

trt_orchard

#blooming period
diel_fs_trt_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ bloom_period_simple, 
                            random = list( ~ 1 | study_ID, 
                                           ~ 1 | effect_ID, 
                                           ~ 1 | phylo), 
                            R = list(phylo = diel_fs_vcv),
                            data = diel.fs)
summary(diel_fs_trt_mod_2)


diel_fs_trt_mod_3 <- rma.mv(yi = yi, V = vi,mods = ~ flower_length_midpoint_mm, 
                            random = list( ~ 1 | study_ID, 
                                           ~ 1 | effect_ID, 
                                           ~ 1 | phylo), 
                            R = list(phylo = diel_fs_vcv),
                            data = diel.fs)
summary(diel_fs_trt_mod_3)

diel_fs_trt_mod_4 <- rma.mv(yi = yi, V = vi,mods = ~ flower_width_midpoint_mm, 
                            random = list( ~ 1 | study_ID, 
                                           ~ 1 | effect_ID, 
                                           ~ 1 | phylo), 
                            R = list(phylo = diel_fs_vcv),
                            data = diel.fs)
summary(diel_fs_trt_mod_4)

diel_fs_trt_mod_5 <- rma.mv(yi = yi, V = vi,mods = ~ flower_width_midpoint_mm, 
                            random = list( ~ 1 | study_ID, 
                                           ~ 1 | effect_ID, 
                                           ~ 1 | phylo), 
                            R = list(phylo = diel_fs_vcv),
                            data = diel.fs)
summary(diel_fs_trt_mod_4)
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
                                            ~ 1 | accepted_name),  # non-phylogenetic species effect
                             #  R = list(phylo = diel_ss_vcv),
                             data = diel.ss)

summary(diel_ss_mod_4)

AIC(diel_ss_mod_4,
    diel_ss_mod_3,
    diel_ss_mod_2,
    diel_ss_mod_1)

###no overall effect
###low effect of phylogeny
ss_model <- mod_results(diel_ss_mod_4, 
                             group = "study_ID")

ss_orchard=orchard_sans_pi(ss_model, xlab = "Standardised Mean Difference")

ss_orchard

####add co-variates
diel_ss_sEnv_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ 
                                    sDTR * 
                                    sDaylength,
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID,
                                                 ~ 1 | phylo), 
                                  R = list(phylo = diel_ss_vcv),
                                  data = diel.ss)

summary(diel_ss_sEnv_mod_1)

diel_ss_sEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ sDTR + I(sDTR^2) + 
                                    sDaylength + I(sDaylength^2),
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID,
                                                 ~ 1 | phylo), 
                                  #  R = list(phylo = diel_ss_vcv),
                                  data = diel.ss)

AIC(diel_ss_sEnv_mod_1,
    diel_ss_sEnv_mod_2)

summary(diel_ss_sEnv_mod_2)


diel_ss_sEnv_mod_3 <- rma.mv(yi = yi, V = vi,mods = ~ sDTR + I(sDTR^2), #+ 
                            #   sDaylength + I(sDaylength^2),
                             random = list( ~ 1 | study_ID, 
                                            ~ 1 | effect_ID,
                                            ~ 1 | phylo), 
                             #  R = list(phylo = diel_ss_vcv),
                             data = diel.ss)

AIC(diel_ss_sEnv_mod_1,
    diel_ss_sEnv_mod_2,
    diel_ss_sEnv_mod_3)

summary(diel_ss_sEnv_mod_3)

dtr_bubble <- mod_results(diel_ss_sEnv_mod_3, mod = "sDTR", 
                          group = "study_ID",
                          weights = "prop")

dtr_bubble_plot=bubble_plot(dtr_bubble, group = "study_ID",  mod = "sDTR", xlab = "Diurnal Temperature Range", legend.pos = "top.left")

dtr_bubble_plot_out=dtr_bubble_plot+
  geom_hline(yintercept = 0)+
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face="bold",size=10))

ggsave(dtr_bubble_plot,file="figures/dtr bubble.jpg")#,width = 5,height=5,units="in")

###TRY with pollination dependency

diel.ss.pd = diel.ss %>% 
  left_join(diel.pd.out) %>% 
  filter(!is.na(pd))

diel.ss.pd.species=diel.ss.pd$phylo
diel.ss.pd.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                          diel.ss.pd.species))
diel_ss_pd_vcv <- vcv(diel.ss.pd.tree, cor = T)

diel_ss_pd_mod_1 <- rma.mv(yi = yi, V = vi,mods = ~ pd,#+I(pd^2), 
                           random = list( ~ 1 | study_ID, 
                                          ~ 1 | effect_ID, 
                                          ~ 1 | phylo), 
                           R = list(phylo = diel_ss_pd_vcv),
                           data = diel.ss.pd)
summary(diel_ss_pd_mod_1)

diel_ss_pd_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ pd+I(pd^2), 
                           random = list( ~ 1 | study_ID, 
                                          ~ 1 | effect_ID, 
                                          ~ 1 | phylo), 
                           R = list(phylo = diel_ss_pd_vcv),
                           data = diel.ss.pd)
summary(diel_ss_pd_mod_2)

AIC(diel_ss_pd_mod_1,
    diel_ss_pd_mod_2)

pd_ss_bubble <- mod_results(diel_ss_pd_mod_2, mod = "pd", 
                         group = "study_ID",
                         weights = "prop" #, by = "treatment_condition"
)

pd_ss_bubble_plot=bubble_plot(pd_ss_bubble, group = "study_ID",  
                           mod = "pd", xlab = "Pollination dependency", 
                           legend.pos = "top.left")

pd_ss_bubble_plot_out=pd_ss_bubble_plot+
  geom_hline(yintercept = 0)+
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face="bold",size=10))

pd_ss_bubble_plot_out

#ggsave(pd_bubble_plot_out,file="pd fruit set.jpg",
#       width=6,height=6,dpi=600)

diel_ss_pd_mod_3 <- rma.mv(yi = yi, V = vi,mods = ~ pd + I(pd^2) + sDTR + I(sDTR^2),# + sDaylength+ I(sDaylength^2), 
                           random = list( ~ 1 | study_ID, 
                                          ~ 1 | effect_ID, 
                                          ~ 1 | phylo), 
                           R = list(phylo = diel_ss_pd_vcv),
                           data = diel.ss.pd)

summary(diel_ss_pd_mod_3)

AIC(diel_ss_pd_mod_1,
    diel_ss_pd_mod_2,
    diel_ss_pd_mod_3)

