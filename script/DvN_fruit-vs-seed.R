
day.night.df2=day.night.df %>% 
  filter(!is.na(Daylength)) %>% 
  filter(!is.na(DTR))

###set up phylo
day.night.df2$phylo=gsub(" ","_",day.night.df2$accepted_name)
diel.both.species=day.night.df2$phylo
diel.both.tree=drop.tip(dn.tree, setdiff(dn.tree$tip.label,
                                       diel.both.species))
diel_both_vcv <- vcv(diel.both.tree, cor = T)
###

diel_both_mod_1 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_effectiveness_metric,
                        random = list( ~ 1 | study_ID), 
                        data = day.night.df2)
summary(diel_both_mod_1)

diel_both_mod_2 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_effectiveness_metric,
                        random = list( ~ 1 | study_ID, 
                                       ~ 1 | effect_ID), 
                        data = day.night.df2)
AIC(diel_both_mod_1,
    diel_both_mod_2)

summary(diel_both_mod_2)

diel_both_mod_3 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_effectiveness_metric,
                        random = list( ~ 1 | study_ID, 
                                       ~ 1 | effect_ID,
                                       ~ 1 | phylo), 
                        R = list(phylo = diel_both_vcv),
                        data = day.night.df2)
summary(diel_both_mod_3)

AIC(diel_both_mod_1,
    diel_both_mod_2,
    diel_both_mod_3)

#phylo does not improve it
diel_both_mod_4 <- rma.mv(yi = yi, V = vi, mods = ~ treatment_effectiveness_metric,
                        random = list( ~ 1 | study_ID, 
                                       ~ 1 | effect_ID,
                                       #~ 1 | phylo,
                                       ~ 1 | accepted_name),  # non-phylogenetic species effect
                      #    R = list(phylo = diel_both_vcv),
                        data = day.night.df2 %>% 
                        filter(!yi>5))

summary(diel_both_mod_4)

AIC(diel_both_mod_4,
    diel_both_mod_3,
    diel_both_mod_2,
    diel_both_mod_1)

both_model <- mod_results(diel_both_mod_4, 
                             group = "study_ID", 
                             mod = "treatment_effectiveness_metric")

both_orchard=orchard_plot(both_model, xlab = "Standardised Mean Difference")

both_orchard

day.night.df2$sDTR=scale(day.night.df2$DTR)
day.night.df2$sDaylength=scale(day.night.df2$Daylength)

####add co-variates
diel_both_sEnv_mod_1 <- rma.mv(yi = yi, V = vi,
                               mods = ~ treatment_effectiveness_metric *
                                 sDTR*
                                 sDaylength,
                               random = list( ~ 1 | study_ID, 
                                              ~ 1 | effect_ID, 
                            #                  ~ 1 | phylo,
                                              ~ 1 | accepted_name), 
                            #   R = list(phylo = diel_both_vcv),
                               data = day.night.df2)

summary(diel_both_sEnv_mod_1)

AIC(diel_both_sEnv_mod_1,
    diel_both_mod_4)

diel_both_sEnv_mod_2 <- rma.mv(yi = yi, V = vi,mods = ~ treatment_effectiveness_metric*sDaylength, 
                                  random = list( ~ 1 | study_ID, 
                                                 ~ 1 | effect_ID, 
                                                 ~ 1 | accepted_name), 
                                  #R = list(phylo = diel_fs_vcv),
                                  data = day.night.df2)
summary(diel_both_sEnv_mod_2)

AIC(diel_both_sEnv_mod_1,
    diel_both_sEnv_mod_2)


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

