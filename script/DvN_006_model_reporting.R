#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "20/11/2024"

###script: 006 - Model reporting

####libraries (CRAN)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(phytools)
library(metafor)
library(ggplot2)
library(cowplot)
library(tidytree)

#non-cran libraries
#to install ggtree
#devtools::install_bioc("ggtree") v.3.12.0
#to install orchaRd
#devtools::install_github("daniel1noble/orchaRd") v.2.0

#load non-cran libraries
library(ggtree)
library(orchaRd)

#additional functions
source('script/DvN_additional_functions.R')

#####GGPLOT theme
plot.theme=theme(
  legend.position = "bottom",
  axis.text = element_text(size = 10),
  text=element_text(family="Times New Roman"),
  axis.title = element_text(size = 12,face="bold"),
  legend.title = element_text(size = 12,face="bold"),
  legend.text = element_text(size = 10),
  plot.title = element_text(size = 12,face="bold", hjust = 0))


#load models
load("output/DvN_meta_analytical_models final.RData")
load("output/DvN_pub_bias_meta_analytical_models final.RData")
load("output/DvN_meta_regression_categorical_trait_models final.RData")
load("output/DvN_meta_regression_metric_models final.RData")
load("output/DvN_meta_regression_continuous_trait_environment_models final.RData")
load("output/DvN_meta_regression_pollination_dependency_models final.RData")
load("output/DvN_pollination_dependency_meta_model final.RData")
load("output/DvN_publication_bias_egger_models final.RData")
load("output/DvN_publication_bias_time_lag_models final.RData")

#load analysis dataframes
load(file="output/DvN_analysis_dataframe final.RData")
load(file="output/DvN_polldep_analysis_dataframe final.RData")

###load effect size and predictor variable dataframe (DvN_001)
load("output/DvN_effect_size_dataframe final.rData")

###load phylogeny
load(file="output/DvN_plant_phylogeny final.RData")

#####
#####Calculate heterogeneity (I^2) for meta-analytic models of each diel pollination difference

## # estimating I2 as measure of heterogeneity
#day vs night (main text)
i2_ma.dn <-i2_ml(overall.list[[1]])
#open vs night (supplementary)
i2_ma.on <- i2_ml(overall.list[[2]])
#open vs day (supplementary)
i2_ma.od <- i2_ml(overall.list[[3]])

#####FIGURE A1-4. Pollination dependency
poll.dep.mod.table = mod_results(pd.full,group="study_ID")

pd.mod.plot=orchard_sans_pi(poll.dep.mod.table,
                            group="study_ID",
                            legend.pos = "none",
                            k.pos="right",
                            xlab="SMD")+
  scale_color_manual(values=c("white","white","white","white","white"))+
  scale_fill_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  plot.theme+
  theme(aspect.ratio = 0.25,
        #axis.text.y = element_text(hjust=1),
        legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(y="Pollination dependency (SMD)")

pd.mod.plot

ggsave(pd.mod.plot,file="figures/Figure A1-4 final.jpg",
       dpi=600,width=5,height=3,units="in")

#####Figure 2 (Figures A2.1 & A2.2)
###Phylogeny + model-estimates (and phylogenetic signal) of species-level diel pollination differences
###
######

trts=unique(diel.all.diffs$treatment)
#
trt.out=trts[1]

###phylo signal and % significant species
phylo.sig.list=list()
diff.zero.list=list()
for (i in 1:length(trts)){
  
  trt.out=trts[i]
  
  #prune tree
  dn.phy.out.species=overall.list[[trt.out]]$data$phylo
  dn.phy.out.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                                  dn.phy.out.species))
  
  #get random intercepts for species
  phylo_re<-data.frame(ranef(overall.list[[trt.out]])$phylo)
  phylo_re$phylo<-rownames(phylo_re)
  colnames(phylo_re)<-c("intrcpt2","se2","pi.lb2","pi.ub2","phylo") #to avoid duplicate colnames later
  
  spp_re<-data.frame(ranef(overall.list[[trt.out]])$accepted_name)
  spp_re$accepted_name<-rownames(spp_re)
  spp_re$phylo<-gsub(" ","_",rownames(spp_re))
  
  #get predictions
  
  if(trt.out %in%"day_night"){
  preds<-predict.rma(overall.list[[trt.out]],newmods=c(1),addx=TRUE)
  }
  if(!trt.out %in%"day_night"){
    preds<-predict.rma(overall.list[[trt.out]],newmods=c(0.5),addx=TRUE)
  }
  #combine
  plotdat<-join(phylo_re,spp_re,by="phylo")
  plotdat$pred<-preds$pred
  plotdat$ci.ub<-preds$ci.ub
  plotdat$ci.lb<-preds$ci.lb
  plotdat$fe_se<-preds$se
  
  #add together intercept and speceis effects
  plotdat$blup<-plotdat$pred+plotdat$intrcpt+plotdat$intrcpt2
  
  #calculate SE
  plotdat$blup.se<-sqrt((plotdat$fe_se^2)+(plotdat$se^2)+(plotdat$se2^2))
  
  #calculate CI
  plotdat$ci.lb<-plotdat$blup-1.96*plotdat$blup.se
  plotdat$ci.ub<-plotdat$blup+1.96*plotdat$blup.se
  
  #overlap zero
  plotdat$diff.zero=ifelse(plotdat$ci.lb>0 | plotdat$ci.ub<0,0,1)
  
  dn.phy.out=data.frame(row.names = plotdat$phylo,mean=plotdat$blup)
  
  svl <- as.matrix(dn.phy.out)[,1]
  
  phylo.sig.list[[trt.out]]=phylosig(dn.phy.out.tree,svl,"lambda",test=T)
  diff.zero.list[[trt.out]]=plotdat
}

#tabulate lambda
lambda.df=data.frame(treatment=trts,
                     lambda=c(phylo.sig.list$day_night$lambda,phylo.sig.list$open_night$lambda,phylo.sig.list$open_day$lambda),
                     P=c(phylo.sig.list$day_night$P,phylo.sig.list$open_night$P,phylo.sig.list$open_day$P))
lambda.df$lambda=round(lambda.df$lambda,3)
lambda.df$P=round(lambda.df$P,3) %>% as.character
lambda.df$P[3]="<0.001"

#plot phylogeny and species estimates
phylo.gg.list=list()
i=1
for (i in 1:length(trts)){
  
  trt.out=trts[i]
  
  dn.phy.out.species=overall.list[[trt.out]]$data$phylo
  dn.phy.out.tree=drop.tip(diel.tree.out, setdiff(diel.tree.out$tip.label,
                                                  dn.phy.out.species))
  
  phylo_re<-data.frame(ranef(overall.list[[trt.out]])$phylo)
  phylo_re$phylo<-rownames(phylo_re)
  colnames(phylo_re)<-c("intrcpt2","se2","pi.lb2","pi.ub2","phylo") #to avoid duplicate colnames later
  
  spp_re<-data.frame(ranef(overall.list[[trt.out]])$accepted_name)
  spp_re$accepted_name<-rownames(spp_re)
  spp_re$phylo<-gsub(" ","_",rownames(spp_re))
  
  plotdat<-join(phylo_re,spp_re,by="phylo")
  
  if(trt.out %in%"day_night"){
    preds<-predict.rma(overall.list[[trt.out]],newmods=c(1),addx=TRUE)
  }
  if(!trt.out %in%"day_night"){
    preds<-predict.rma(overall.list[[trt.out]],newmods=c(0.5),addx=TRUE)
  }
  plotdat$pred<-preds$pred
  plotdat$ci.ub<-preds$ci.ub
  plotdat$ci.lb<-preds$ci.lb
  plotdat$fe_se<-preds$se
  
  plotdat$blup<-plotdat$pred+plotdat$intrcpt+plotdat$intrcpt2
  plotdat$blup.se<-sqrt((plotdat$fe_se^2)+
                          (plotdat$se^2)+(plotdat$se2^2))
  
  plotdat$ci.lb<-plotdat$blup-1.96*plotdat$blup.se
  plotdat$ci.ub<-plotdat$blup+1.96*plotdat$blup.se
  plotdat$diff.zero=ifelse(plotdat$ci.lb>0 | plotdat$ci.ub<0,0,1)
  
  dn.phy.out=data.frame(row.names = plotdat$phylo,mean=plotdat$blup)
  
  svl <- as.matrix(dn.phy.out)[,1]
  
  lambda.out=lambda.df %>% filter(treatment%in%trt.out)
  
  #for colouring branches by SMD (not used)
  fit <- phytools::fastAnc(dn.phy.out.tree, svl, vars=TRUE, CI=TRUE)
  
  td <- data.frame(node = nodeid(dn.phy.out.tree, names(svl)),
                   trait = svl)
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  
  
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(dn.phy.out.tree, d, by = 'node')
  tree@phylo$tip.label=gsub("_"," ",tree@phylo$tip.label)
  tree@phylo$tip.label=ifelse(tree@phylo$tip.label %in% "Ipomoea aff-marcellia",
                              "Ipomoea aff. marcellia",tree@phylo$tip.label)
  
  title=paste0(str_to_sentence(word(trt.out,1,sep="_"))," vs. ",word(trt.out,2,sep="_")," pollination ","(λ = ",lambda.out$lambda,", P = ",lambda.out$P,")")
  
  ###plot
  
  p1 <- ggtree(tree,# aes(color=trait), 
               ladderize = FALSE, continuous = 'colour', 
               size=0.5) +
    geom_tiplab(geom="text",
                inherit.aes=F,size=0,family="Times",fontface=2.5,
                align=T,alpha=0,
                col="black")+
    xlim(0,170) + 
    theme(legend.position =c(.05, .1),
          aspect.ratio=3,
          legend.background = element_blank(),
          text = element_text(family="Times"),
          title = element_text(face="bold"))#+
  # ggtitle(title)
  
  #ADD FAMILY LABELS
  
  pX=p1
  pX$data=pX$data %>% left_join(diel.all.diffs %>% distinct(accepted_name,plant_family),
                                by=c("label"="accepted_name")) 
  
  plant_families=unique(pX$data$plant_family) %>% droplevels() %>% as.character() %>% na.omit()
  
  plant.family.node=list()
  
  for(i in plant_families){
    out.pX=pX %>% as.treedata %>% as_tibble %>% filter(plant_family%in%i)
    
    if(nrow(out.pX)==1){
      plant.family.node[[i]]=out.pX$node
      
    }
    if(nrow(out.pX)>1 & !i=="NA"){
      pX.out=pX
      fam.x=pX.out$data %>% filter(plant_family%in%i)
      
      plant.family.node[[i]]=MRCA(tree, fam.x$node)[1]
    }
  }
  
  for(i in plant_families){
    
    p1=p1 + geom_cladelab(node=plant.family.node[[i]], label=i, angle=0, 
                          fontsize=2, offset=0, vjust=.5,family="Times") 
  }
  
  crop.ids=diel.all.diffs %>% 
    distinct(accepted_name,plant_crop) #%>% 
  #mutate(crop_out=ifelse(accepted_name%in%"Selenicereus costaricensis","y",crop_out))
  
  crop.ids$crop=revalue(crop.ids$plant_crop,
                        c("y"="Crop","n"="Wild"))
  
  dot.plot.c=plotdat %>%  
    left_join(crop.ids) %>% 
    mutate(accepted_name=ifelse(accepted_name %in% "Ipomoea aff-marcellia",
                                "Ipomoea aff. marcellia",accepted_name)) %>% 
    mutate(label=factor(accepted_name,
                        levels=rev(get_taxa_name(p1)))) %>% 
    ggplot(aes(y=label,x=blup,fill=crop,alpha=as.factor(diff.zero)))+
    annotate("rect",xmin=preds$ci.lb,
             xmax=preds$ci.ub,
             ymin=0,ymax=Inf,alpha=0.5,
             fill="lightgrey")+
    geom_vline(xintercept = preds$pred,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="solid",col="grey",size=0.5)+
    geom_errorbarh(aes(xmin=blup-(blup.se*1.96),
                       xmax=blup+(blup.se*1.96)),
                   height=0.5,linewidth=0.35,show.legend = F)+
    scale_alpha_manual(values=c(1,0.5))+
    geom_point(pch=21,col="black",show.legend = F,alpha=1)+
    theme_bw()+
    scale_fill_manual(values=c("white","darkgrey"))+
    theme(aspect.ratio = 3.25,
          text=element_text(family="Times New Roman"),
          axis.text.y = element_text(size=4,face="italic"),
          axis.text.x = element_text(size=6),#,face="italic"),
          axis.title=element_text(face="bold"),
          plot.margin = margin(0,0,0,0),
          axis.ticks.length = unit(0, "pt"))+
    xlab("Diel pollination difference (SMD)")+
    ylab(NULL)+
    scale_y_discrete(position = "right")
  
  if(trt.out=="day_night"){
    phylo.gg.list[[trt.out]]=cbind(ggplotGrob(p1),ggplotGrob(dot.plot.c+ 
                                                               draw_image("images/day.png",
                                                                          y = 135.5,x=-3.5,scale=4)+
                                                               draw_image("images/night.png",
                                                                          y = 135.5,x=2.5,scale=4)))
  }else{
    
    phylo.gg.list[[trt.out]]=cbind(ggplotGrob(p1),ggplotGrob(dot.plot.c))
  }
}
#save plots
ggsave(phylo.gg.list[["day_night"]],
       file="figures/Figure 2 final.jpg",
       dpi=600,height=10,width=8.3,units="in")
ggsave(phylo.gg.list[["open_day"]],
       file="figures/Figure A2-1 final.jpg",
       dpi=600,height=10,width=8.3,units="in")
ggsave(phylo.gg.list[["open_night"]],
       file="figures/Figure A2-2 final.jpg",
       dpi=600,height=10,width=8.3,units="in")


######Figure 3 (Figures A2.3 & 2.4)
#by effectiveness metric
#####

#Effectiveness metric models, generate dataframe and orchard plot using orchaRd and simplified function
eff.mod.tables=lapply(eff.list,function (x) mod_results(x,group="study_ID",mod = "treatment_effectiveness_metric",
                                                        at = list(exposure = c(1))))

eff.mod.tables=list(day_night=mod_results(eff.list[["day_night"]],
                                          group="study_ID",
                                          mod = "treatment_effectiveness_metric",
                                          at = list(exposure = c(1))),
                    open_day=mod_results(eff.list[["open_day"]],
                                          group="study_ID",
                                          mod = "treatment_effectiveness_metric",
                                          at = list(exposure = c(0.5))),
                    open_night=mod_results(eff.list[["open_night"]],
                                         group="study_ID",
                                         mod = "treatment_effectiveness_metric",
                                         at = list(exposure = c(0.5))))

eff.mod.plots=lapply(eff.mod.tables,function (x) orchard_sans_pi(x,group="study_ID",legend.pos = "none",
                                                                 k.pos="right",mod = "treatment_effectiveness_metric",xlab="SMD",
                                                                 tree.order=(c("Seed set","Seed mass","Pollen deposition","Fruit set",
                                                                               "Fruit mass"))))

dn.eff.plot=eff.mod.plots[[1]]+
  scale_color_manual(values=c("white","white","white","white","white"))+
  scale_fill_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  plot.theme+
  theme(aspect.ratio = 1,
        axis.text.y = element_text(hjust=1),
        legend.position="none")+
  labs(y="Diel pollination difference (SMD)")

ggsave(dn.eff.plot+ 
         draw_image("images/noun-day-1327879.png",
                    y = -5,x=0.15,scale=1.5)+
         draw_image("images/noun-night-1330487.png",
                    y = 4.1,x=0.1,scale=1.5), 
       file="figures/Figure 3 final.jpg", width=5, height=4, units="in", dpi=600)

###supplementary (Figure A2.1)
#Open pollination comparisons - effectiveness metric plot

eff.mod.plots2=lapply(eff.mod.tables[2:3],function (x) orchard_sans_pi(x,group="study_ID",legend.pos = "none",
                                                                 k.pos="left",mod = "treatment_effectiveness_metric",xlab="SMD",
                                                                 tree.order=(c("Seed set","Seed mass","Pollen deposition","Fruit set",
                                                                               "Fruit mass"))))

eff.mod.tables[["open_day"]]
eff.mod.tables[["open_night"]]

no.eff.plot=eff.mod.plots2[["open_night"]]+
  scale_color_manual(values=c("white","white","white","white","white"))+
  scale_fill_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  plot.theme+
  theme(aspect.ratio = 1,
        plot.title =element_text(size=10),
        axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none")+
  ggtitle("B) Night vs. open pollination")+
  labs(y="Diel pollination difference (SMD)")+
  xlab(NULL)+
  ylab(NULL)+
  labs(y="Diel pollination difference (SMD)")

do.eff.plot=eff.mod.plots2[["open_day"]]+
  scale_color_manual(values=c("white","white","white","white","white"))+
  scale_fill_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  plot.theme+
  theme(aspect.ratio = 1,
        plot.title =element_text(size=10),
        axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        axis.text.y = element_text(hjust=1,size=8),
        legend.position="none")+
  ggtitle("A) Day vs. open pollination")+
  labs(y="Diel pollination difference (SMD)")

open.eff.plot.out=cbind(ggplotGrob(do.eff.plot),ggplotGrob(no.eff.plot))

ggsave(open.eff.plot.out, 
       file="figures/Figure A2-3 final.jpg", 
       width=6, height=5, units="in", dpi=600)

###########
#Figure 4, Table S2, Figure A2.4
#Significance of plant functional traits and environmental variables on diel pollination differences 
###########

#rename categorical trait list to add "linear" for grouping
cat.trait.mod.list2=cat.trait.mod.list
names(cat.trait.mod.list2) = paste(names(cat.trait.mod.list2),"linear")

#combine all model lists
all.trait.mod.list <- c(cat.trait.mod.list, cont.env.trait.mod.list,poll.dep.mods)

#extract AIC, R2 and Q test from list of metafor models and put into a table
mod.table=list()

for (i in 1:length(all.trait.mod.list)){
  trait=word(names(all.trait.mod.list)[i],1)
  treatment=word(names(all.trait.mod.list)[i],2)
  linear=word(names(all.trait.mod.list)[i],3)
  delta=round(AIC(all.trait.mod.list[[i]])-AIC(overall.list[[treatment]]),3)
  QM=round(all.trait.mod.list[[i]]$QM,3)
  QMdf=round(all.trait.mod.list[[i]]$QMdf[1],3)
  QMp=all.trait.mod.list[[i]]$QMp
  QE=round(all.trait.mod.list[[i]]$QE,3)
  QEdf=round(all.trait.mod.list[[i]]$QEdf,3)
  QEp=all.trait.mod.list[[i]]$QEp
  r2m=round(r2_ml(all.trait.mod.list[[i]])[1] %>% as.numeric,3)
  r2c=round(r2_ml(all.trait.mod.list[[i]])[2] %>% as.numeric,3)
  
  mod.table[[i]]=data.frame(trait,treatment,linear,delta,QM,QMdf,QMp,QE,QEdf,QEp,r2m,r2c)
}

#combine to df
mod.table.df=data.table::rbindlist(mod.table)

###remove delta for poll.dep
mod.table.df=mod.table.df %>% 
  mutate(delta=ifelse(trait%in%"poll.dep",NA,delta))
#adjust p-values for multiple testing within each comparison
mod.table.df=mod.table.df %>% 
  group_by(treatment) %>% 
  mutate(n=length(QMp)) %>% 
  mutate(Qm.adj=p.adjust(QMp,method="BY")) %>% 
  mutate(Qm.adj.r=round(Qm.adj,4)) %>% 
  ungroup()

round(mod.table.df$Qm.adj,4)
#assign group (Trait or environment)
mod.table.df$env=ifelse(mod.table.df$trait%in%c("sDTR","sDaylength","sElevation"),
                        "Environment",
                        "Trait")

mod.table.df$trait=ifelse(mod.table.df$linear%in%"quadratic",paste(mod.table.df$trait,"^2"),mod.table.df$trait)

#mod.table.df$sign=ifelse(mod.table.df$delta<(-2),"*","")
mod.table.df$sign=ifelse(mod.table.df$Qm.adj<(0.05),"*","")


mod.table.df$treatment=revalue(mod.table.df$treatment,
                               c("day_night"="Day vs. night",
                                 "open_night"="B) Night vs. open pollination",
                                 "open_day"="A) Day vs. open pollination"))

mod.table.df$trait=revalue(mod.table.df$trait,
                           c("sDTR"="DTR",
                             "sDTR ^2"="DTR^2",
                             "sDaylength"="Daylength",
                             "sDaylength ^2"="Daylength^2",
                             "sElevation"="Elevation",
                             "sElevation ^2"="Elevation^2",
                             "flower_symmetry"="Flower symmetry",
                             "lifespan" ="Lifespan",
                             "life_form"       ="Life form",
                             "ps_pathway"     ="PS pathway",         
                             "flower_shape"       ="Flower shape",
                             "breeding_system"      ="Breeding system",
                             "bloom_period_simple"   ="Bloom period",
                             "nectar"           ="Nectar",
                             "odour"             ="Odour",     
                             "color"     ="Flower colour",          
                             "sStyle_length"="Style length",
                             "sPlant_height"="Plant height",
                             "poll.dep ^2" = "Pollination dependency^2",
                             "poll.dep" = "Pollination dependency"))

write.csv(mod.table.df,"output/Model summary table final.csv",row.names=FALSE)

##Figure 4
#Ordered bar graph of Marginal R^2 for each trait or environmental variable

mod.plot=mod.table.df %>% 
  filter(treatment%in%"Day vs. night") %>% 
  ggplot(aes(x=reorder(trait,r2m),y=r2m,#fill=env,
             alpha=sign))+
  geom_bar(stat="identity",fill="black")+
  #facet_wrap(~treatment)+
  labs(y=expression(R^2 * phantom() [M]),
       x=NULL)+
  #scale_fill_manual(values=c("#C582B2", "#51806a"),name="Variable type")+
  theme_bw()+
  coord_flip()+
  scale_x_discrete(labels=c("Daylength",
                            "Elevation",
                            "Nectar",
                            "Style length",
                            expression(Daylength^2),
                            "Plant height",
                            "Lifespan",
                            "DTR",
                            "Flower symmetry",
                            "Pollination dependency",
                            "Breeding system",
                            "Life form",
                            expression(DTR^2),
                            "Flower shape",
                            expression(Pollination~dependency^2),
                            "Photosynthetic pathway",
                            "Flower colour",
                            "Odour",
                            "Anthesis time",
                            expression(Elevation^2)))+#labels = c(""))
  scale_alpha_manual(values=c(0.5,1),guide="none")+
  plot.theme+
  theme(aspect.ratio=2,
        legend.position = "none")

ggsave(mod.plot,file="figures/Figure 4 final.jpg",width=5,height=5,unit="in",dpi=600)

####Figure A2.4
###Ordered bar graph of Marginal R^2 for each trait or environmental variable for both open pollination comparisons

#day vs open pollination
mod.bars2=mod.table.df %>% 
  filter(treatment%in%"A) Day vs. open pollination") %>% 
  ggplot(aes(x=reorder(trait,r2m),y=r2m,#fill=env,
             alpha=sign))+
  geom_bar(stat="identity",fill="black")+
  facet_wrap(~treatment,scales = "free")+
  labs(y=expression(R^2 * phantom() [M]),
       x=NULL)+
  #scale_fill_manual(values=c("#C582B2", "#51806a"),name="Variable type")+
  
  theme_bw()+
  coord_flip()+
  scale_x_discrete(labels=c("DTR",
                            "Flower symmetry",
                            "Daylength",
                            expression(Daylength^2),
                            expression(DTR^2),
                            "Elevation",
                            "Plant height",
                            "Breeding system",
                            "Lifespan",
                            "Nectar",
                            "Style length",
                            "Odour",
                            expression(Elevation^2),
                            "PS pathway",
                            "Anthesis time",
                            "Flower shape",
                            "Life form",
                            "Flower colour"))+#labels = c(""))
  scale_alpha_manual(values=c(0.5,1),guide="none")+
  theme(legend.position = "none",
        aspect.ratio=2,
        text=element_text(family="Times New Roman"),
        strip.background = element_blank(),
        strip.text = element_text(face="bold",hjust=0,size=8),
        title=element_text(face="bold",size=8),
        axis.text.y=element_text(size=7),
        axis.title=element_text(face="bold"))

#night vs open pollination
mod.bars3=mod.table.df %>% 
  filter(treatment%in%"B) Night vs. open pollination") %>% 
  ggplot(aes(x=reorder(trait,r2m),y=r2m,#fill=env,
             alpha=sign))+
  geom_bar(stat="identity",fill="black")+
  facet_wrap(~treatment,scales = "free")+
  labs(y=expression(R^2 * phantom() [M]),
       x=NULL)+
  #scale_fill_manual(values=c("#C582B2", "#51806a"),name="Variable type")+
  
  theme_bw()+
  coord_flip()+
  scale_x_discrete(labels=c("Daylength",
                            "Lifespan",
                            expression(Daylength^2),
                            "DTR",
                            "Nectar",
                            "PS pathway",
                            
                            "Flower symmetry",
                            "Style length",
                            expression(DTR^2),
                            "Breeding system",
                            "Plant height",
                            "Anthesis time",
                            "Elevation",
                            "Life form",
                            "Flower shape",
                            "Flower colour",
                            expression(Elevation^2),
                            "Odour"))+
  scale_alpha_manual(values=c(0.5,1),guide="none")+
  theme(legend.position = "right",
        aspect.ratio=2,
        text=element_text(family="Times New Roman"),
        strip.background = element_blank(),
        strip.text = element_text(face="bold",hjust=0,size=8),
        title=element_text(face="bold",size=8),
        axis.text.y=element_text(size=7),
        axis.title=element_text(face="bold"))

#combine and save as single plot
ggsave(cbind(ggplotGrob(mod.bars2),ggplotGrob(mod.bars3)),
       file="figures/Figure A2-4 final.jpg",width=6,height=4,unit="in",dpi=600)

###########

###########
###Figure 5, Figure A2.5
#Diel pollination differences in relation to elevation (m)  
######
#Figure 5

ele.dn.bub=mod_results(cont.env.trait.mod.list[["sElevation day_night quadratic"]], mod = "sElevation", at=list(exposure=1),
                       group = "study_ID",
                       weights = "prop")

e.dn.data_trim <- ele.dn.bub$data

#back-transform elevation from scale
ele.scale=es.list$dvn_effects %>% 
  filter(treatment%in%c("day_night","open_day","open_night"))%>%
  filter(yi > -10 & yi < 10) %>% 
  select(elevation) %>% 
  mutate(ele.scale=scale(elevation))

ele.scale=ele.scale$ele.scale

e.dn.data_trim$moderator=e.dn.data_trim$moderator*attr(ele.scale,"scaled:scale")+attr(ele.scale,"scaled:center")

e.dn.data_trim$scale <- (1/sqrt(e.dn.data_trim[,"vi"]))

e.legend <- "Precision (1/SE)"

e.dn.group_num <- as.vector(by(e.dn.data_trim, e.dn.data_trim[,"stdy"], function(x) base::length(base::unique(x[,"stdy"]))))

e.dn.dat_text <- data.frame(K = nrow(e.dn.data_trim), G = length(unique(e.dn.data_trim$stdy)))


e.dn.plot <-ggplot2::ggplot() +
  # putting bubbles
  ggplot2::geom_point(data = e.dn.data_trim, 
                      ggplot2::aes(x = moderator, y = yi, size = scale, fill = condition), 
                      fill = "darkgrey",shape = 21, col="white",alpha = 0.5,show.legend = F) +
  # confidence interval
  ggplot2::geom_ribbon(data = ele.dn.bub$mod_table%>% mutate(moderator=moderator*attr(ele.scale,"scaled:scale")+
                                                               attr(ele.scale,"scaled:center")), 
                       ggplot2::aes(x = moderator, ymin = lowerCL,ymax=upperCL,fill=condition), #method =  "loess", formula = y~x,se = FALSE,lty = "dashed", 
                       lwd = 1, alpha=0.25,fill = "black",
                       colour = NA) +
  # main line
  ggplot2::geom_smooth(data = ele.dn.bub$mod_table%>% mutate(moderator=moderator*attr(ele.scale,"scaled:scale")+
                                                               attr(ele.scale,"scaled:center")), 
                       ggplot2::aes(x = moderator, y = estimate,col=condition), 
                       col="black",method =  "loess", formula = y~x, se = FALSE, lwd = 1) +
  ggplot2::labs(x = "Elevation (m.a.s.l.)", y = "Diel pollination difference (SMD)", size = legend, parse = TRUE) +
  ggplot2::guides(fill = "none", colour = "none") +
  ggplot2::geom_hline(yintercept=0,linetype="dashed") +
  
  # themses
  ggplot2::theme_bw()  +
  ggplot2::theme(aspect.ratio = 1)+
  plot.theme +
  ggplot2::geom_text(data = e.dn.dat_text,
                     mapping = ggplot2::aes(x = Inf, y = -Inf),
                     label =  paste("italic(k)==", e.dn.dat_text$K,
                                    "~","(", e.dn.dat_text$G, ")"),
                     family="Times New Roman",
                     parse = TRUE,
                     hjust   = 1.1,
                     vjust   = -0.5)+ 
  draw_image("images/day.png",
             y = -5.075,x=0,scale=300) +
  draw_image("images/night.png",
             y = 3.75,x=20,scale=300)

e.dn.plot

ggsave(e.dn.plot,file="figures/Figure 5 final.jpg", width=5, height=5, units="in", dpi=600)

#supplementary elevation plots

ele.on.bub=mod_results(cont.env.trait.mod.list[["sElevation open_night quadratic"]], mod = "sElevation", at=list(exposure=0.5),
                       group = "study_ID",#by="treatment",
                       #at=list(treatment="open_day"),
                       weights = "prop")

e.on.data_trim <- ele.on.bub$data

e.on.data_trim$moderator=e.on.data_trim$moderator*attr(ele.scale,"scaled:scale")+attr(ele.scale,"scaled:center")

e.on.data_trim$scale <- (1/sqrt(e.on.data_trim[,"vi"]))

e.on.group_num <- as.vector(by(e.on.data_trim, e.on.data_trim[,"stdy"], function(x) base::length(base::unique(x[,"stdy"]))))

e.on.dat_text <- data.frame(K = nrow(e.on.data_trim), G = length(unique(e.on.data_trim$stdy)))

e.on.plot <-ggplot2::ggplot() +
  # putting bubbles
  ggplot2::geom_point(data = e.on.data_trim, ggplot2::aes(x = moderator, y = yi, size = scale, fill = condition), fill = "darkgrey",shape = 21, col="white",alpha = 0.5,show.legend = F) +
  # confidence interval
  ggplot2::geom_ribbon(data = ele.on.bub$mod_table%>% mutate(moderator=moderator*attr(ele.scale,"scaled:scale")+attr(ele.scale,"scaled:center")), ggplot2::aes(x = moderator, ymin = lowerCL,ymax=upperCL,fill=condition), #method =  "loess", formula = y~x,se = FALSE,lty = "dashed", 
                       lwd = 1, alpha=0.25,fill = "black",
                       colour = NA) +
  # main line
  ggplot2::geom_smooth(data = ele.on.bub$mod_table%>% mutate(moderator=moderator*attr(ele.scale,"scaled:scale")+attr(ele.scale,"scaled:center")), 
                       ggplot2::aes(x = moderator, y = estimate,col=condition), 
                       col="black",method =  "loess", formula = y~x, se = FALSE, lwd = 1) +
  ggplot2::labs(x = "Elevation (m)", y = "Diel pollination difference (SMD)", size = legend, parse = TRUE) +
  ggplot2::guides(fill = "none", colour = "none") +
  ggplot2::geom_hline(yintercept=0,linetype="dashed") +
  
  # themses
  ggplot2::theme_bw()  +
  ggplot2::theme(aspect.ratio = 1,
                 text=element_text(family="Times New Roman"),
                 strip.background = element_blank(),
                 strip.text = element_text(face="bold"),
                 title=element_text(face="bold",size=9),
                 axis.title=element_text(face="bold"),
                 axis.ticks.y=element_blank())+
  ggplot2::geom_text(data = e.on.dat_text,
                     mapping = ggplot2::aes(x = Inf, y = -Inf),
                     label =  paste("italic(k)==", e.on.dat_text$K,
                                    "~","(", e.on.dat_text$G, ")"),
                     family="Times New Roman",
                     parse = TRUE,
                     hjust   = 1.1,
                     vjust   = -0.5,
                     size=3)+
  ylim(-6,3.5)

ggsave(e.on.plot,file="figures/Figure A2-5 final.jpg", width=5, height=5, units="in", dpi=600)


#####
#Figure 6, Figures A2.6 and A2.7
#Trait effects on diel pollination differences
######

#Figure 6: effect of odour, bloom period (anthesis time), color and pollination dependency

#odour
odour_tab_out <- mod_results(cat.trait.mod.list[["odour day_night"]],
                             at=list(exposure=1),
                             group = "study_ID", 
                             mod = "odour")

odour_orchard=orchard_sans_pi(odour_tab_out, xlab = "Standardised Mean Difference",
                              legend.pos ="bottom.out",k.pos = c("left"),
                              tree.order = c("Present","Absent"))

#colour
colour_tab_out<- mod_results(cat.trait.mod.list[["color day_night"]], at=list(exposure=1),
                             group = "study_ID", 
                             mod = "color")

colour_orchard=orchard_sans_pi(colour_tab_out, xlab = "Diel pollination difference (SMD)",
                               legend.pos ="bottom.out",k.pos = c("left"),
                               tree.order = c("Brown_yellow","White","Red","Blue_purple","Pink","Orange",
                                              "Green"))

#bloom period
bloom_tab_out<- mod_results(cat.trait.mod.list[["bloom_period_simple day_night"]], at=list(exposure=1),
                            group = "study_ID", 
                            mod = "bloom_period_simple")

bloom_orchard=orchard_sans_pi(bloom_tab_out, xlab = "Diel pollination difference (SMD)",
                              legend.pos ="bottom.out",k.pos = c("left"),
                              tree.order = c("Night","Day","Both"))

###odour
odo.dn.plot=odour_orchard+
  scale_color_manual(values=c("white","white","white","white","white"))+
  scale_fill_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  plot.theme+
  theme(aspect.ratio = 0.5,
        text=element_text(family="Times New Roman"),
        title=element_text(face="bold",size=9),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        #subtitle=element_text(face="bold",size=8),
        legend.position = "none",
        #axis.text.y = element_text(size=8),
        axis.title=element_blank()
  )+
  labs(title = NULL,
       subtitle = "i) Flower odour")+ 
  draw_image("images/day.png",
             y = -5.5,x=0.15,scale=0.75)+
  draw_image("images/night.png",
             y = 3.8,x=0.15,scale=0.75)

#####colour
col.dn.plot=colour_orchard+
  scale_color_manual(values=c("white","white","white","white","white","white","white"))+
  scale_fill_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  plot.theme+
  theme(aspect.ratio = 1,
        text=element_text(family="Times New Roman"),
        title=element_text(face="bold",size=9),
        legend.position = "none",
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=9),
        axis.title.y=element_blank(),)+
  labs(title=NULL,subtitle = "iii) Flower colour")+
  scale_x_discrete(labels=c("Yellow","White","Red","Purple","Pink","Orange","Green"))+ 
  draw_image("images/day.png",
             y = -5.5,x=0.25,scale=0.75)+
  draw_image("images/night.png",
             y = 3.8,x=0.25,scale=0.75)

#####bloom period
bp.dn.plot=bloom_orchard+
  scale_color_manual(values=c("white","white","white","white","white","white","white"))+
  scale_fill_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  theme(aspect.ratio = 0.5,
        text=element_text(family="Times New Roman"),
        title=element_text(face="bold",size=9),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.y=element_blank(),)+
  
  labs(title=NULL,subtitle = "ii) Anthesis time")+ 
  draw_image("images/day.png",
             y = -5.5,x=0.25,scale=0.75)+
  draw_image("images/night.png",
             y = 3.8,x=0.25,scale=0.75)

ggsave(rbind(ggplotGrob(odo.dn.plot),ggplotGrob(bp.dn.plot),
             ggplotGrob(col.dn.plot)
             ),
       file="figures/Figure 6 final.jpg",width=4,height=7,dpi=600,units="in")


#Figure A2.6 - night vs open pollination differences in relation to traits


#odour
od_no_tab_out <- mod_results(cat.trait.mod.list[["odour open_night"]], at=list(exposure=0.5),
                             group = "study_ID", 
                             mod = "odour")

od_no_orchard=orchard_sans_pi2(od_no_tab_out, xlab = "Diel pollination difference (SMD)",
                               legend.pos ="bottom.out",k.pos = c("left"),
                               tree.order = c("Present","Absent"))

od.no.plot=od_no_orchard+
  scale_color_manual(values=c("white","white","white","white","white"))+
  scale_fill_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  theme(aspect.ratio = 0.5,
        text=element_text(family="Times New Roman",size=8),
        title=element_text(face="bold",size=9),
        legend.position = "none",
        axis.text.y = element_text(family="Times New Roman",size=8),
        axis.text.x=element_text(family="Times New Roman",size=8))+#,
        #axis.title.y=element_blank())+
  labs(title=NULL)+
  xlab("Odour")#,
     #  subtitle = "i) Odour")

ggsave(od.no.plot,
       file="figures/Figure A2-6 final.jpg",
       width=5,height=4,device="jpg",dpi=600)

######
####Funnel plots (Figure A1-3, A2-7)
######

#Day vs night pollination (Figure A1-3)

est1=pub.bias.ma.list[["day_night"]]$b[1]
se1=pub.bias.ma.list[["day_night"]]$se[1]

se.seq1=seq(0,max(sqrt(pub.bias.ma.list[["day_night"]]$data$vi)),length.out=100)

ll95.1=est1-1.96*se.seq1
ul95.1=est1+1.96*se.seq1

meanll95.1=est1-1.96*se1
meanul95.1=est1+1.96*se1

df.fun1=data.frame(ll95.1,ul95.1,se.seq1,est1,meanll95.1,meanul95.1)

dvn.funnel=ggplot(data=pub.bias.ma.list[["day_night"]]$data,
                  aes(y=yi,
                      x=sqrt(vi)))+
  geom_point(pch=21,fill="grey90",col="black",alpha=0.5,size=2)+
  geom_line(aes(x = se.seq1, y = ll95.1), linetype = 'dashed', data = df.fun1) +
  geom_line(aes(x = se.seq1, y = ul95.1), linetype = 'dashed', data = df.fun1) +
  
  geom_segment(aes(x = min(se.seq1), 
                   y = pub.bias.ma.list[["day_night"]]$b[1], xend = max(se.seq1), 
                   yend = pub.bias.ma.list[["day_night"]]$b[1]), linetype='dashed', data=df.fun1)+
  xlab("Standard error")+
  ylab("Diel pollination difference (SMD)")+
  scale_x_reverse()+
  coord_flip()+
  theme_bw()+
  theme(aspect.ratio = 1,
        text=element_text(family="Times New Roman",size=9),
        title=element_text(face="bold",size=10))

ggsave(dvn.funnel,file="figures/Figure A1-5 final.jpg",width=4,height=4,units="in",dpi=600)

####Figure A2.6

#open vs day pollination
est2=pub.bias.ma.list[["open_day"]]$b[1]
se2=pub.bias.ma.list[["open_day"]]$se[1]

se.seq2=seq(0,max(sqrt(pub.bias.ma.list[["open_day"]]$data$vi)),length.out=100)

ll95.2=est2-1.96*se.seq2
ul95.2=est2+1.96*se.seq2

meanll95.2=est2-1.96*se2
meanul95.2=est2+1.96*se2

df.fun2=data.frame(ll95.2,ul95.2,se.seq2,est2,meanll95.2,meanul95.2)

dvo.funnel=ggplot(data=pub.bias.ma.list[["open_day"]]$data,
                  aes(y=yi,
                      x=sqrt(vi)))+
  geom_point(pch=21,fill="grey90",col="black",alpha=0.5,size=2)+
  geom_line(aes(x = se.seq2, y = ll95.2), linetype = 'dashed', data = df.fun2) +
  geom_line(aes(x = se.seq2, y = ul95.2), linetype = 'dashed', data = df.fun2) +
  geom_segment(aes(x = min(se.seq2), 
                   y = pub.bias.ma.list[["open_day"]]$b[1], xend = max(se.seq2), 
                   yend = pub.bias.ma.list[["open_day"]]$b[1]), linetype='dashed', data=df.fun2)+
  xlab("Standard error")+
  ylab("Diel pollination difference (SMD)")+
  scale_x_reverse()+
  coord_flip()+
  theme_bw()+
  ggtitle("A) Day vs. open pollination")+
  theme(aspect.ratio = 1,
        text=element_text(family="Times New Roman",size=8),
        title=element_text(face="bold",size=8))

#open vs night
est3=pub.bias.ma.list[["open_night"]]$b[1]
se3=pub.bias.ma.list[["open_night"]]$se[1]

se.seq3=seq(0,max(sqrt(pub.bias.ma.list[["open_night"]]$data$vi)),length.out=100)

ll95.3=est3-1.96*se.seq3
ul95.3=est3+1.96*se.seq3

meanll95.3=est3-1.96*se3
meanul95.3=est3+1.96*se3

df.fun3=data.frame(ll95.3,ul95.3,se.seq3,est3,meanll95.3,meanul95.3)

nvo.funnel=ggplot(data=pub.bias.ma.list[["open_night"]]$data,
                  aes(y=yi,
                      x=sqrt(vi)))+
  geom_point(pch=21,fill="grey90",col="black",alpha=0.5,size=2)+
  geom_line(aes(x = se.seq3, y = ll95.3), linetype = 'dashed', data = df.fun3) +
  geom_line(aes(x = se.seq3, y = ul95.3), linetype = 'dashed', data = df.fun3) +
  geom_segment(aes(x = min(se.seq3), 
                   y = pub.bias.ma.list[["open_night"]]$b[1], xend = max(se.seq3), 
                   yend = pub.bias.ma.list[["open_night"]]$b[1]), linetype='dashed', data=df.fun3)+
  xlab(NULL)+
  ylab("Diel pollination difference (SMD)")+
  scale_x_reverse()+
  coord_flip()+
  theme_bw()+
  theme(aspect.ratio = 1,
        text=element_text(family="Times New Roman",size=8),
        title=element_text(face="bold",size=8))+
  ggtitle("B) Night vs. open pollination")

open.funs=cbind(ggplotGrob(dvo.funnel),ggplotGrob(nvo.funnel))

ggsave(open.funs,file="figures/Figure A2-7 final.jpg",width=6, height=4,units="in",dpi=600)

######reporting values from Egger tests and time-lag models (publication bias)

#main text
summary(egger.list[["day_night"]])

#supplement
summary(egger.list[["open_day"]])
summary(egger.list[["open_night"]])

#main text
summary(time.lag.list[["day_night"]])

#supplement
summary(time.lag.list[["open_day"]])
summary(time.lag.list[["open_night"]])

#####additional diel quotient plots
ratio.bub=mod_results(overall.list[[1]], mod = "exposure", 
                      group = "study_ID",
                      weights = "prop")

data_trim <- ratio.bub$data

data_trim$scale <- (1/sqrt(data_trim[,"vi"]))

legend <- "Precision (1/SE)"

group_num <- as.vector(by(data_trim, data_trim[,"stdy"], function(x) base::length(base::unique(x[,"stdy"]))))

dat_text <- data.frame(K = nrow(data_trim), G = length(unique(data_trim$stdy)))

ratio.plot <-ggplot2::ggplot() +
  # putting bubbles
  ggplot2::geom_point(data = data_trim, ggplot2::aes(x = moderator, y = yi, size = scale, fill = condition), fill = "darkgrey",shape = 21, col="white",alpha = 0.5,show.legend = F) +
  # confidence interval
  ggplot2::geom_ribbon(data = ratio.bub$mod_table, ggplot2::aes(x = moderator, ymin = lowerCL,ymax=upperCL,fill=condition), #method =  "loess", formula = y~x,se = FALSE,lty = "dashed", 
                       lwd = 1, alpha=0.25,fill = "black",
                       colour = NA) +
  # main line
  ggplot2::geom_smooth(data = ratio.bub$mod_table, 
                       ggplot2::aes(x = moderator, y = estimate,col=condition), 
                       col="black",method =  "loess", formula = y~x, se = FALSE, lwd = 1) +
  ggplot2::labs(x = "Diel quotient", y = "Diel pollination difference (SMD)", size = legend, parse = TRUE) +
  ggplot2::guides(fill = "none", colour = "none") +
  ggplot2::geom_hline(yintercept=0,linetype="dashed") +
  
  # themses
  ggplot2::theme_bw()  +
  plot.theme+
  ggplot2::theme(aspect.ratio = 1)+
  ggplot2::geom_text(data = dat_text,
                     mapping = ggplot2::aes(x = Inf, y = -Inf),
                     label =  paste("italic(k)==", dat_text$K,
                                    "~","(", dat_text$G, ")"),
                     family="Times New Roman",
                     parse = TRUE,
                     hjust   = 1.1,
                     vjust   = -0.5,
                     size=3)+ 
  draw_image("images/day.png",
             y = -5.075,x=-0.14,scale=0.75) +
  draw_image("images/night.png",
             y = 3.75,x=-0.15,scale=0.75)

ratio.plot

ggsave(ratio.plot,file="figures/Figure A1-3 final.jpg", width=5, height=5, units="in", dpi=600)
