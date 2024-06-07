#title: "Pollination across the diel cycle: a global meta-analysis"
#authors: L Kendall & C.C. Nicholson
#date: "05/06/2024"

###script: 003A - Assessment of correlation among predictor variables

#additional libraries
library(tidyverse)
library(rcompanion)
require(corrr)
library(ggcorrplot)


# Calculate a pairwise association between all variables in a data-frame. In particular nominal vs nominal with Chi-square, numeric vs numeric with Pearson correlation, and nominal vs numeric with ANOVA.
# Adopted from https://stackoverflow.com/a/52557631/590437
mixed_assoc = function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
  df_comb = expand.grid(names(df), names(df),  stringsAsFactors = F) %>% set_names("X1", "X2")
  
  is_nominal = function(x) class(x) %in% c("factor", "character")
  # https://community.rstudio.com/t/why-is-purr-is-numeric-deprecated/3559
  # https://github.com/r-lib/rlang/issues/781
  is_numeric <- function(x) { is.integer(x) || is_double(x)}
  
  f = function(xName,yName) {
    x =  pull(df, xName)
    y =  pull(df, yName)
    
    result = if(is_nominal(x) && is_nominal(y)){
      # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
      cv = cramerV(as.character(x), as.character(y), bias.correct = adjust_cramersv_bias)
      data.frame(xName, yName, assoc=cv, type="cramersV")
      
    }else if(is_numeric(x) && is_numeric(y)){
      correlation = cor(x, y, method=cor_method, use="complete.obs")
      data.frame(xName, yName, assoc=correlation, type="correlation")
      
    }else if(is_numeric(x) && is_nominal(y)){
      # from https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
      r_squared = summary(lm(x ~ y))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else if(is_nominal(x) && is_numeric(y)){
      r_squared = summary(lm(y ~x))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else {
      warning(paste("unmatched column type combination: ", class(x), class(y)))
    }
    
    # finally add complete obs number and ratio to table
    result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
  }
  
  # apply function to each variable combination
  map2_df(df_comb$X1, df_comb$X2, f)
}

#traits
traits.cor.out=mixed_assoc(traits %>% select(-c(bloom_period,plant_family,old_name,plant_species)))%>%
  #if assoc is negative, make it positive
  mutate(assoc=ifelse(assoc<0,assoc*-1,assoc))

##environment (study-site level)
environment.cors=environment %>% select(study_ID,Daylength,elevation,DTR) %>% 
  distinct(study_ID,Daylength,elevation,DTR) %>%
  select(-study_ID) %>% 
  mixed_assoc()%>%
  #if assoc is negative, make it positive
  mutate(assoc=ifelse(assoc<0,assoc*-1,assoc))

####between environment and traits (species-site level)
env.cors=es.list$dvn_effects %>% select(-c(plant_family, bloom_period)) %>%
  select(accepted_name,elevation,DTR,Daylength,
         flower_symmetry:plant_height_midpoint_m) %>% #colnames
  distinct(accepted_name,elevation,DTR,Daylength,flower_symmetry,
  lifespan,life_form,ps_pathway,flower_shape,breeding_system,
  bloom_period_simple,nectar,odour,color,flower_width_midpoint_mm,
  flower_length_midpoint_mm,style_length_midpoint_mm,plant_height_midpoint_m)

env.cors.out=mixed_assoc(env.cors %>% select(-accepted_name))

env.cors.out2=env.cors.out %>% 
  #filter correlations that contain elevation, DTR and Daylength
  filter(x %in% c("elevation","DTR","Daylength") | y %in% c("elevation","DTR","Daylength")) %>% 
  #filter out correlations betweeen elevation, DTR and Daylength
  filter(!(x %in% c("elevation","DTR","Daylength") & y %in% c("elevation","DTR","Daylength"))) %>%
  #if assoc is negative, make it positive
  mutate(assoc=ifelse(assoc<0,assoc*-1,assoc))

traits.cor.all=traits.cor.out %>%
  rbind(env.cors.out2) %>% 
  rbind(environment.cors)

traits.cor.all$x=revalue(traits.cor.all$x,
                         c("DTR"="DTR",
                           "Daylength"="Daylength",
                           "elevation"="Elevation",
                           "flower_symmetry"="Flower symmetry",
                           "lifespan" ="Lifespan",
                           "life_form"       ="Life form",
                           "ps_pathway"     ="PS pathway",         
                           "flower_shape"       ="Flower shape",
                           "breeding_system"      ="Breeding system",
                           "bloom_period_simple"   ="Anthesis time",
                           "nectar"           ="Nectar",
                           "odour"             ="Odour",     
                           "color"     ="Flower colour",  
                           "flower_width_midpoint_mm" = "Flower width",
                           "flower_length_midpoint_mm"="Flower length",
                           "style_length_midpoint_mm"="Style length",
                           "plant_height_midpoint_m"="Plant height"))


traits.cor.all$y=revalue(traits.cor.all$y,
                         c("DTR"="DTR",
                           "Daylength"="Daylength",
                           "elevation"="Elevation",
                           "flower_symmetry"="Flower symmetry",
                           "lifespan" ="Lifespan",
                           "life_form"       ="Life form",
                           "ps_pathway"     ="PS pathway",         
                           "flower_shape"       ="Flower shape",
                           "breeding_system"      ="Breeding system",
                           "bloom_period_simple"   ="Anthesis time",
                           "nectar"           ="Nectar",
                           "odour"             ="Odour",     
                           "color"     ="Flower colour",  
                           "flower_width_midpoint_mm" = "Flower width",
                           "flower_length_midpoint_mm"="Flower length",
                           "style_length_midpoint_mm"="Style length",
                           "plant_height_midpoint_m"="Plant height"))

traits.cor.mat=traits.cor.all%>% 
  select(x, y, assoc) %>%
  spread(y, assoc) 

rownames(traits.cor.mat)=traits.cor.mat$x

traits.cor.mat=traits.cor.mat %>% select(-x) %>% as.matrix()

mat.ord=traits.cor.mat[c(3,4,5,7,10,17,15,1,2,6,8:9,11:14,16),c(3,4,5,7,10,17,15,1,2,6,8:9,11:14,16)]

mat.plot.out=ggcorrplot(mat.ord,type="upper",lab = TRUE,lab_size=2,legend.title = "Association")+ 
    scale_fill_gradient2(low = "white", high = "black", 
                         breaks=c(0, 0.25,0.5,0.75, 1), limit=c(0, 1))+
  theme_bw()+
  guides(fill=guide_legend(title="Association"))+
  plot.theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")+
  theme(axis.title=element_blank(),
        axis.ticks = element_blank())

ggsave(mat.plot.out,
       file="figures/Figure A1-2 June5.jpg",
       width=6,height=6,dpi=600)

