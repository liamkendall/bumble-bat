########################################################################################################################################################
#SCRIPT TO CALCULATE THE PHLOGENETIC INFORMED PRINCIPAL COMPONENT ANALYSIS
#Just filled rows for the different quantitative variables (no imputed data)
########################################################################################################################################################
#LOAD LIBRARIES
library(phytools) #ppca
library(ape) #for phylogenetic distance
library(dplyr) #data processing
library(rtrees) #for phylogenetic distancelibrary(MASS)
library(reshape2) #data processing
library(viridis) #COLOUR GGPLOT
library(MASS) #I think I used it for the kernel density of the plotting function; no longer used but I leave in case its handy later on
library(ggplot2) #plotting
library(broman) #crayon colours

#################################################################################################################################
#4) CALCULATE PPCA
########################################################################################################################################################

dn.pca.tree=dn.tree
dn.pca.tree=drop.tip(dn.tree, setdiff(
  dn.tree$tip.label,rownames(traits.pca)))

traits.pca=traits %>% select_if(is.numeric) %>% 
  mutate_if(is.numeric,scale)

rownames(traits.pca)=gsub(" ","_",rownames(traits.pca))

setdiff(
  rownames(traits.pca),dn.pca.tree$tip.label)
#opposite
setdiff(
  dn.pca.tree$tip.label,rownames(traits.pca))

dvn_phyl_pca <- phyl.pca(dn.tree, traits.pca,
                                           method="lambda",mode="cov")

##save pca as rdata in data folder
save(dvn_phyl_pca,file="data/dvn_phyl_pca.rdata")

#CALL the output PC for simplicity
PC <- dvn_phyl_pca
#CHECK CONTENT
#EIGENVALUES
PC$Eval
#PC score (POINTS)
PC$S

traits.out=traits

traits.out$PC1 <- PC$S[,1]
traits.out$PC2 <- PC$S[,2]
traits.out$PC3 <- PC$S[,3]
traits.out$PC4 <- PC$S[,4]

#PC loadings (ARROWS)
PC$L

nrow(PC$S)

percentage <- round(diag(PC$Eval) / sum(PC$Eval) * 100, 2) #calculate percentage

sum(percentage[1] + percentage[2]) #73.62% of the variance is explained by the first two principal components
#first three
sum(percentage[1] + percentage[2] + percentage[3]) #94.8% of the variance is explained by the first three principal components
########################################################################################################################################################
#4) PLOT PPCA
########################################################################################################################################################

#scale numeric variables in trait df
dat_cleaning_5 <- traits.pca 

#plot same order as the other plots
#PC$S[,2] <- -PC$S[,2]


all_spp <- function(PC, x="PC1", y="PC2") {
 # x="PC1"
#  y="PC2"
  # PC being a prcomp object
  data <- data.frame(PC$S)
  plot <- ggplot(data, aes_string(x=x, y=y)) #generate plot
  dat <- data.frame(x = data[,x], y = data[,y])
  
  #######
  #DENSITY FUNCTION
  #######
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  dat$density <- get_density(dat$x, dat$y, h = c(2, 2), n = 1000) #obtain density
  
  
  plot <- plot+stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',bins=8) + 
    scale_fill_continuous(low="green",high="red",guide=FALSE,lim=c(0.01,0.082),breaks = c(0.02, 0.035, 0.082), 
                          labels = c("Low", "Medium", "High")) + theme(legend.position = "none")+
    guides(fill = guide_legend(override.aes = list(alpha = 0.6),title="Kernel density"))+
    scale_alpha(guide = 'none')
  
  plot <- plot + geom_point(data=dat, aes(x, y),size=0.95, color="black")
  
  ########
  #ADD ARROWS 
  ########
  datapc <- data.frame(PC$L) #CREATE DATAFRAME WITH LOADINGS
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  # add arrows
  plot <- plot + geom_segment(data=datapc,linejoin="round", lineend="round",
                              aes(x=0, y=0, xend=v1, yend=-v2),size=1.8, 
                              arrow=arrow(length=unit(0.5,"cm")), alpha=0.8, colour=c("black"))
  
  
  #ADD THE OTHER DIRECTION OF THE SEGMENT BECAUSE LOOKS COOL
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, 
                                               xend=-v1, yend=v2),size=1.6, 
                              arrow=arrow(length=unit(0,"cm")),
                              linetype=2, alpha=0.5, colour=c("black"))
  
  #Add axis with perctentage
  percentage <- round(diag(PC$Eval) / sum(PC$Eval) * 100, 2) #calculate percentage
  
  plot <- plot + xlab(paste("PC1 ", "(",(percentage[1]),"%",")", sep = "")) #XLAB
  plot <- plot + ylab(paste("PC2 ", "(",(percentage[2]),"%",")", sep = "")) #YLAB
  
  
  #CHANGE THEME
  
  plot <- plot + theme_bw() 
  
  #ADD LABELS
  rownames(PC$L) <- c("Flower width", "Flower length", "Style length", "Plant height")
  PCAloadings <- data.frame(Variables = rownames(PC$L), PC$L)
  plot <- plot + annotate("text", x = (PCAloadings$PC1*3), 
                          y = -(PCAloadings$PC2*3),#+c(0,0,0,0,0.4,0)),
                          label = PCAloadings$Variables, color="black",fontface =2,size=4)
  
  plot <- plot + theme_bw() +#ylim(-4,4) + xlim(-4,4) +  
    theme(legend.position = "none",
          aspect.ratio = 1)# c(0.095, 0.11)) +ggtitle("") 
  
  
  
  plot
  
}

all_spp(PC, x="PC1", y="PC2")
all_spp(PC, x="PC2", y="PC3")

