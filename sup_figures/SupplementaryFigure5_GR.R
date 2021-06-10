
#############################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#############################################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R") ## paths are defined based on the user name
  
  library(imager)
  
  set.seed(19)
}

###########################################################################################

if(prepare){
  shh.img=load.image("img/SHH_expression.png")
}

###########################################################################################

pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_S5.pdf", sep=""), width=6.85, height=3.5)

par(mar=c(0.25, 1, 0.15, 0.25))
plot(shh.img, axes=F, ylim=c(210, 20), xlim=c(0, 400))

###########################################################################################

dev.off()

###########################################################################################
