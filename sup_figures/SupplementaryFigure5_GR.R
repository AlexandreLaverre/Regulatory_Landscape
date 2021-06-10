
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

pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_S5.pdf", sep=""), width=6.85, height=2.5)

par(mar=c(0.25, 1, 0.15, 0.25))
plot(shh.img, axes=F, ylim=c(300, 10), xlim=c(0, 1100))

mtext("SHH", side=3, font=3, adj=1, line=-1, at=400, cex=0.7)
mtext("expression pattern (HPA RNA-seq data)", side=3, font=1, adj=0, line=-1, at=410, cex=0.7)
mtext("expression level (TPM)", side=2, line=0, cex=0.7)

###########################################################################################

dev.off()

###########################################################################################
