######################################################################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

######################################################################################################################

if(load){
  sp="human"
  
  load(paste(pathFigures, "RData/", sp, ".cells.types.parallel.trends.RData", sep=""))

  load=FALSE
}

######################################################################################################################

cells <- c("Bcell", "ESC", "adipo")
cell.colors=c("firebrick1", "forestgreen", "navy")
names(cell.colors) = cells

######################################################################################################################

pdf(file=paste(pathFigures, "/SupplementaryMaterialFigure34.pdf", sep=""), width=6.85, height=3)

m=matrix(rep(NA, 1*14), nrow=1)
m[1,]=c(rep(1, 2), rep(2,4),  rep(3,4),  rep(4,4))
layout(m)

######################################################################################################################

par(mar=c(0,0,0,0))
plot.new()

######################################################################################################################

## distribution of expression correlations - 100 bootstrap replicats

par(mar=c(4.1, 1.25, 1.5, 1))
boxplot(exp.corr[cells], col="white", border=cell.colors, axes=F, names=rep("", 3), outline=F, horizontal=T, boxwex=0.55, notch=T)

axis(side=2, las=2, cex.axis=1.1, mgp=c(3, 0.5, 0), labels=rep("", 3), at=1:3)
mtext(c("B cells", "ESC", "pre-adipocytes"), side=2, line=0.75, cex=0.8, las=2, at=1:3)

axis(side=1, cex.axis=1.1)
mtext("expression correlation", side=1, line=2.85, cex=0.8)

box()

mtext("a", side=3, at=0.745, line=0.4, font=2, cex=1.1)

######################################################################################################################

## distribution of dN dS - top 25% expressed genes
par(mar=c(4.1, 1.25, 1.5, 1))
boxplot(dNdS[cells], col="white", border=cell.colors, axes=F, names=rep("", 3), outline=F, horizontal=T, boxwex=0.55, notch=T)

axis(side=2, las=2, cex.axis=1.1, mgp=c(3, 0.5, 0), labels=rep("", 3), at=1:3)

axis(side=1, cex.axis=1.1)
mtext("dN/dS ratio", side=1, line=2.85, cex=0.8)

box()

mtext("b", side=3, at=0.0, line=0.4, font=2, cex=1.1)

######################################################################################################################

par(mar=c(4.1, 1.25, 1.5, 1))
boxplot(enh.alignment[cells], col="white", border=cell.colors, axes=F, names=rep("", 3), outline=F, horizontal=T, boxwex=0.55, notch=T)

axis(side=2, las=2, cex.axis=1.1, mgp=c(3, 0.5, 0), labels=rep("", 3), at=1:3)

axis(side=1, cex.axis=1.1)
mtext("enhancer alignment score", side=1, line=2.85, cex=0.8)

box()

mtext("c", side=3, at=0.0, line=0.4, font=2, cex=1.1)

######################################################################################################################


dev.off()

##############################################################################################################
