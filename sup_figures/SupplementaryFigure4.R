######################################################################################################################

library(Hmisc)

options(stringsAsFactors = FALSE)

#######################################################################################################################

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
  
  load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.RData", sep=""))
  
  if (sp == "human"){
    sp_name="Human"
  }else{
    sp_name="Mouse"
  }
  
  bord.measure = c("#AF46B4", "#4BB446")
  names(bord.measure) = c("Spearman's rho", "1-Euclidean distance")
  col.measure = c("#AF46B41A", "#4BB4461A")
  names(col.measure) = c("Spearman's rho", "1-Euclidean distance")
  nb = 1

  load=FALSE
}

######################################################################################################################

clearboxplot <- function(measure, divergences, plotlabel){
  
  if (measure == "classTau"){
    labels=c("broad", "", "", "narrow")
    xlab="expression breadth"
  } else{
    labels=c("low", "", "", "high")
    xlab="mean expression level"
  }
  
  xlim=c(0.5, length(levels(expdiv[[measure]]))+0.5)
  
  xpos=seq(1, length(levels(expdiv[[measure]])), 1)
  names(xpos) = levels(expdiv[[measure]])
  smallx=c(-0.15, 0.15)
  names(smallx)=c("Spearman's rho", "1-Euclidean distance")
  
  for (divergence in divergences){
    if (any(grepl("Spearman", divergence))){
      name="Spearman's rho"
      ylim = c(-0.5, 1)
    } else{
      name="1-Euclidean distance"
      ylim = c(0.55, 1)
    }
    
    if (name=="Spearman's rho"){
      par(mar=c(1.5, 3.75, 2.6, 1.1))
    } else{
      par(mar=c(4, 3.75, 0.6, 1.1))
    }

    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    for(class in levels(expdiv[[measure]])){
      x=xpos[class]+smallx[name]
    
      boxplot(expdiv[which(expdiv[[measure]] == class), divergence], at=x, border=bord.measure[name], col=col.measure[name],
             boxwex=0.5, axes=F, add=T, notch=T, outline=F)
    }
    
    abline(v=xpos[1:length(levels(expdiv[[measure]]))-1]+0.5, lty=3, col="gray40")
    axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
    mtext(name, side=2, line=2.5, cex=CEXLAB)
    
    if (name=="Spearman's rho"){
      mtext(plotlabel, side=3, line=1, at=-0.1, font=2, cex=1.2)
    }else{
      axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(expdiv[[measure]]))), cex.axis=0.8)
      mtext(labels, at=xpos, side=1, line=1, cex=CEXLAB)
      
      mtext(xlab, side=1, line=2.5, cex=CEXLAB)
    }
  }
}

######################################################################################################################
######################################################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure4.pdf", sep=""), width = 6.85)

m=matrix(rep(NA, 10*2), nrow=10)
for(i in 1:3){
  m[i,]=c(1,3)
}

for(i in 4:6){
  m[i,]=c(2,4)
}
for(i in 7:10){
  m[i,]=c(5,6)
}

layout(m)

CEXLAB = 0.9
CEXstats = 0.8
expdiv$EuclideanSimilarity <- 1-expdiv$EuclideanDistance

################ Cofounding factors ################ 
# Gene expression level
expdiv$classRPKM <- cut2(expdiv$Human_MeanRPKM, g=4, include.lowest=T) 
clearboxplot("classRPKM", c("CorrelationSpearman","EuclideanSimilarity"), "a")

# Specificity
clearboxplot("classTau", c("CorrelationSpearman","EuclideanSimilarity"), "b")

################ Correlation between Spearman and Euclidean  ################ 
par(mar=c(4.5, 5, 3, 5)) # bottom, left, top, right

smoothScatter(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman, xlab="", ylab="", cex.axis=1.1)

R=cor(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman,method="pearson")
rho=cor(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman, method="spearman")
x=expdiv$EuclideanSimilarity
y=expdiv$CorrelationSpearman
abline(lm(y~x), col="red")

mtext(paste("Pearson's R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("1-Euclidean distance", side=1, line=2.5, cex=CEXLAB)
mtext("Spearman's rho", side=2, line=2.5, cex=CEXLAB)
mtext("c", side=3, line=2, at=-0.01, font=2, cex=1.2)

### Correlation between measures
smoothScatter(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman, xlab="", ylab="", cex.axis=1.1)
R=cor(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman,method="pearson")
rho=cor(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman, method="spearman")
x=expdiv$CorrectedEuclideanSimilarity
y=expdiv$CorrectedSpearman
abline(lm(y~x), col="red")

mtext(paste("Pearson's R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("1-Euclidean distance (corrected)", side=1, line=2.5, cex=CEXLAB)
mtext("Spearman's rho (corrected)", side=2, line=2.5, cex=CEXLAB)
mtext("d", side=3, line=2, at=-0.8, font=2, cex=1.2)

##############################################################################################

dev.off()

##############################################################################################

