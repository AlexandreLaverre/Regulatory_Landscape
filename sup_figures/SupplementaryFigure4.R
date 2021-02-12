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
  
  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.RData", sep=""))

  expdiv$EuclideanSimilarity <- 1-expdiv$EuclideanDistance
  
  expdiv$ClassRPKM <- cut(expdiv$MeanRPKM, breaks=quantile(expdiv$MeanRPKM, p=seq(from=0, to=1, length=6)), include.lowest=T)

  expdiv$MeanTau=(expdiv$TauHuman+expdiv$TauMouse)/2
  expdiv$ClassTau <- cut(expdiv$MeanTau, breaks=quantile(expdiv$MeanTau, p=seq(from=0, to=1, length=6)), include.lowest=T)
  
  if (sp == "human"){
    sp_name="Human"
  }else{
    sp_name="Mouse"
  }
  

  load=FALSE
}

######################################################################################################################
######################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

######################################################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure4.pdf", sep=""), width = 6.85, height=9)

m=matrix(rep(NA, 9*2), nrow=9)
for(i in 1:3){
  m[i,]=c(1,2)
}

for(i in 4:6){
  m[i,]=c(3,4)
}
for(i in 7:9){
  m[i,]=c(5,6)
}

layout(m)

cex.lab= 0.8
CEXstats = 0.8

######################################################################################################################
######################################################################################################################

## uncorrected Spearman's rho against class RPKM

par(mar=c(4.1, 4.5, 2.1, 1.5))
boxplot(expdiv$CorrelationSpearman~expdiv$ClassRPKM, outline=F, names=rep("", 5), border="black", col="white", notch=T, boxwex=0.75, xlab="", ylab="", axes=F)
axis(side=1, at=1:5, labels=c("low", "", "medium", "", "high"), mgp=c(3, 0.75, 0), cex.axis=1.3)
axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)

mtext("expression level class", side=1, line=2.5, cex=cex.lab)
mtext("Spearman's rho", side=2, line=3, cex=cex.lab)

mtext("a", side=3, at=-0.65, font=2, cex=1.1, line=1)

######################################################################################################################

## uncorrected Spearman's rho against class Tau

par(mar=c(4.1, 5.5, 2.1, 0.5))
boxplot(expdiv$CorrelationSpearman~expdiv$ClassTau, outline=F, names=rep("", 5), border="black", col="white", notch=T, boxwex=0.75, xlab="", ylab="", axes=F)
axis(side=1, at=1:5, labels=c("low", "", "medium", "", "high"), mgp=c(3, 0.75, 0), cex.axis=1.3)
axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)

mtext("expression specificity class", side=1, line=2.5, cex=cex.lab)
mtext("Spearman's rho", side=2, line=3, cex=cex.lab)

mtext("b", side=3, at=-0.65, font=2, cex=1.1, line=1)

######################################################################################################################

## uncorrected Euclidean similarity against class RPKM

par(mar=c(4.6, 4.5, 2.1, 1.5))
boxplot(expdiv$EuclideanSimilarity~expdiv$ClassRPKM, outline=F, names=rep("", 5), border="black", col="white", notch=T, boxwex=0.75, xlab="", ylab="", axes=F)
axis(side=1, at=1:5, labels=c("low", "", "medium", "", "high"), mgp=c(3, 0.75, 0), cex.axis=1.3)
axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)

mtext("expression level class", side=1, line=2.5, cex=cex.lab)
mtext("1-Euclidean distance", side=2, line=3, cex=cex.lab)

mtext("c", side=3, at=-0.65, font=2, cex=1.1, line=1)

######################################################################################################################


## uncorrected Euclidean similarity against class Tau

par(mar=c(4.6, 5.5, 2.1, 0.5))
boxplot(expdiv$EuclideanSimilarity~expdiv$ClassTau, outline=F, names=rep("", 5), border="black", col="white", notch=T, boxwex=0.75, xlab="", ylab="", axes=F)
axis(side=1, at=1:5, labels=c("low", "", "medium", "", "high"), mgp=c(3, 0.75, 0), cex.axis=1.3)
axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)

mtext("expression specificity class", side=1, line=2.5, cex=cex.lab)
mtext("1-Euclidean distance", side=2, line=3, cex=cex.lab)

mtext("d", side=3, at=-0.65, font=2, cex=1.1, line=1)

######################################################################################################################

## correlation between Spearman's rho and Euclidean similarity, before correction ################ 

par(mar=c(4.1, 4.5, 2.1, 1.5))

smoothScatter(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman, xlab="", ylab="", cex.axis=1.1)

R=cor(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman,method="pearson")
rho=cor(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman, method="spearman")
x=expdiv$EuclideanSimilarity
y=expdiv$CorrelationSpearman
abline(lm(y~x), col="red")

mtext(paste("Pearson's R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("1-Euclidean distance", side=1, line=2.5, cex=CEXLAB)
mtext("Spearman's rho", side=2, line=2.75, cex=CEXLAB)
mtext("e", side=3, line=2, at=-0.01, font=2, cex=1.2)

##############################################################################################

## after correction
smoothScatter(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman, xlab="", ylab="", cex.axis=1.1)
R=cor(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman,method="pearson")
rho=cor(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman, method="spearman")
x=expdiv$CorrectedEuclideanSimilarity
y=expdiv$CorrectedSpearman
abline(lm(y~x), col="red")

mtext(paste("Pearson's R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("1-Euclidean distance (corrected)", side=1, line=2.5, cex=CEXLAB)
mtext("Spearman's rho (corrected)", side=2, line=2.75, cex=CEXLAB)
mtext("f", side=3, line=2, at=-0.78, font=2, cex=1.2)

##############################################################################################

dev.off()

##############################################################################################

