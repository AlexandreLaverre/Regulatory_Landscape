######################################################################################################################
library(Hmisc)
options(stringsAsFactors = FALSE)
source("../main_figures/parameters.R") ## pathFinalData are defined based on the user name

sp="human"

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.Rdata", sep=""))
load(paste(pathFigures, "RData/data.", sp, ".regland.conservation.RData", sep=""))

if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}

bord.measure = c("#AF46B4", "#4BB446")
names(bord.measure) = c("Spearman", "Euclidean")
col.measure = c("#AF46B41A", "#4BB4461A")
names(col.measure) = c("Spearman", "Euclidean")

clearboxplot <- function(measure, divergences, ylim, ylab){
  
  if (measure == "classTau"){labels=c("generalist", "", "", "specialist"); xlab="Expression Specificity"
  }else{labels=c("low", "", "", "high"); xlab="Mean Expression Level"}
  
  xlim=c(0.5, length(levels(expdiv[[measure]]))+0.5)
  
  xpos=seq(1, length(levels(expdiv[[measure]])), 1)
  names(xpos) = levels(expdiv[[measure]])
  smallx=c(-0.15, 0.15)
  names(smallx)=c("Spearman", "Euclidean")
  
  par(mar=c(4.5, 3.75, 3, 1.1))
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for (divergence in divergences){
    if (any(grepl("Spearman", divergence))){name="Spearman"}else{name="Euclidean"}
    
    for(class in levels(expdiv[[measure]])){
      x=xpos[class]+smallx[name]
      b = boxplot(expdiv[which(expdiv[[measure]] == class), divergence], plot=F)
      med=b$stats[3]
      ci=as.numeric(b$conf)
      points(x, med, pch=20, col=bord.measure[name], cex=1)
      segments(x, ci[1], x, ci[2],  col=bord.measure[name])
      
      # boxplot(expdiv[which(expdiv[[measure]] == class), divergence], at=x, border=bord.measure[name], col=col.measure[name],
      #         boxwex=0.5, axes=F, add=T, notch=T, outline=F)
    }
  }
  
  abline(v=xpos[1:length(levels(expdiv[[measure]]))-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(expdiv[[measure]]))), cex.axis=0.8)
  mtext(labels, at=xpos, side=1, line=1, cex=CEXLAB)
  mtext(xlab, side=1, line=2.5, cex=CEXLAB)
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
  mtext(ylab, side=2, line=2.5, cex=CEXLAB)
  
}

######################################################################################################################
pdf(file=paste(pathFigures, "ExtendedFigure4.pdf", sep=""), width = 8.5)

par(mfrow=c(2,2))
CEXLAB = 0.9
CEXstats = 0.8
expdiv$EuclideanSimilarity <- 1-expdiv$EuclideanDistance

################ Cofounding factors ################ 
# Gene expression level
expdiv$classRPKM <- cut2(expdiv$Human_MeanRPKM, g=4, include.lowest=T) 
clearboxplot("classRPKM", c("CorrelationSpearman","EuclideanSimilarity"), c(0.4, 1),  "Profile Similarity")

legend("left", legend=c("1 - Euclidean Distance", "Spearman's rho"), pch=20, col=rev(bord.measure), cex=1, 
       box.col="white", bg="white",  inset=c(0.01, 0.01))

mtext("a", side=3, line=1, at=0.1, font=2, cex=1.2)

# Specificity
clearboxplot("classTau", c("CorrelationSpearman","EuclideanSimilarity"), c(0.3, 1),  "Profile Similarity")
mtext("b", side=3, line=1, at=0.1, font=2, cex=1.2)

################ Correlation between Spearman and Euclidean  ################ 
smoothScatter(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman, xlab="", ylab="")

R=cor(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman,method="pearson")
rho=cor(expdiv$EuclideanSimilarity, expdiv$CorrelationSpearman, method="spearman")
x=expdiv$EuclideanSimilarity
y=expdiv$CorrelationSpearman
abline(lm(y~x), col="red")

mtext(paste("R2 = ", round(summary(lm(y~x))$r.squared, digits=2), ", Pearson R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("Euclidean Similarity", side=1, line=2.5, cex=CEXLAB)
mtext("Spearman's rho", side=2, line=2.5, cex=CEXLAB)
mtext("c", side=3, line=1, at=0.08, font=2, cex=1.2)

### Correlation between measures
smoothScatter(expdiv$CorrectedEuclideanSimilarityExpTau, expdiv$CorrectedSpearmanExpTau, xlab="", ylab="")
R=cor(expdiv$CorrectedEuclideanSimilarityExpTau, expdiv$CorrectedSpearmanExpTau,method="pearson")
rho=cor(expdiv$CorrectedEuclideanSimilarityExpTau, expdiv$CorrectedSpearmanExpTau, method="spearman")
x=expdiv$CorrectedEuclideanSimilarityExpTau
y=expdiv$CorrectedSpearmanExpTau
abline(lm(y~x), col="red")

mtext(paste("R2 = ", round(summary(lm(y~x))$r.squared, digits=2), ", Pearson R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("Corrected Euclidean Similarity", side=1, line=2.5, cex=CEXLAB)
mtext("Corrected Spearman's rho", side=2, line=2.5, cex=CEXLAB)
mtext("d", side=3, line=1, at=-0.72, font=2, cex=1.2)

dev.off()

