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
names(bord.measure) = c("Spearman's rho", "Euclidean Similarity")
col.measure = c("#AF46B41A", "#4BB4461A")
names(col.measure) = c("Spearman's rho", "Euclidean Similarity")
nb = 1

clearboxplot <- function(measure, divergences, nb){
  
  if (measure == "classTau"){labels=c("generalist", "", "", "specialist"); xlab="Expression Specificity"
  }else{labels=c("low", "", "", "high"); xlab="Mean Expression Level"}
  
  xlim=c(0.5, length(levels(expdiv[[measure]]))+0.5)
  
  xpos=seq(1, length(levels(expdiv[[measure]])), 1)
  names(xpos) = levels(expdiv[[measure]])
  smallx=c(-0.15, 0.15)
  names(smallx)=c("Spearman's rho", "Euclidean Similarity")
 
  for (divergence in divergences){
    if (any(grepl("Spearman", divergence))){name="Spearman's rho";  ylim = c(-0.5, 1)}else{name="Euclidean Similarity";  ylim = c(0.55, 1)}
    
    if (name=="Spearman's rho"){par(mar=c(0.6, 3.75, 3, 1.1)); }else{par(mar=c(4.5, 3.75, 0.6, 1.1))}

    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    for(class in levels(expdiv[[measure]])){
      x=xpos[class]+smallx[name]
      # b = boxplot(expdiv[which(expdiv[[measure]] == class), divergence], plot=F)
      # med=b$stats[3]
      # ci=as.numeric(b$conf)
      # points(x, med, pch=20, col=bord.measure[name], cex=1)
      # segments(x, ci[1], x, ci[2],  col=bord.measure[name])
      
      boxplot(expdiv[which(expdiv[[measure]] == class), divergence], at=x, border=bord.measure[name], col=col.measure[name],
             boxwex=0.5, axes=F, add=T, notch=T, outline=F)
    }
    
    abline(v=xpos[1:length(levels(expdiv[[measure]]))-1]+0.5, lty=3, col="gray40")
    axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
    mtext(name, side=2, line=2.5, cex=CEXLAB)
    
    if (name=="Spearman's rho"){
      mtext(nb, side=3, line=1, at=0.1, font=2, cex=1.2)
    }else{
      axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(expdiv[[measure]]))), cex.axis=0.8)
      mtext(labels, at=xpos, side=1, line=1, cex=CEXLAB)
      mtext(xlab, side=1, line=2.5, cex=CEXLAB)
    }
  }
  
}

######################################################################################################################
pdf(file=paste(pathFigures, "ExtendedFigure4.pdf", sep=""), width = 6.85)

m=matrix(rep(NA, 3*2), nrow=3)
m[1,]=c(1,3)
m[2,]=c(2,4)
m[3,]=c(5,6)

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
mtext("c", side=3, line=2, at=-0.01, font=2, cex=1.2)

### Correlation between measures
smoothScatter(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman, xlab="", ylab="")
R=cor(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman,method="pearson")
rho=cor(expdiv$CorrectedEuclideanSimilarity, expdiv$CorrectedSpearman, method="spearman")
x=expdiv$CorrectedEuclideanSimilarity
y=expdiv$CorrectedSpearman
abline(lm(y~x), col="red")

mtext(paste("R2 = ", round(summary(lm(y~x))$r.squared, digits=2), ", Pearson R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("Corrected Euclidean Similarity", side=1, line=2.5, cex=CEXLAB)
mtext("Corrected Spearman's rho", side=2, line=2.5, cex=CEXLAB)
mtext("d", side=3, line=2, at=-0.85, font=2, cex=1.2)

dev.off()

