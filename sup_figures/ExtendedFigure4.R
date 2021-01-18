######################################################################################################################
options(stringsAsFactors = FALSE)
source("../main_figures/parameters.R") ## pathFinalData are defined based on the user name

sp="human"

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.Rdata", sep=""))
load(paste(pathFigures, "RData/data.", sp, ".regland.conservation.RData", sep=""))

load(paste(pathFigures, "RData/data.", sp, ".common.cells.expdiv.Rdata", sep=""))
load(paste(pathFigures, "RData/data.", sp, ".common.cells.regland.conservation.RData", sep=""))

if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}

cells <- c("Bcell", "ESC", "adipo")
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

clearboxplot <- function(measure, divergence, ylim, ylab){
  
  if (measure == "classTau"){labels=c("generalist", "", "", "specialist"); xlab="Expression Specificity"
  }else{labels=c("low", "", "", "high"); xlab="Mean Expression Level"}
  
  xlim=c(0.5, length(levels(expdiv[[measure]]))+0.5)
  
  xpos=seq(1, length(levels(expdiv[[measure]])), 1)
  names(xpos) = levels(expdiv[[measure]])
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  par(mar=c(4.5, 3.75, 3, 1.1))
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(class in levels(expdiv[[measure]])){
    x=xpos[class]
    boxplot(expdiv[which(expdiv[[measure]] == class), divergence], at=x, axes=F, add=T, notch=T, outline=F)
  }
  
  abline(v=xpos[1:length(levels(expdiv[[measure]]))-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(expdiv[[measure]]))), cex.axis=0.8)
  mtext(labels, at=xpos, side=1, line=1, cex=CEXLAB)
  mtext(xlab, side=1, line=2.5, cex=CEXLAB)
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
  mtext(ylab, side=2, line=2.5, cex=CEXLAB)
  
}

######################################################################################################################
#pdf(file=paste(pathFigures, "ExtendedFigure4.pdf", sep=""), width = 8.5)

par(mfrow=c(3,3))

######################################################################################################################
############################################ Gene expression level  ##################################################
#### Divergence (Common cells types)
CEXLAB = 0.9
CEXstats = 0.8

for (cell in cells){
  expdiv_cells[[paste0(cell, "_log2Expression")]] = log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1)
  expdiv_cells[[paste0(cell ,"_classRPKM")]]=cut2(expdiv_cells[[paste0(cell, "_log2Expression")]], g=4, include.lowest=T)
  
  smoothScatter(expdiv_cells[[paste0(cell, "_log2Expression")]], expdiv_cells[[paste0(cell, "_ExpressionConservation")]], xlab="", ylab="")

  R=cor(expdiv_cells[[paste0(cell, "_log2Expression")]], expdiv_cells[[paste0(cell, "_ExpressionConservation")]],method="pearson")
  rho=cor(expdiv_cells[[paste0(cell, "_log2Expression")]], expdiv_cells[[paste0(cell, "_ExpressionConservation")]], method="spearman")
  x=expdiv_cells[[paste0(cell, "_log2Expression")]]
  y=expdiv_cells[[paste0(cell, "_ExpressionConservation")]]
  abline(lm(y~x), col="red")
  mtext(paste("R2 = ", round(summary(lm(y~x))$r.squared, digits=2), ", Pearson R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
        side=3, line=0.5, cex=CEXstats)

  if (cell == "Bcell"){mtext("Gene Expression conservation", side=2, line=2.5, cex=CEXLAB)}
  mtext(paste0("Mean Expression Level in ", cell), side=1, line=2.5, cex=CEXLAB)
  
}

#### Spearman's rho (Cardoso-Moreira)
expdiv$classRPKM <- cut2(expdiv$Human_MeanRPKM, g=4, include.lowest=T) 
clearboxplot("classRPKM", "CorrelationSpearman", c(-0.5, 1),  "Spearman's rho")

#### Euclidean Similarity (Cardoso-Moreira)
expdiv$EuclideanSimilarity <- 1-expdiv$EuclideanDistance
clearboxplot("classRPKM", "EuclideanSimilarity", c(0.6, 1),  "Euclidean Similarity")

### Correlation between Spearman and Euclidean 
smoothScatter(expdiv$ResidualExpEuclideanSimilarity, expdiv$ResidualExpSpearman, xlab="", ylab="")

R=cor(expdiv$ResidualExpEuclideanSimilarity, expdiv$ResidualExpSpearman,method="pearson")
rho=cor(expdiv$ResidualExpEuclideanSimilarity, expdiv$ResidualExpSpearman, method="spearman")
x=expdiv$ResidualExpEuclideanSimilarity
y=expdiv$ResidualExpSpearman
abline(lm(y~x), col="red")

mtext(paste("R2 = ", round(summary(lm(y~x))$r.squared, digits=2), ", Pearson R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("Residual Euclidean Similarity", side=1, line=2.5, cex=CEXLAB)
mtext("Residual Spearman's rho", side=2, line=2.5, cex=CEXLAB)

######################################################################################################################
############################################ Specificity  # ##########################################################
### Spearman
clearboxplot("classTau", "ResidualExpSpearman", c(-1, 0.6),  "Residual Spearman's rho")

### Euclidean 
clearboxplot("classTau", "ResidualExpEuclideanSimilarity", c(-0.3, 0.15),  "Residual Euclidean Similarity")

### Correlation between measures
smoothScatter(expdiv$CorrectedExpTauEuclideanSimilarity, expdiv$CorrectedExpTauSpearman, xlab="", ylab="")
R=cor(expdiv$CorrectedExpTauEuclideanSimilarity, expdiv$CorrectedExpTauSpearman,method="pearson")
rho=cor(expdiv$CorrectedExpTauEuclideanSimilarity, expdiv$CorrectedExpTauSpearman, method="spearman")
x=expdiv$CorrectedExpTauEuclideanSimilarity
y=expdiv$CorrectedExpTauSpearman
abline(lm(y~x), col="red")

mtext(paste("R2 = ", round(summary(lm(y~x))$r.squared, digits=2), ", Pearson R = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""),
      side=3, line=0.5, cex=CEXstats)

mtext("Corrected Euclidean Similarity", side=1, line=2.5, cex=CEXLAB)
mtext("Corrected Spearman's rho", side=2, line=2.5, cex=CEXLAB)

#dev.off()
