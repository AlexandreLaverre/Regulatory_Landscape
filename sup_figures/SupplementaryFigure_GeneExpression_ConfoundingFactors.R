######################################################################################################################

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalData are defined based on the user name

sp="human"

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/Fig6_", sp, "_common_cells.Rdata", sep=""))
load(paste(pathFigures, "RData/Fig6_", sp, "_all_samples_CardosoMoreira.Rdata", sep=""))

if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}
DivergenceMeasure = "CorrelationSpearman" # "EuclidianSimilarity" or "CorrelationSpearman" or SpearmanResidual

cells <- c("Bcell", "ESC", "adipo")
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

######################################################################################################################
##### A - Gene Expression Specificity (Cardoso-Moreira) and gene expression profil similarity (Spearman) ####
xpos=seq(1, length(levels(expdiv_all$classTau)), 1)
names(xpos) = levels(expdiv_all$classTau)

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]

if (DivergenceMeasure == "CorrelationSpearman"){ylim=c(-0.1, 1)}else{ylim=c(0.7, 1.2)}
xlim=c(0.5, length(levels(expdiv_all$classTau))+0.5)

par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(class in levels(expdiv_all$classTau)){
  x=xpos[class]
  boxplot(expdiv_all[which(expdiv_all$classTau == class), DivergenceMeasure], at=x, axes=F, add=T, notch=T, outline=F)
}

abline(v=xpos[1:length(levels(expdiv_all$classTau))-1]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(expdiv_all$classTau))), cex.axis=0.8)
mtext(c("generalist", "", "", "specialist"), at=xpos, side=1, line=1, cex=0.7)
mtext("Expression Specificity", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Expression Similarity (Spearman's coefficient)", side=2, line=2.5, cex=0.9)

mtext("A", side=3, at=0.5, font=2, cex=1.1, line=1.5)

##### B - Gene expression level (Cardoso-Moreira) and gene expression profil similarity (Spearman) ####

##### C - Gene expression level (common cell types) and gene expression level evolution  ####