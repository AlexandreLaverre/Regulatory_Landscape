######################################################################################################################
library(Hmisc)
setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")

source("parameters.R") ## pathFinalData are defined based on the user name

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
sp="human"
DivergenceMeasure = "EuclidianSimilarity" # "EuclidianSimilarity" or "CorrelationSpearman"

#########################################################################
# Expression divergence

expdiv=read.table(paste0(pathFinalData, "SupplementaryDataset6/expression_divergence/ExpressionDivergence_CardosoMoreira2019_correlations.txt"), h=T, stringsAsFactors=F, sep="\t")
rownames(expdiv)=expdiv$IDHuman

expdiv$classTau=cut2(expdiv$TauHuman, g=4, include.lowest=T)
expdiv$EuclidianSimilarity = 1-expdiv$ResidualExpressionDivergence

## select genes: protein-coding genes, in the PCHiC data, and in expression stats
annot=gene.annot[[sp]]
genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
genes=intersect(rownames(expdiv), genes)
expdiv_all = expdiv[genes,]


#########################################################################
# Regulatory Divergence

maxdist = "25000"

regland = list()
expdiv = list()
for (enh in enhancer.datasets[[sp]]){
  file = paste("/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset6_regulatory_landscape_evolution/human/",
               enh, "_original_evolution_summary_all_", maxdist, "_0.8.txt", sep="")
  
  regland_enh = read.table(file, h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  
  regland_enh <- regland_enh[which(regland_enh$nb_total >= 5),] 
  common=intersect(rownames(expdiv_all), rownames(regland_enh))
  regland_enh=regland_enh[common,]
  
  regland_enh$ratio_cons_seq = regland_enh$nb_seq_conserv/regland_enh$nb_total
  regland_enh$ratio_cons_synt = ifelse(regland_enh$nb_seq_conserv > 0, regland_enh$nb_synt2M_conserv/regland_enh$nb_seq_conserv, NA)
  regland_enh$ratio_cons_int = ifelse(regland_enh$nb_seq_conserv > 0, regland_enh$nb_contact_conserv/regland_enh$nb_seq_conserv, 0)
  
  regland_enh$class_nb_contact=cut2(regland_enh$nb_total, g=5, include.lowest=T)
  regland_enh$class_cons_seq=cut(regland_enh$ratio_cons_seq, breaks=c(0, 0.10, 0.25, 0.5, 0.75, 1), include.lowest=T)
  regland_enh$class_cons_synt=cut(regland_enh$ratio_cons_synt, breaks=c(0, 0.75, 0.99, 1), include.lowest=T)
  regland_enh$class_cons_int=cut(regland_enh$ratio_cons_int,  breaks=c(0, 0.01, 0.25, 0.5, 0.75, 1), include.lowest=T)
  regland_enh$class_align_score=cut(regland_enh$med_align_score, breaks=seq(0,1,0.1), include.lowest=T)
  
  
  regland[[enh]] = regland_enh
  expdiv[[enh]] = expdiv_all[common,]
}

##################################################################

#pdf(file=paste(pathFigures, "Figure7_Spearman.pdf", sep=""), width=6.85, height=5.5)

m=matrix(rep(NA, 3*10), nrow=3)
m[1,]=c(rep(1, 5), rep(2,5))
m[2,]=c(rep(3, 5), rep(4,5))
m[3,]=c(rep(5, 5), rep(6,5))
layout(m)

##################################################################

## A - classes of expression specificity

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
mtext(xpos, at=xpos, side=1, line=1, cex=0.8)
mtext("expression specificity quartile", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Expression Similarity", side=2, line=2.5, cex=0.9)

mtext("A", side=3, at=0.5, font=2, cex=1.1, line=1.5)


##################################################################

## B - Specificity vs nb enhancers

ylim=c(0, 50)
par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  
  for(class in levels(expdiv[[enh]]$classTau)){
    this.genes=rownames(expdiv[[enh]][which(expdiv[[enh]]$classTau == class),])
    
    x=xpos[class]+smallx[enh]
    
    b=boxplot(regland[[enh]][this.genes, "nb_total"], plot=FALSE)
    med=median(regland[[enh]][this.genes, "nb_total"])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:length(levels(expdiv_all$classTau))-1]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(expdiv_all$classTau))), cex.axis=0.8)
mtext(xpos, at=xpos, side=1, line=1, cex=0.8)
mtext("expression specificity quartile", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("nb. contacted enhancers", side=2, line=2.5, cex=0.9)

legend("topright", legend=label.enhancers[enhancer.datasets[[sp]]], lty=1, col=col.enhancers[enhancer.datasets[[sp]]], bty="o", cex=0.8, inset=c(-0.01, -0.25), seg.len=1,  xpd=NA, box.col="white", bg="white")
mtext("B", side=3, at=0.5, font=2, cex=1.1, line=1.5)

#########################################################################
# C - Divergence vs nb enhancers
if (DivergenceMeasure == "CorrelationSpearman"){ylim=c(0.45, 0.65)}else{ylim=c(1,1.05)}

xlim=c(0.5, length(levels(regland[[enh]]$class_nb_contact))+0.5)

par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  message(enh)
  for(class in levels(regland[[enh]]$class_nb_contact)){
    this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_nb_contact == class),])
    
    xpos=seq(1, length(levels(regland[[enh]]$class_nb_contact)), 1)
    names(xpos) = levels(regland[[enh]]$class_nb_contact)
    x=xpos[class]+smallx[enh]
    
    b=boxplot(expdiv[[enh]][this.genes, DivergenceMeasure], plot=FALSE)
    med=median(expdiv[[enh]][this.genes, DivergenceMeasure])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[enh]]$class_nb_contact))), cex.axis=0.8)
mtext(xpos, at=xpos, side=1, line=1, cex=0.8)
mtext("nb. contacted enhancers quartile", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Expression Similarity", side=2, line=2.5, cex=0.9)

mtext("C", side=3, at=0.5, font=2, cex=1.1, line=1.5)

#########################################################################
# D - Divergence vs enhancers conserved in sequences
if (DivergenceMeasure == "CorrelationSpearman"){ylim=c(0.4, 0.9)}else{ylim=c(0.9,1.05)}

xlim=c(0.5, length(levels(regland[[enh]]$class_cons_seq))+0.5)

par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  message(enh)
  for(class in levels(regland[[enh]]$class_cons_seq)){
    this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_cons_seq == class),])
    
    xpos=seq(1, length(levels(regland[[enh]]$class_cons_seq)), 1)
    names(xpos) = levels(regland[[enh]]$class_cons_seq)
    x=xpos[class]+smallx[enh]
    
    b=boxplot(expdiv[[enh]][this.genes, DivergenceMeasure], plot=FALSE)
    med=median(expdiv[[enh]][this.genes, DivergenceMeasure])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[enh]]$class_cons_seq))), cex.axis=0.8)
mtext(c("<10%", "10-25%", "25-50%", "50-75%", ">75%"), at=xpos, side=1, line=1, cex=0.8)
mtext("conserved enhancer in sequence", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Expression Similarity", side=2, line=2.5, cex=0.9)

mtext("D", side=3, at=0.5, font=2, cex=1.1, line=1.5)

#########################################################################
# E - Divergence vs enhancers conserved in synteny

if (DivergenceMeasure == "CorrelationSpearman"){ylim=c(0.45, 0.6)}else{ylim=c(1,1.05)}

xlim=c(0.5, length(levels(regland[[enh]]$class_cons_synt))+0.5)

par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  message(enh)
  for(class in levels(regland[[enh]]$class_cons_synt)){
    this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_cons_synt == class),])
    
    xpos=seq(1, length(levels(regland[[enh]]$class_cons_synt)), 1)
    names(xpos) = levels(regland[[enh]]$class_cons_synt)
    x=xpos[class]+smallx[enh]
    
    b=boxplot(expdiv[[enh]][this.genes, DivergenceMeasure], plot=FALSE)
    med=median(expdiv[[enh]][this.genes, DivergenceMeasure])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[enh]]$class_cons_synt))), cex.axis=0.8)
mtext(c("<75%", "75-99%", ">99%"), at=xpos, side=1, line=1, cex=0.8)
mtext("conserved enhancer in synteny", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Expression Similarity", side=2, line=2.5, cex=0.9)

mtext("E", side=3, at=0.5, font=2, cex=1.1, line=1.5)

#########################################################################
# F - Divergence vs enhancers conserved in contact

if (DivergenceMeasure == "CorrelationSpearman"){ylim=c(0.5, 0.65)}else{ylim=c(1,1.05)}

xlim=c(0.5, length(levels(regland[[enh]]$class_cons_int))+0.5)

par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  message(enh)
  for(class in levels(regland[[enh]]$class_cons_int)){
    this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_cons_int == class),])
    
    xpos=seq(1, length(levels(regland[[enh]]$class_cons_int)), 1)
    names(xpos) = levels(regland[[enh]]$class_cons_int)
    x=xpos[class]+smallx[enh]
    
    b=boxplot(expdiv[[enh]][this.genes, DivergenceMeasure], plot=FALSE)
    med=median(expdiv[[enh]][this.genes, DivergenceMeasure])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[enh]]$class_cons_int))), cex.axis=0.8)
mtext(c("<1%", "1-25%", "25-50%", "50-75%", ">75%"), at=xpos, side=1, line=1, cex=0.8)
mtext("conserved enhancer in contact", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Expression Similarity", side=2, line=2.5, cex=0.9)

mtext("F", side=3, at=0.5, font=2, cex=1.1, line=1.5)
