######################################################################################################################
library(Hmisc)
#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  
  source("parameters.R")
}

##############################################################################
if(load){
  sp="human"
  load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
  load(paste(pathFigures, "RData/Fig6_", sp, "_SomaticOrgans_CardosoMoreira.Rdata", sep=""))
  if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}

  load=FALSE
}

################################################################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure34.pdf", sep=""), width = 8.5)

par(mfrow=c(2,3))

################################################################################################################################
############################## Cardoso-Moreira  - Euclidean Similarity ##########################################################

data = "original"
nb = 1
expdiv$EuclideanSimilarity = 1-expdiv$EuclideanDistance

for (DivergenceMeasure in c("EuclideanSimilarity", "CorrectedEuclideanSimilarity")){
  
  if (DivergenceMeasure == "EuclideanSimilarity"){xlab="Euclidean Similarity"; ylim=c(0.9, 0.95)
  }else{xlab="Residual Euclidean Similarity"; ylim=c(0, 0.02)}
  
  #### A - Gene expression profil similarity and enhancers conserved in sequences ####
  par(mai = c(0.5, 0.5, 0.5, 0.1)) # bottom, left, top, right
  xpos=seq(1, length(levels(expdiv$classTau)), 1)
  names(xpos) = levels(expdiv$classTau)
  
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  xlim=c(0.5, length(levels(regland[[data]][["ENCODE"]]$class_align_score))+0.5)
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(class in levels(regland[[data]][[enh]]$class_align_score)){
      this.genes=rownames(regland[[data]][[enh]][which(regland[[data]][[enh]]$class_align_score == class),])
      
      xpos=seq(1, length(levels(regland[[data]][[enh]]$class_align_score)), 1)
      names(xpos) = levels(regland[[data]][[enh]]$class_align_score)
      x=xpos[class]+smallx[enh]
      
      b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=FALSE)
      med=median(expdiv[this.genes, DivergenceMeasure])
      ci=as.numeric(b$conf)
      
      points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[data]][[enh]]$class_align_score))), cex.axis=0.8)
  mtext(xpos, at=xpos, side=1, line=1, cex=0.8)
  mtext("Alignment score", side=1, line=2.5, cex=0.9)
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
  mtext(xlab, side=2, line=2.5, cex=0.9)
  
  mtext(letters[nb], side=3, at=0.45, font=2, cex=1.1, line=0.5)
  legend("topleft", legend=label.enhancers[enhancer.datasets[[sp]]], pch=20,
         col=col.enhancers[enhancer.datasets[[sp]]], cex=1,
         bty="o", box.col="white", bg="white",  inset=c(0.01, 0.01))
  
  nb = nb +1
  
  #### B - Gene expression profil similarity and enhancers conserved in synteny ####
  par(mai = c(0.5, 0.2, 0.5, 0)) # bottom, left, top, right
  
  xlim=c(0.5, length(levels(regland[[data]][[enh]]$class_cons_synt))+0.5)
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(class in levels(regland[[data]][[enh]]$class_cons_synt)){
      this.genes=rownames(regland[[data]][[enh]][which(regland[[data]][[enh]]$class_cons_synt == class),])
      
      xpos=seq(1, length(levels(regland[[data]][[enh]]$class_cons_synt)), 1)
      names(xpos) = levels(regland[[data]][[enh]]$class_cons_synt)
      x=xpos[class]+smallx[enh]
      
      b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=FALSE)
      med=median(expdiv[this.genes, DivergenceMeasure])
      ci=as.numeric(b$conf)
      
      points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:2]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[data]][[enh]]$class_cons_synt))), cex.axis=0.8)
  mtext(c("<75%", "75-99%", ">99%"), at=xpos, side=1, line=1, cex=0.8)
  mtext("Synteny Conservation", side=1, line=2.5, cex=0.9)
  mtext(letters[nb], side=3, at=0.45, font=2, cex=1.1, line=0.5)
  
  #### C - Gene expression profil similarity and enhancers conserved in contact ####
  xlim=c(0.5, length(levels(regland[[data]][[enh]]$class_cons_int))+0.5)
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(class in levels(regland[[data]][[enh]]$class_cons_int)){
      this.genes=rownames(regland[[data]][[enh]][which(regland[[data]][[enh]]$class_cons_int == class),])
      
      xpos=seq(1, length(levels(regland[[data]][[enh]]$class_cons_int)), 1)
      names(xpos) = levels(regland[[data]][[enh]]$class_cons_int)
      x=xpos[class]+smallx[enh]
      
      b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=FALSE)
      med=median(expdiv[this.genes, DivergenceMeasure])
      ci=as.numeric(b$conf)
      
      points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[data]][[enh]]$class_cons_int))), cex.axis=0.8)
  mtext(c("<1%", "1-25%", "25-50%", "50-75%", ">75%"), at=xpos, side=1, line=1, cex=0.8)
  mtext("Contact Conservation", side=1, line=2.5, cex=0.9)
  mtext(letters[nb], side=3, at=0.45, font=2, cex=1.1, line=0.5)
  
}

dev.off()
