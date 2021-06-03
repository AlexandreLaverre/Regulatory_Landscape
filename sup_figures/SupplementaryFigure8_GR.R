###################################################################################

library(data.table)

###################################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")

  pathSequenceConservation=paste(pathFinalData, "SupplementaryDataset7/", sep="")
}

###################################################################################

if(load){
  load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))
  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

  enh="ENCODE"
  
  load=FALSE
}

###################################################################################

for(ref in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref)
 
  frag.obs=fragment.statistics[[ref]][["original"]]
  frag.sim=fragment.statistics[[ref]][["simulated"]]
  
  frag.cons=fread(paste(pathSequenceConservation, ref, "/sequence_conservation/restriction_fragments/AlignmentStatistics_Excluding_Exons_", ref,"2", tg,".txt", sep=""), h=T, stringsAsFactors=F)
  class(frag.cons)="data.frame"
  
  rownames(frag.cons)=frag.cons[,paste("ID", ref, sep=".")]

  ## enhancers

  enh.obs=enhancer.statistics[[ref]][[enh]][["original"]]
  enh.sim=enhancer.statistics[[ref]][[enh]][["simulated"]]
  
  enh.cons=fread(paste(pathSequenceConservation, ref, "/sequence_conservation/enhancers/", enh,"/AlignmentStatistics_Excluding_Exons_", ref,"2", tg,".txt", sep=""), h=T, stringsAsFactors=F)
  class(enh.cons)="data.frame"

  rownames(enh.cons)=enh.cons[,paste("ID", ref, sep=".")]
  
 
  ## is lifted
  frag.obs$is_lifted=rep("no",nrow(frag.obs))
  frag.obs$is_lifted[which(rownames(frag.obs)%in%rownames(frag.cons))]="yes"
  
  frag.sim$is_lifted=rep("no",nrow(frag.sim))
  frag.sim$is_lifted[which(rownames(frag.sim)%in%rownames(frag.cons))]="yes"

  p.lifted.frag.obs=prop.test(length(which(frag.obs$is_lifted=="yes")), nrow(frag.obs))
  p.lifted.frag.sim=prop.test(length(which(frag.sim$is_lifted=="yes")), nrow(frag.sim))

  ## pc repeats

  frag.obs$pcrep=100*frag.obs$repeat_bp/frag.obs$length
  frag.sim$pcrep=100*frag.sim$repeat_bp/frag.sim$length

  ## pc ungapped and pc identical

  frag.cons$pcungapped=100*frag.cons$FilteredUngappedLength/frag.cons$FilteredAlignmentLength
  frag.cons$pcidentical=100*frag.cons$FilteredIdenticalLength/frag.cons$FilteredUngappedLength

  frag.cons$pcungapped[which(frag.cons$FilteredAlignmentLength<10)]=NA
  frag.cons$pcidentical[which(frag.cons$FilteredUngappedLength<10)]=NA

  ## same for enhancers

  
  ## is lifted
  enh.obs$is_lifted=rep("no",nrow(enh.obs))
  enh.obs$is_lifted[which(rownames(enh.obs)%in%rownames(enh.cons))]="yes"
  
  enh.sim$is_lifted=rep("no",nrow(enh.sim))
  enh.sim$is_lifted[which(rownames(enh.sim)%in%rownames(enh.cons))]="yes"

  p.lifted.enh.obs=prop.test(length(which(enh.obs$is_lifted=="yes")), nrow(enh.obs))
  p.lifted.enh.sim=prop.test(length(which(enh.sim$is_lifted=="yes")), nrow(enh.sim))

  ## pc repeats

  enh.obs$pcrep=100*enh.obs$repeat_bp/enh.obs$length
  enh.sim$pcrep=100*enh.sim$repeat_bp/enh.sim$length

  ## pc ungapped and pc identical

  enh.cons$pcungapped=100*enh.cons$FilteredUngappedLength/enh.cons$FilteredAlignmentLength
  enh.cons$pcidentical=100*enh.cons$FilteredIdenticalLength/enh.cons$FilteredUngappedLength

  enh.cons$pcungapped[which(enh.cons$FilteredAlignmentLength<10)]=NA
  enh.cons$pcidentical[which(enh.cons$FilteredUngappedLength<10)]=NA


  ## actual plot

  if(ref=="human"){
    pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_S7.pdf", sep=""), width=6.85, height=6.5)
  }

  if(ref=="mouse"){
    pdf(paste(pathFigures, "GenomeResearch_Figures/SupplementaryMaterialFigure24.pdf", sep=""), width=6.85, height=6.5)
  }
  
  m=matrix(rep(NA, 2*19), nrow=2)
  m[1,]=c(rep(1,4), rep(2,6), rep(3,4), rep(4,4), rep(5,1))
  m[2,]=c(rep(6,4), rep(7,6), rep(8,4), rep(9,4), rep(10, 1))
  
  layout(m)
  
  ###################################################################################
  
  ## percentage lifted fragments
  
  ylim=c(0, 100*max(c(p.lifted.frag.obs$conf.int,p.lifted.frag.sim$conf.int)))
  ylim[2]=ylim[2]+diff(ylim)/5
  xlim=c(0.25, 2.75)
  
  par(mar=c(6.1, 4.5, 2.5, 1.5))
  
  b=barplot(100*c(p.lifted.frag.obs$estimate, p.lifted.frag.sim$estimate),xlim=xlim, ylim=ylim, axes=F, xlab="", ylab="", space=0.5, names=rep("",2), col=dataset.colors, border=dataset.colors)
  segments(b[1], 100*p.lifted.frag.obs$conf.int[1], b[1], 100*p.lifted.frag.obs$conf.int[2])
  
  segments(b[2], 100*p.lifted.frag.sim$conf.int[1], b[2], 100*p.lifted.frag.sim$conf.int[2])
  
  axis(side=1, cex.axis=0.95, at=b, labels=rep("", 2), mgp=c(3, 0.5,0))
  mtext(c("PCHi-C", "simulated"), side=1, line=0.75, las=2, cex=0.75, at=b)
  
  axis(side=2, cex.axis=0.95, mgp=c(3, 0.75,0))
  mtext("% successful liftOver", side=2, line=2.5, cex=0.75)
  
  mtext("A", side=3, line=1, font=2, cex=1, at=-1.5)
  
 ###################################################################################
  ##, repetitive elements
  
  ylim=c(0, 100)
  xlim=c(0.75, 3)
  
  par(mar=c(5.1, 4.5, 2.5, 1))
  
  plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
  
  boxplot(frag.obs$pcrep[which(frag.obs$is_lifted=="yes")], at=1, col="white", border=dataset.colors["Original"], boxwex=0.65, add=T, axes=F, notch=T, outline=F)
  boxplot(frag.sim$pcrep[which(frag.sim$is_lifted=="yes")], at=1.5, col="white", border=dataset.colors["Simulated"], boxwex=0.65, add=T, axes=F, notch=T, outline=F)
  
  boxplot(frag.obs$pcrep[which(frag.obs$is_lifted=="no")], at=2.25, col="white", border=dataset.colors["Original"], boxwex=0.65, add=T, axes=F, notch=T, outline=F)
  boxplot(frag.sim$pcrep[which(frag.sim$is_lifted=="no")], at=2.75, col="white", border=dataset.colors["Simulated"], boxwex=0.65, add=T, axes=F, notch=T, outline=F)
  
  
  axis(side=1, cex.axis=0.95, at=c(1.25, 2.75), labels=rep("", 2), mgp=c(3, 0.5,0))
  mtext(c("lifted", "not lifted"), side=1, at=c(1.25, 2.75), line=0.75, cex=0.75, las=2)
  
  axis(side=2, cex.axis=0.95, mgp=c(3, 0.75,0))
  mtext("% length covered by repeats", side=2, line=2.5, cex=0.75)
  
  mtext("B", side=3, line=1, font=2, cex=1, at=0)
  
 ###################################################################################

## pc ungapped for lifted fragments

  ylim=c(0, 100)
  xlim=c(0.5, 2.5)
  
  par(mar=c(6.1, 3.5, 2.5, 1))
  
  plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
  
  boxplot(frag.cons$pcungapped[which(rownames(frag.cons)%in%rownames(frag.obs))],at=1, col="white", border=dataset.colors["Original"], boxwex=0.95, add=T, axes=F, notch=T, outline=F)
  boxplot(frag.cons$pcungapped[which(rownames(frag.cons)%in%rownames(frag.sim))],at=2, col="white", border=dataset.colors["Simulated"], boxwex=0.95, add=T, axes=F, notch=T, outline=F)
  
  
  axis(side=1, cex.axis=0.95, at=c(1, 2), labels=rep("", 2), mgp=c(3, 0.5,0))
  mtext(c("PCHi-C", "simulated"), side=1, line=0.75, las=2, cex=0.75, at=c(1,2))
  
  axis(side=2, cex.axis=0.95, mgp=c(3, 0.75,0))
  mtext("% aligned sequence", side=2, line=2.5, cex=0.75)
  
  mtext("C", side=3, line=1, font=2, cex=1, at=-0.5)
  
###################################################################################
  
  ## pc identical for lifted fragments
  
  ylim=c(61, 77)
  xlim=c(0.5, 2.5)
  
  par(mar=c(6.1, 3.5, 2.5, 1))
  
  plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
  
  boxplot(frag.cons$pcidentical[which(rownames(frag.cons)%in%rownames(frag.obs))],at=1, col="white", border=dataset.colors["Original"], boxwex=0.95, add=T, axes=F, notch=T, outline=F)
  boxplot(frag.cons$pcidentical[which(rownames(frag.cons)%in%rownames(frag.sim))],at=2, col="white", border=dataset.colors["Simulated"], boxwex=0.95, add=T, axes=F, notch=T, outline=F)
  

  axis(side=1, cex.axis=0.95, at=c(1, 2), labels=rep("", 2), mgp=c(3, 0.5,0))
  mtext(c("PCHi-C", "simulated"), side=1, line=0.75, las=2, cex=0.75, at=c(1,2))
  
  axis(side=2, cex.axis=0.95, mgp=c(3, 0.75,0))
  mtext("% identical sequence", side=2, line=2.5, cex=0.75)
  
  mtext("D", side=3, line=1, font=2, cex=1, at=-0.5)
  
###################################################################################
  
## legend
  par(mar=c(6.1,0,2.5,0))
  plot.new()
  mtext("restriction fragments", side=2, line=-1, cex=0.75)
  
###################################################################################
###################################################################################
  
  ## percentage lifted enhancers
  
  ylim=c(0, 100*max(c(p.lifted.enh.obs$conf.int,p.lifted.enh.sim$conf.int)))
  ylim[2]=ylim[2]+diff(ylim)/5
  xlim=c(0.25, 2.75)
  
  par(mar=c(6.1, 4.5, 2.5, 1.5))
  
  b=barplot(100*c(p.lifted.enh.obs$estimate, p.lifted.enh.sim$estimate),xlim=xlim, ylim=ylim, axes=F, xlab="", ylab="", space=0.5, names=rep("",2), col=dataset.colors, border=dataset.colors)
  segments(b[1], 100*p.lifted.enh.obs$conf.int[1], b[1], 100*p.lifted.enh.obs$conf.int[2])
  
  segments(b[2], 100*p.lifted.enh.sim$conf.int[1], b[2], 100*p.lifted.enh.sim$conf.int[2])
  
  axis(side=1, cex.axis=0.95, at=b, labels=rep("", 2), mgp=c(3, 0.5,0))
  mtext(c("PCHi-C", "simulated"), side=1, line=0.75, las=2, cex=0.75, at=b)
  
  axis(side=2, cex.axis=0.95, mgp=c(3, 0.75,0))
  mtext("% successful liftOver", side=2, line=2.5, cex=0.75)
  
  mtext("E", side=3, line=1, font=2, cex=1, at=-1.5)
  
###################################################################################
  ##, repetitive elements
  
  ylim=c(0, 100)
  xlim=c(0.75, 3.25)
  
  par(mar=c(5.1, 4.5, 2.5, 1))
  
  plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
  
  boxplot(enh.obs$pcrep[which(enh.obs$is_lifted=="yes")], at=1, col="white", border=dataset.colors["Original"], boxwex=0.65, add=T, axes=F, notch=T, outline=F)
  boxplot(enh.sim$pcrep[which(enh.sim$is_lifted=="yes")], at=1.5, col="white", border=dataset.colors["Simulated"], boxwex=0.65, add=T, axes=F, notch=T, outline=F)
  
  boxplot(enh.obs$pcrep[which(enh.obs$is_lifted=="no")], at=2.5, col="white", border=dataset.colors["Original"], boxwex=0.65, add=T, axes=F, notch=T, outline=F)
  boxplot(enh.sim$pcrep[which(enh.sim$is_lifted=="no")], at=3, col="white", border=dataset.colors["Simulated"], boxwex=0.65, add=T, axes=F, notch=T, outline=F)
  
  
  axis(side=1, cex.axis=0.95, at=c(1.25, 2.75), labels=rep("", 2), mgp=c(3, 0.5,0))
  mtext(c("lifted", "not lifted"), side=1, at=c(1.25, 2.75), line=0.75, cex=0.75, las=2)
  
  axis(side=2, cex.axis=0.95, mgp=c(3, 0.75,0))
  mtext("% length covered by repeats", side=2, line=2.5, cex=0.75)
  
  
  mtext("F", side=3, line=1, font=2, cex=1, at=0)
###################################################################################
  
  ## pc ungapped for lifted enhancers
  
  ylim=c(0, 100)
  xlim=c(0.5, 2.5)
  
  par(mar=c(6.1, 3.5, 2.5, 1))
  
  plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
  
  boxplot(enh.cons$pcungapped[which(rownames(enh.cons)%in%rownames(enh.obs))],at=1, col="white", border=dataset.colors["Original"], boxwex=0.95, add=T, axes=F, notch=T, outline=F)
  boxplot(enh.cons$pcungapped[which(rownames(enh.cons)%in%rownames(enh.sim))],at=2, col="white", border=dataset.colors["Simulated"], boxwex=0.95, add=T, axes=F, notch=T, outline=F)
  
  
  axis(side=1, cex.axis=0.95, at=c(1, 2), labels=rep("", 2), mgp=c(3, 0.5,0))
  mtext(c("PCHi-C", "simulated"), side=1, line=0.75, las=2, cex=0.75, at=c(1,2))
  
  axis(side=2, cex.axis=0.95, mgp=c(3, 0.75,0))
  mtext("% aligned sequence", side=2, line=2.5, cex=0.75)
  
  
  mtext("G", side=3, line=1, font=2, cex=1, at=-0.5)
  
###################################################################################
  
  ## pc identical for lifted enhancers
  
  ylim=c(50, 90)
  xlim=c(0.5, 2.5)
  
  par(mar=c(6.1, 3.5, 2.5, 1))
  
  plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)
  
  boxplot(enh.cons$pcidentical[which(rownames(enh.cons)%in%rownames(enh.obs))],at=1, col="white", border=dataset.colors["Original"], boxwex=0.95, add=T, axes=F, notch=T, outline=F)
  boxplot(enh.cons$pcidentical[which(rownames(enh.cons)%in%rownames(enh.sim))],at=2, col="white", border=dataset.colors["Simulated"], boxwex=0.95, add=T, axes=F, notch=T, outline=F)
  
  
  axis(side=1, cex.axis=0.95, at=c(1, 2), labels=rep("", 2), mgp=c(3, 0.5,0))
  mtext(c("PCHi-C", "simulated"), side=1, line=0.75, las=2, cex=0.75, at=c(1,2))
  
  axis(side=2, cex.axis=0.95, mgp=c(3, 0.75,0))
  mtext("% identical sequence", side=2, line=2.5, cex=0.75)
  
  
  mtext("H", side=3, line=1, font=2, cex=1, at=-0.5)
    
###################################################################################
  
  ## legend
  
  par(mar=c(6.1,0,2.5,0))
  plot.new()
  
  mtext("enhancers", side=2, line=-1, cex=0.75)
  
 #################################################################################
  
  dev.off()
 ##################################################################################
  
}

##################################################################################
