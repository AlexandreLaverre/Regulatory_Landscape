##########################################################################

objects=ls()

if(!"pathFigures"%in%objects){

  source("../main_figures/parameters.R")
 
## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

}

##########################################################################

  
for(ref_sp in c( "human", "mouse")){
  
  if(ref_sp=="human"){
    target_species=c("macaque", "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken") ## other species already done
    pdfname=paste(pathFigures, "GenomeResearch_Figures/SupplementaryMaterialFigure10.pdf", sep="")
  }
  
  if(ref_sp=="mouse"){
    target_species=c("rat", "rabbit", "human", "macaque", "dog", "cow", "elephant", "opossum", "chicken")
    pdfname=paste(pathFigures, "GenomeResearch_Figures/SupplementaryMaterialFigure11.pdf", sep="")
  }
  
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", ref_sp, ".RData", sep=""))
  
  selenh="ENCODE"
  
  align_enhancer_obs=list_align_enh[[selenh]][["enh_align_obs"]]
  align_enhancer_simul=list_align_enh[[selenh]][["enh_align_simul"]]
  
##########################################################################
  
  pdf(pdfname, width=6.85, height=10)
  
  ## layout
  
  m=matrix(1:8, nrow=4)
  layout(m)
  
###########################################################################
  
  labels=toupper(letters[1:length(target_species)])
  names(labels)=target_species
  
  par(mar=c(4.1, 4.5, 2.1, 1.5))
  
  for(other_sp in target_species){
  
    nbclasses=length(levels(align_enhancer_obs$dist_class))
    xpos=1:nbclasses
    
    xlim=c(-0.5, max(xpos)+1)
    
  ## axis position
    class_leg <- c("0", "0.5", "1", "1.5", "2")
    xax=seq(from=0, to=max(xpos)+1, by=10)
    
    load(paste(pathFigures, "RData/data.bootstrap.enhancer.conservation.distance.",ref_sp,".",other_sp,".RData",sep=""))
  
    ylim=range(c(ci.low.obs, ci.high.obs, ci.low.simul, ci.high.simul))
    
    dy=diff(ylim)/20
    ylim=ylim+c(-dy, dy)
    
    plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")
    
    ## lines(xpos, mean.val.obs, col=dataset.colors["Original"])
    points(xpos, mean.val.obs, col=dataset.colors["Original"], pch=20)
    segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])
    
    ##lines(xpos, mean.val.simul, col=dataset.colors["Simulated"])
    points(xpos, mean.val.simul, col=dataset.colors["Simulated"], pch=20)
    segments(xpos, ci.low.simul, xpos, ci.high.simul, col=dataset.colors["Simulated"])
    
    axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
    mtext("distance to promoters (Mb)", side=1, line=2.2, cex=0.8)
    
    axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
    mtext("% aligned sequence", side=2, line=3, cex=0.8)
    
    mtext(paste(ref_sp, " vs. ", other_sp,sep=""), side=3, cex=0.8)
    
    mtext(labels[other_sp], side=3, line=1, at=-7.5, font=2, cex=1.2)
  }
  
###########################################################################
  
  ## legend
  
  legend("topleft", legend = c("PCHi-C data", "simulated data"), col=dataset.colors, lty=1, bty='n',  inset=c(0, 0.01), xpd=NA, cex=1.25)
    
###########################################################################
  
dev.off()

###########################################################################

}
