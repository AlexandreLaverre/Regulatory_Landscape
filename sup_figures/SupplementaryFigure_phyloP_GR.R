#########################################################################################################################

objects=ls()

if(!"pathFigures"%in%objects){
  
  source("parameters.R")
  
  library(ape)
  library(vioplot)
  
  set.seed(19)
  
  load=T
  prepare=T
}

#########################################################################################################################

if(load){
  ref_sp = "human"
  target_sp = "mouse"
  
  enhancers = enhancer.datasets[[ref_sp]]
  
  selenh="ENCODE"
  
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.and.phyloP.", ref_sp, ".RData", sep=""))
  
  
  ## enhancer statistics
  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
  
  enh.stats.obs=enhancer.statistics[[ref_sp]][[selenh]][["original"]]
  enh.stats.simul=enhancer.statistics[[ref_sp]][[selenh]][["simulated"]]
  
  ## enhancer alignment
  align_enhancer_obs=list_align_enh[[selenh]][["enh_align_obs"]]
  align_enhancer_sim=list_align_enh[[selenh]][["enh_align_simul"]]
  
  ## fragment statistics
  load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))
  
  frag.stats.obs=fragment.statistics[[ref_sp]][["original"]]
  frag.stats.simul=fragment.statistics[[ref_sp]][["simulated"]]
  
  load=F
}

#########################################################################################################################

if(prepare){
  
  ## add nb genes 
  
  frag_align_obs$nb_genes_500kb=frag.stats.obs[frag_align_obs$ID, "nb_genes_500kb"]
  frag_align_simul$nb_genes_500kb=frag.stats.simul[frag_align_simul$ID, "nb_genes_500kb"]
  
  align_enhancer_obs$nb_genes_500kb=enh.stats.obs[align_enhancer_obs$ID, "nb_genes_500kb"]
  align_enhancer_sim$nb_genes_500kb=enh.stats.simul[align_enhancer_sim$ID, "nb_genes_500kb"]
  
  ## add nb genes class
  
  frag_align_obs$class_genes_500kb=cut(frag_align_obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag_align_obs$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))
  
  frag_align_simul$class_genes_500kb=cut(frag_align_simul$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag_align_simul$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))
  
  align_enhancer_obs$class_genes_500kb=cut(align_enhancer_obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(align_enhancer_obs$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))
  
  align_enhancer_sim$class_genes_500kb=cut(align_enhancer_sim$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(align_enhancer_sim$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))
  
  ## add repeats
  
  frag.stats.obs$pcrepeat=100*frag.stats.obs$repeat_bp/frag.stats.obs$length
  frag.stats.simul$pcrepeat=100*frag.stats.simul$repeat_bp/frag.stats.simul$length
  
  frag_align_obs$pcrepeat=frag.stats.obs[frag_align_obs$ID, "pcrepeat"]
  frag_align_simul$pcrepeat=frag.stats.simul[frag_align_simul$ID, "pcrepeat"]
  
  enh.stats.obs$pcrepeat=100*enh.stats.obs$repeat_bp/enh.stats.obs$length
  enh.stats.simul$pcrepeat=100*enh.stats.simul$repeat_bp/enh.stats.simul$length
  
  align_enhancer_obs$pcrepeat=enh.stats.obs[align_enhancer_obs$ID, "pcrepeat"]
  align_enhancer_sim$pcrepeat=enh.stats.simul[align_enhancer_sim$ID, "pcrepeat"]
  
  prepare=FALSE
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#########################################################################################################################

pdf(paste(pathFigures, "GenomeResearch_Figures/Supplemental_Fig_phyloP.pdf", sep=""), width=6.85, height=7)

par(mfrow=c(2,2))
par(mai = c(0.5, 0.1, 0.3, 0.1)) #bottom, left, top and right

######################## Conserved sequence human vs distance to promoters  ########################

par(mai = c(0.8, 0.6, 0.2, 0.2)) #bottom, left, top and right
par(mar=c(4.1, 4.5, 2, 1.5))

nbclasses=length(levels(frag_align_obs$dist_class))
xpos=1:nbclasses

xlim=c(-0.5, max(xpos)+1)

## axis position
class_leg <- c("0", "0.5", "1", "1.5", "2")
xax=seq(from=0, to=max(xpos)+1, by=10)

labels=c("A", "B")
names(labels)=c("restriction fragments", "enhancers")

for (ref_sp in c("human", "mouse")){
  for(type in c("fragments", "enhancer")){
    
    load(paste(pathFigures, "RData/data.bootstrap.",type,".conservation.distance.",ref_sp,".phyloP.RData", sep=""))
    
    ylim=range(c(ci.low.obs, ci.high.obs, ci.low.simul, ci.high.simul))
    
    dy=diff(ylim)/20
    ylim=ylim+c(-dy, dy)
    
    plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")
    
    points(xpos, mean.val.obs, col=dataset.colors["Original"], pch=20)
    segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])
    
    points(xpos, mean.val.simul, col=dataset.colors["Simulated"], pch=20)
    segments(xpos, ci.low.simul, xpos, ci.high.simul, col=dataset.colors["Simulated"])
    
    axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
    
    if(type=="enhancer"){
      mtext("distance to promoters (Mb)", side=1, line=2.2, cex=0.8)
    }
    
    if(type=="fragments"){
      mtext("distance to baits (Mb)", side=1, line=2.2, cex=0.8)
    }
    
    axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
    mtext("phyloP score", side=2, line=3, cex=0.8)
    
    mtext(paste(ref_sp,", ", type, sep=""), side=3, cex=0.8, line=1)
    
    mtext(labels[type], side=3, line=1, at=-7.75, font=2, cex=1.1)
  }
}



#######################################################################################################
# 
# ## sequence conservation for gene classes
# 
# par(mai = c(0.8, 0.6, 0.2, 0.2)) #bottom, left, top and right
# 
# par(mar=c(5.1, 4.5, 1, 1.5))
# 
# nbclasses=length(levels(frag_align_obs$class_genes_500kb))
# xpos=1:nbclasses
# 
# xlim=c(0.5, max(xpos)+0.5)
# 
# smallx=c(-0.1, 0.1)
# names(smallx)=c("obs", "sim")
# 
# ## axis position
# 
# xax=xpos
# class_leg=levels(frag_align_obs$class_genes_500kb)
# 
# labels=c("F", "G")
# names(labels)=c("restriction fragments", "enhancers")
# 
# 
# for(type in c("restriction fragments", "enhancers")){
#   
#   load(paste(pathFigures, "RData/data.bootstrap.conservation.gene.density.",type,".",ref_sp,".RData", sep=""))
#   
#   ## plot
#   
#   ylim=range(c(ci.low.obs, ci.high.obs, ci.low.simul, ci.high.simul, ci.high.simul.norep, ci.high.obs.norep))
#   
#   dy=diff(ylim)/20
#   ylim=ylim+c(-dy, dy)
#   
#   plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")
#   
#   points(xpos+smallx["obs"], mean.val.obs, col=dataset.colors["Original"], pch=20, cex=1.35)
#   segments(xpos+smallx["obs"], ci.low.obs, xpos+smallx["obs"], ci.high.obs, col=dataset.colors["Original"])
#   
#   points(xpos+smallx["sim"], mean.val.simul, col=dataset.colors["Simulated"], pch=20, cex=1.35)
#   segments(xpos+smallx["sim"], ci.low.simul, xpos+smallx["sim"], ci.high.simul, col=dataset.colors["Simulated"])
#   
#   points(xpos+smallx["obs"], mean.val.obs.norep, bg="white", col=dataset.colors["Original"], pch=21)
#   segments(xpos+smallx["obs"], ci.low.obs.norep, xpos+smallx["obs"], ci.high.obs.norep, col=dataset.colors["Original"])
#   
#   points(xpos+smallx["sim"], mean.val.simul.norep, bg="white", col=dataset.colors["Simulated"], pch=21)
#   segments(xpos+smallx["sim"], ci.low.simul.norep, xpos+smallx["sim"], ci.high.simul.norep, col=dataset.colors["Simulated"])
#   
#   ## axes
#   axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1, las=2)
#   mtext("number of genes within 500kb", side=1, line=3.85, cex=0.8)
#   
#   axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
#   mtext("% aligned sequence", side=2, line=3, cex=0.8)
#   
#   abline(v=sort(xpos)[-length(xpos)]+0.5, lty=2, col="gray40")
#   
#   mtext(labels[type], side=3, line=1, at=-0.75, font=2, cex=1.2)
#   
#   if(type=="enhancers"){
#     legend("topright", box.col="white", bg="white", pch=21, pt.bg=c("black", "white"), legend=c("all data", "without repeats"),xpd=NA, inset=c(0.01, -0.05), cex=1.1)
#   }
#   
# }
# 

#######################################################################################################

dev.off()

#######################################################################################################
