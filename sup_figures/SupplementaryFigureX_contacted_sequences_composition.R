#################################################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
}

source("parameters.R") ## paths are defined based on the user name

if(load){
  sp = "mouse"
  load(paste(pathFigures, "RData/data.features.coverage.", sp, ".Rdata", sep=""))
  
  features = c("GC", "all_exon", "repeat")
  sequences <- c("fragment", "ENCODE")
  ylab = c("% length covered\n by repeated elements", "GC rate", "% length covered\n by exons")
}

###########################################################################################################################
##########################  featancer proportion according to distance ############################################

pdf(paste(pathFigures, "SupplementaryFigureX_contacted_sequence_composition_", sp, ".pdf", sep=""), width=6.85, height=5.5)

par(mai = c(0.5, 0.5, 0.5, 0.2)) # bottom, left, top, right

par(mfrow=c(length(sequences),3))

nb=1

for (seq in sequences){
  for (feat in features){
    
    ymax=max(c(feat_prop_dist[[seq]][["obs"]][[paste0(feat,"_conflow")]], feat_prop_dist[[seq]][["obs"]][[paste0(feat,"_confup")]],  feat_prop_dist[[seq]][["simul"]][[paste0(feat,"_conflow")]], feat_prop_dist[[seq]][["simul"]][[paste0(feat,"_confup")]]), na.rm =T)
    ymin=min(c(feat_prop_dist[[seq]][["obs"]][[paste0(feat,"_conflow")]], feat_prop_dist[[seq]][["obs"]][[paste0(feat,"_confup")]],  feat_prop_dist[[seq]][["simul"]][[paste0(feat,"_conflow")]], feat_prop_dist[[seq]][["simul"]][[paste0(feat,"_confup")]]), na.rm =T)
    
    ymax=ymax*1.02
    
    
    plot(feat_prop_dist[[seq]][["obs"]][[feat]], col=dataset.colors["Original"], main="", type="n", xlab="", ylab="",  axes=F, ylim=c(ymin,ymax))
    
    lines(feat_prop_dist[[seq]][["obs"]][[feat]], col=dataset.colors["Original"])
    lines(feat_prop_dist[[seq]][["simul"]][[feat]], col=dataset.colors["Simulated"])
    
    
    xpos=1:length(feat_prop_dist[[seq]][["obs"]][[feat]])
    
    segments(xpos, feat_prop_dist[[seq]][["obs"]][[paste0(feat,"_conflow")]], xpos, feat_prop_dist[[seq]][["obs"]][[paste0(feat,"_confup")]], col=dataset.colors["Original"])
    segments(xpos, feat_prop_dist[[seq]][["simul"]][[paste0(feat,"_conflow")]], xpos, feat_prop_dist[[seq]][["simul"]][[paste0(feat,"_confup")]], col=dataset.colors["Simulated"])
    
    
    class_leg <- c("0", "0.5", "1", "1.5", "2")
    axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
    axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)
    
    if (nb == 1){
      legend("topright", legend=c("PCHi-C data", "simulated data"), 
             col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1.1, xpd=NA)}
    
    mtext(ylab[nb%%3+1], side=2, cex=0.85, line=2)
    
    mtext("distance to promoters (Mb)", side=1, line=2, cex=0.85)

    
    if (nb %in% c(2,5,8,11,14)){
      if (seq == "fragment"){seq.lab="restriction fragments"}else{seq.lab=seq}
      mtext(seq.lab, side=3, cex=0.8, line=1)
    }
    
    
    mtext(letters[nb], side=3, line=1.45, at=-5.75, font=2, cex=1.2)
    
    nb = nb+1
  }
  
  }
  

dev.off()

