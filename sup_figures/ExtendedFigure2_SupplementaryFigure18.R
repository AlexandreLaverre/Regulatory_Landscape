#################################################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  source("../main_figures/parameters.R") ## paths are defined based on the user name
  load=T
}

##################################################################################################################

features = c("repeat", "all_exon", "GC", "genes")
sequences <- c("fragment", "ENCODE")

ylab = c("% length covered by repeats", "% length covered by exons", "GC content", "number of genes within 500kb")
names(ylab) = features

##################################################################################################################
##########################  enhancer proportion according to distance ############################################

for (sp in c("human", "mouse")){
  load(paste(pathFigures, "RData/data.features.coverage.", sp, ".Rdata", sep=""))
  
  if (sp == "human"){
    pdf.name="ExtendedFigure2.pdf"
  } else{
    pdf.name="SupplementaryFigure18.pdf"
  }
  
  pdf(paste(pathFigures, pdf.name, sep=""), width=4.49, height=10)
  
  par(mai = c(0.5, 0.5, 0.5, 0.2)) # bottom, left, top, right
  par(mfcol=c(length(sequences),2))
  
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
      
      
      class_leg <- c("0.05", "0.5", "1", "1.5", "2")
      axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1)
      axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1, las=2)
      
      if (nb == 1){
        legend("topleft", legend=c("PCHi-C data", "simulated data"), 
               col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1, xpd=NA, inset=c(0.01, -0.15))
      }
      
      mtext(ylab[feat], side=2, cex=0.75, line=2.25)
      
      mtext("distance to promoters (Mb)", side=1, line=2, cex=0.75)
      
      
      if (nb %in% c(2,5,8,11,14)){
        if (seq == "fragment"){
          seq.lab="restriction fragments"
        }else{
          seq.lab=paste(seq, "enhancers")
        }
        mtext(seq.lab, side=3, cex=0.75, line=1)
      }
      
      mtext(letters[nb], side=3, line=1.45, at=-9.75, font=2, cex=1)
      
      nb = nb+1
    }
    
  }
  
  dev.off()
  
}

###########################################################################################################################
