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

ylab = c("% length covered by repeats", "% length covered by exons", "GC content", "nb. genes within 500kb")
names(ylab) = features

##################################################################################################################
##########################  enhancer proportion according to distance ############################################

for (sp in c("human", "mouse")){
  load(paste(pathFigures, "RData/data.features.coverage.", sp, ".RData", sep=""))
  
  if (sp == "human"){
    pdf.name="SupplementaryFigure2.pdf"
  } else{
    pdf.name="SupplementaryMaterialFigure19.pdf"
  }
  
  pdf(paste(pathFigures, pdf.name, sep=""), width=6.85, height=4.5)
   
  m=matrix(rep(NA, 82), nrow=2)
  m[1,]=c(rep(1:4, each=10), rep(5, 1))
  m[2,]=c(rep(6:9, each=10), rep(10,1))

  layout(m)
    
  nb=1
  
  for (seq in sequences){

    par(mar = c(3, 3.1, 1.5, 0.1)) # bottom, left, top, right
    
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
      axis(side=1, at=c(0,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=0.9)
      axis(side=2, mgp=c(3, 0.65, 0), cex.axis=0.9, las=2)
      
      if (nb == 4){
        legend("topright", legend=c("PCHi-C data", "simulated data"), 
               col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1, xpd=NA, inset=c(0.01, -0.05))
      }
      
      mtext(ylab[feat], side=2, cex=0.7, line=2)

      if(seq=="fragment"){
         mtext("distance to baits (Mb)", side=1, line=1.75, cex=0.7)
      } else{
        mtext("distance to promoters (Mb)", side=1, line=1.75, cex=0.7)
      }
      mtext(letters[nb], side=3, line=0, at=-10.75, font=2, cex=0.9)
      
      nb = nb+1
    }

    par(mar = c(0, 0, 0, 0)) # bottom, left, top, right

    plot.new()
    
    if(seq=="fragment"){
      mtext("restriction fragments", side=2, line=-1, cex=0.7)
    } else{
      mtext("enhancers", side=2, line=-1, cex=0.7)
    }
    
  }
  
  dev.off()
  
}

###########################################################################################################################
