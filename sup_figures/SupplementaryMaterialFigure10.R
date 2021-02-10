#################################################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  prepare=T
}

source("../main_figures/parameters.R") ## paths are defined based on the user name

###########################################################################################################################
##########################  Enhancer proportion according to distance ############################################

pdf(paste(pathFigures, "SupplementaryMaterialFigure10.pdf", sep=""), width=6.85, height=7.5)

par(mai = c(0.5, 0.5, 0.1, 0.2)) # bottom, left, top, right

par(mfrow=c(3,2))

nb=1

for (sp in c("human", "mouse")){
  
  load(paste(pathFigures, "RData/data.enhancer.coverage.", sp, ".RData", sep=""))
  
  for (enh in enhancer.datasets[[sp]]){
    ymax=max(c(enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], enh_prop_dist[["obs"]][[paste0(enh,"_confup")]],  enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], enh_prop_dist[["simul"]][[paste0(enh,"_confup")]]))
    ymin=min(c(enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], enh_prop_dist[["obs"]][[paste0(enh,"_confup")]],  enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], enh_prop_dist[["simul"]][[paste0(enh,"_confup")]]))
    
    ymax=ymax*1.1
    
    par(mar=c(3.1, 4.5, 2.75, 2))
    
    plot(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"], main="", type="n", xlab="",ylab="",  axes=F, ylim=c(ymin,ymax))
    
    lines(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"])
    lines(enh_prop_dist[["simul"]][[enh]], col=dataset.colors["Simulated"])
        
    xpos=1:length(enh_prop_dist[["obs"]][[enh]])
    
    segments(xpos, enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], xpos, enh_prop_dist[["obs"]][[paste0(enh,"_confup")]], col=dataset.colors["Original"])
    segments(xpos, enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], xpos, enh_prop_dist[["simul"]][[paste0(enh,"_confup")]], col=dataset.colors["Simulated"])
    
    
    class_leg <- c("0.05", "0.5", "1", "1.5", "2") ## first class has an average of 50kb
    
    axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
    axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)
    
    if (nb == 1){
      legend("topright", legend=c("PCHi-C data", "simulated data"), 
             col=dataset.colors[c("Original", "Simulated")], lty=1, seg.len=1, bty='n', cex=1.1, inset=c(0.05, 0.05), xpd=NA)
    }
    
    
    mtext("% length covered\n by enhancers", side=2, cex=0.85, line=2, at=(ymin+ymax)*0.9/2)
    mtext("distance to promoters (Mb)", side=1, line=2, cex=0.85)
    
    mtext(letters[nb], side=3, line=1, at=-7.75, font=2, cex=1.2)
    mtext(paste(sp, enh, sep=" "), side=3, cex=0.8, line=0)
    
    nb = nb+1
  }
}

dev.off()


############################################################################################################
