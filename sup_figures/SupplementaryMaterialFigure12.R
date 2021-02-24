#################################################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  prepare=T
}

source("../main_figures/parameters.R") ## paths are defined based on the user name

###########################################################################################################################
##########################  Enhancer proportion according to distance ############################################

pdf(paste(pathFigures, "SupplementaryMaterialFigure12.pdf", sep=""), width=6.85, height=7.5)

par(mai = c(0.5, 0.5, 0.1, 0.2)) # bottom, left, top, right

par(mfrow=c(3,2))

nb=1

for (sp in c("human", "mouse")){

  enhancers = enhancer.datasets[[sp]]
  
  all.enh.prop=data.frame(id=character(0), data=numeric(0), conf_up=numeric(0), conf_low=numeric(0))
  all.enh.prop.nbcell=list("obs"=list(), "simul"=list())
  all.enh.prop.dist=list("obs"=list(), "simul"=list())
  
  for(enh in enhancers){
    load(paste(pathFigures, "RData/data.enhancer.coverage.", sp,".",enh,".RData", sep=""))
    
    all.enh.prop=rbind(all.enh.prop, enh_prop, stringsAsFactors=F)
    
    all.enh.prop.nbcell[["obs"]]=c(all.enh.prop.nbcell[["obs"]], enh_prop_nb_cell[["obs"]])
    all.enh.prop.nbcell[["simul"]]=c(all.enh.prop.nbcell[["simul"]], enh_prop_nb_cell[["simul"]])
    
    all.enh.prop.dist[["obs"]]=c(all.enh.prop.dist[["obs"]], enh_prop_dist[["obs"]])
    all.enh.prop.dist[["simul"]]=c(all.enh.prop.dist[["simul"]], enh_prop_dist[["simul"]])
  }
  
  enh_prop=all.enh.prop
  enh_prop_dist=all.enh.prop.dist
  enh_prop_nb_cell=all.enh.prop.nbcell
  
  
  for (enh in enhancers){
    ymax=max(c(enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], enh_prop_dist[["obs"]][[paste0(enh,"_confup")]],  enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], enh_prop_dist[["simul"]][[paste0(enh,"_confup")]]))
    ymin=min(c(enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], enh_prop_dist[["obs"]][[paste0(enh,"_confup")]],  enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], enh_prop_dist[["simul"]][[paste0(enh,"_confup")]]))
    
    ymax=ymax*1.1
    
    par(mar=c(3.1, 4.5, 2.75, 2))
    
    plot(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"], main="", type="n", xlab="",ylab="",  axes=F, ylim=c(ymin,ymax))
    
    ## lines(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"])
    ## lines(enh_prop_dist[["simul"]][[enh]], col=dataset.colors["Simulated"])

    points(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"], pch=20)
    points(enh_prop_dist[["simul"]][[enh]], col=dataset.colors["Simulated"], pch=20)

          
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
    mtext("distance to baits (Mb)", side=1, line=2, cex=0.85)
    
    mtext(letters[nb], side=3, line=1, at=-7.75, font=2, cex=1.2)
    mtext(paste(sp, enh.syn[enh], sep=" "), side=3, cex=0.8, line=0)
    
    nb = nb+1
  }
}

dev.off()


############################################################################################################
