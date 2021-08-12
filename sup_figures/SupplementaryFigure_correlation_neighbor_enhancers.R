#################################################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
}

source("../main_figures/parameters.R") ## paths are defined based on the user name

#################################################################################################################

if(load){
 
  obs=list()
  neighbors=list()
  
  for(ref_sp in c("human", "mouse")){
    load(paste(pathFigures, "RData/data.promoter.enhancer.correlation.neighbor.enhancers.",ref_sp,".RData", sep=""))
    
    obs[[ref_sp]]=correl_activity[["obs"]]
    neighbors[[ref_sp]]=correl_activity[["neighbors"]]
  }
  
}

#################################################################################################################

pdf(paste(pathFigures, "GenomeResearch_Figures/SupplementaryFigure_correlation_neighbor_enhancers.pdf", sep=""), width=6.85, height=3.2)

m=matrix(c(1,2, 3), nrow=1)
layout(m)

############################################################################################################################
############################################## correlation gene expression and enhancers activity ##########################

for(sp in c("human", "mouse")){
  for(enh in c("ENCODE", "FANTOM5")){
    
    if(length(obs[[sp]][[enh]])>0){
      
      ymin=min(c(obs[[sp]][[paste0(enh,"_conflow")]], obs[[sp]][[paste0(enh,"_confup")]], neighbors[[sp]][[paste0(enh,"_conflow")]], neighbors[[sp]][[paste0(enh,"_confup")]]))
      ymax=max(c(obs[[sp]][[paste0(enh,"_conflow")]], obs[[sp]][[paste0(enh,"_confup")]], neighbors[[sp]][[paste0(enh,"_conflow")]], neighbors[[sp]][[paste0(enh,"_confup")]]))
      
      ylim=c(ymin, ymax)
      
      par(mar=c(3.1, 4.5, 1, 1))
      
      plot(obs[[sp]][[enh]], ylab="", main="", las=2, ylim=ylim, axes=F, pch=20, col=dataset.colors[["Original"]])
      points(neighbors[[sp]][[enh]], col="black", pch=20)
      
      xpos=1:length(obs[[sp]][[enh]])
      
      segments(xpos, obs[[sp]][[paste0(enh,"_conflow")]], xpos, obs[[sp]][[paste0(enh,"_confup")]], col=dataset.colors["Original"])
      segments(xpos, neighbors[[sp]][[paste0(enh,"_conflow")]], xpos, neighbors[[sp]][[paste0(enh,"_confup")]], col="black")
      
      
      class_leg <- c("0", "0.5", "1", "1.5", "2")
      axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=0.9)
      axis(side=2, mgp=c(3, 0.65, 0), cex.axis=0.9, las=2)
      
      mtext("Spearman's rho", side=2, cex=0.85, line=3, at=(ymin+ymax)*0.9/2)
      mtext("distance to promoters (Mb)", side=1, line=2, cex=0.85)
    }
  }
}

##################################################################################################################################

dev.off()


##################################################################################################################################
