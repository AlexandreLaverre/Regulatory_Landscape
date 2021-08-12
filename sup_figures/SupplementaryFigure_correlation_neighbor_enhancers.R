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

pdf(paste(pathFigures, "GenomeResearch_Figures/SupplementaryFigure_correlation_neighbor_enhancers.pdf", sep=""), width=6.85, height=2.9)

m=matrix(c(1,2, 3), nrow=1)
layout(m)

############################################################################################################################
############################################## correlation gene expression and enhancers activity ##########################

labels=c("A", "B", "C")
index=1

for(sp in c("human", "mouse")){
  for(enh in c("ENCODE", "FANTOM5")){
    
    if(length(obs[[sp]][[enh]])>0){
      
      ymin=min(c(obs[[sp]][[paste0(enh,"_conflow")]], obs[[sp]][[paste0(enh,"_confup")]], neighbors[[sp]][[paste0(enh,"_conflow")]], neighbors[[sp]][[paste0(enh,"_confup")]]))
      ymax=max(c(obs[[sp]][[paste0(enh,"_conflow")]], obs[[sp]][[paste0(enh,"_confup")]], neighbors[[sp]][[paste0(enh,"_conflow")]], neighbors[[sp]][[paste0(enh,"_confup")]]))
      
      ylim=c(ymin, ymax)
      xlim=c(0.5, length(obs[[sp]][[enh]])+0.5)
      
      par(mar=c(3.5, 4.5, 2.1, 1))
      
      plot(obs[[sp]][[enh]],  main="", las=2, ylim=ylim, axes=F, pch=20, col=dataset.colors[["Original"]], xlab="", ylab="")
      points(neighbors[[sp]][[enh]], col="black", pch=20)
      
      xpos=1:length(obs[[sp]][[enh]])
      
      segments(xpos, obs[[sp]][[paste0(enh,"_conflow")]], xpos, obs[[sp]][[paste0(enh,"_confup")]], col=dataset.colors["Original"])
      segments(xpos, neighbors[[sp]][[paste0(enh,"_conflow")]], xpos, neighbors[[sp]][[paste0(enh,"_confup")]], col="black")
            
      class_leg <- c("0", "0.25", "0.5", "0.75", "1")
      axis(side=1, at=c(0.5, 8.5, 19.5, 29.5, 39.5), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=0.95)
      axis(side=2, mgp=c(3, 0.65, 0), cex.axis=0.95, las=2)
      
      mtext("Spearman's rho", side=2, cex=0.8, line=2.5, at=(ymin+ymax)*0.9/2)
      mtext("distance to promoters (Mb)", side=1, line=2, cex=0.8)
      mtext(paste(sp, enh), side=3, line=0.5, cex=0.8)

      mtext(labels[index], side=3, line=0.5, font=2, at=-10)
      
      index=index+1
    }
  }
}

##################################################################################################################################

dev.off()

##################################################################################################################################
