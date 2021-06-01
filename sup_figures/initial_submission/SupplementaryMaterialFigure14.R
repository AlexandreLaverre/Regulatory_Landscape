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
  ref_sp = "human"
  
  load(paste(pathFigures, "RData/data.promoter.enhancer.correlation.", ref_sp, ".RData", sep=""))
  
  enhancers = setdiff(enhancer.datasets[[ref_sp]], "ENCODE") ## ENCODE shown in main figure
}

#################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#################################################################################################################

pdf(paste(pathFigures, "SupplementaryMaterialFigure14.pdf", sep=""), width=3.4, height=6.85)

############################################################################################################################
############################################## correlation gene expression and enhancers activity ##########################

par(mfrow=c(3,1))

labels=c("a", "b", "c")
names(labels)=enhancers

for (enh in enhancers){
  ymin=min(c(correl_activity[["obs"]][[paste0(enh,"_conflow")]], correl_activity[["obs"]][[paste0(enh,"_confup")]],  correl_activity[["simul"]][[paste0(enh,"_conflow")]], correl_activity[["simul"]][[paste0(enh,"_confup")]]))
  ymax=max(c(correl_activity[["obs"]][[paste0(enh,"_conflow")]], correl_activity[["obs"]][[paste0(enh,"_confup")]],  correl_activity[["simul"]][[paste0(enh,"_conflow")]], correl_activity[["simul"]][[paste0(enh,"_confup")]]))
  
  ylim=c(ymin, ymax+0.01)
  
  par(mar=c(3.1, 4.5, 3, 1))
  
  plot(correl_activity[["obs"]][[enh]], type="n", ylab="", main="", las=2, ylim=ylim, axes=F)
  
  ## lines(correl_activity[["obs"]][[enh]], col=dataset.colors["Original"])
  ## lines(correl_activity[["simul"]][[enh]], col=dataset.colors["Simulated"])

  points(correl_activity[["obs"]][[enh]], col=dataset.colors["Original"], pch=20)
  points(correl_activity[["simul"]][[enh]], col=dataset.colors["Simulated"], pch=20)

  
  xpos=1:length(correl_activity[["obs"]][[enh]])
  
  segments(xpos, correl_activity[["obs"]][[paste0(enh,"_conflow")]], xpos, correl_activity[["obs"]][[paste0(enh,"_confup")]], col=dataset.colors["Original"])
  segments(xpos, correl_activity[["simul"]][[paste0(enh,"_conflow")]], xpos, correl_activity[["simul"]][[paste0(enh,"_confup")]], col=dataset.colors["Simulated"])
  
  
  
  class_leg <- c("0.05", "0.5", "1", "1.5", "2") ## first class is from 25kb to 75kb
  axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
  axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)
  
  mtext("Spearman's rho", side=2, line=3, cex=0.8)
  mtext("distance to promoters (Mb)", side=1, line=2, cex=0.8)
  
  mtext(enh.syn[[enh]], side=3, line=0, cex=0.8)
  
  if (enh == "FANTOM5"){
    legend("topright", legend=c("PCHi-C data", "simulated data"), col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1.1, inset=c(0.05, 0.1), xpd=NA)
    
  }

  mtext(labels[enh], side=3, line=1, at=-7.7, font=2, cex=1.2)
}

##################################################################################################################################

dev.off()


##################################################################################################################################
