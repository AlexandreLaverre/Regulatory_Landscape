###################################################################################

library(data.table)

###################################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T

  enh="ENCODE"
  source("../main_figures/parameters.R")
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "SupplementaryMaterialFigure25.pdf", sep=""), width=6.85, height=8.5)

layout(matrix(1:6, nrow=3, byrow=F))

###################################################################################

labels=c("a", "b", "c", "d", "e", "f")

i=1

###################################################################################

for(ref in c("human", "mouse")){

  tg=setdiff(c("human", "mouse"), ref)
  
  load(paste(pathFigures, "/RData/data.synteny.conservation.", ref, ".RData", sep=""))
  
  synteny=conserv_synteny[[enh]]
  
  obs=synteny[[tg]][["synt_obs"]]
  sim=synteny[[tg]][["synt_simul"]]

  #################################
  
  prop.cons=list()
  conf.low=list()
  conf.high=list()
  
  obs$class_genes=cut(obs$fr_length_genes_inbetween, breaks=seq(from=0, to=1, length=4), include.lowest=T)
  sim$class_genes=cut(sim$fr_length_genes_inbetween, breaks=seq(from=0, to=1, length=4), include.lowest=T)

  classes=levels(obs$class_genes)
  
  for(class in classes){
    prop.cons[[class]]=list()
    conf.low[[class]]=list()
    conf.high[[class]]=list()
    
    this.obs=obs[which(obs$class_genes==class),]
    this.sim=sim[which(sim$class_genes==class),]
    
    this.prop.cons.obs=tapply(this.obs$target_dist, this.obs$class_dist, function(x) length(which(x<=maxDistanceSyntenyTarget))/length(x))
    this.conf.low.cons.obs=tapply(this.obs$target_dist, this.obs$class_dist, function(x) prop.test(length(which(x<=maxDistanceSyntenyTarget)), length(x))$conf.int[1])
    this.conf.high.cons.obs=tapply(this.obs$target_dist, this.obs$class_dist, function(x) prop.test(length(which(x<=maxDistanceSyntenyTarget)), length(x))$conf.int[2])
  
    this.prop.cons.sim=tapply(this.sim$target_dist, this.sim$class_dist, function(x) length(which(x<=maxDistanceSyntenyTarget))/length(x))
    this.conf.low.cons.sim=tapply(this.sim$target_dist, this.sim$class_dist, function(x) prop.test(length(which(x<=maxDistanceSyntenyTarget)), length(x))$conf.int[1])
    this.conf.high.cons.sim=tapply(this.sim$target_dist, this.sim$class_dist, function(x) prop.test(length(which(x<=maxDistanceSyntenyTarget)), length(x))$conf.int[2])
    
    prop.cons[[class]]=list("obs"=this.prop.cons.obs, "sim"=this.prop.cons.sim)
    conf.low[[class]]=list("obs"=this.conf.low.cons.obs, "sim"=this.conf.low.cons.sim)
    conf.high[[class]]=list("obs"=this.conf.high.cons.obs, "sim"=this.conf.high.cons.sim)
  }
 
###################################################################################
  
  nbkept=25
  
  for(class in classes){
    
    par(mar=c(4.5, 3.5, 2.5, 1.1))
    
    this.prop.cons=prop.cons[[class]]
    this.conf.low=conf.low[[class]]
    this.conf.high=conf.high[[class]]
    
    xpos=1:nbkept
    smallx=c(0, 0)
    names(smallx)=c("obs", "sim")
    
    ylim=c(65, 100)
    
    plot(1, type="n", xlim=c(0.5, nbkept+0.5), ylim=ylim, axes=F, xlab="", ylab="")
    
    points(xpos+smallx["obs"], 100*this.prop.cons[["obs"]][1:nbkept], col=dataset.colors["Original"], pch=20)
    segments(xpos+smallx["obs"], 100*this.conf.low[["obs"]][1:nbkept], xpos+smallx["obs"], 100*this.conf.high[["obs"]][1:nbkept], col=dataset.colors["Original"])
    
    points(xpos+smallx["sim"], 100*this.prop.cons[["sim"]][1:nbkept], col=dataset.colors["Simulated"], pch=20)
    segments(xpos+smallx["sim"], 100*this.conf.low[["sim"]][1:nbkept], xpos+smallx["sim"], 100*this.conf.high[["sim"]][1:nbkept], col=dataset.colors["Simulated"])
    
    axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=0.8)
    mtext("% pairs in conserved synteny", side=2, line=2.25, cex=0.75)
    
    
    xax=c(1, 5, 10, 15, 20, 25)
    axlab=as.character(c(0.05, 0.25, 0.5, 0.75, 1, 1.5))
    
    axis(side=1, mgp=c(3, 0.65, 0), at=xax, labels=axlab, cex.axis=0.95)
    mtext("distance to promoters (Mb)", side=1, line=2.1, cex=0.75)
    
    mtext(paste(ref, " vs. ", tg, ", gene density ", class, sep=""), side=3, line=0.5, cex=0.75)
    
    mtext(labels[i], side=3, at=-4, font=2, line=1)

    i=i+1
  }
}

###################################################################################

dev.off()

###################################################################################
