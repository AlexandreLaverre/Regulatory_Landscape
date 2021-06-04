## same as Figure 1 for mouse
objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R") ## paths are defined based on the user name
  
  library(bootBCa, lib=pathRlibs)
  
  set.seed(19)
}

###########################################################################################

## load all necessary data, scripts and libraries for the figure

if(load){
  ref_sp="mouse"
  
  load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.samples.cumulative.interactions.RData", sep=""))

  ## observed and simulated contacts - already bait-other, in the right distance range
  load(paste(pathFigures, "RData/data.bootstrap.nb.cell.types.",ref_sp,".RData",sep=""))

 
  nb_cell_max = 5
  load=FALSE ## we only do this once
}

############################################################################################

## prepare data for figure

if(prepare){
 ## already prepared data for bootrap intervals
    
  ## we finish preparing the data
  prepare=FALSE
}

#############################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "GenomeResearch_Figures/SupplementaryMaterialFigure3.pdf", sep=""), width=6.85, height=5.5)

par(mfrow=c(2,2))
par(mai = c(1, 0.8, 0.5, 0.1)) # internal margins
par(mar = c(3.5, 4.5, 3.1, 1)) # external margins

#################### a - Histogram with number of cell types in which an interaction is observed #####################

b=barplot(as.matrix(pc_nb_celltypes_matrix), beside=T, xlab='',
          names=rep("", dim(pc_nb_celltypes_matrix)[2]), ylim=c(0,80), space=c(0.4,1),
          ylab="", border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
  lwd=1.5,  mgp=c(3, 0.75, 0), axes=F)

axis(side=1, cex.axis=1, at=apply(b, 2, mean), labels=rep("", 6))
axis(side=2, cex.axis=1, las=2)

mtext(colnames(nb_celltypes_matrix), at=apply(b, 2, mean), side=1, line=0.5, cex=0.85)

## axis labels
mtext("number of cell types", side=1, line=2.15, cex=0.8)
mtext("% of interactions", side=2, line=2.5, cex=0.8)

## legend 
legend("topright", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n',
       cex=1.1, inset=c(0.05, -0.2), xpd=NA, title=ref_sp)

##plot label
mtext("A", side=3, line=1.5, at=-3.9, font=2, cex=1.1)

################################################################################################

#################### b - Distribution of number of cell types according to distance ############

ylim=c(0.5, max(c(as.numeric(mean_nb_celltypes_dist["Original",]), as.numeric(mean_nb_celltypes_dist["Simulated",]))))
ylim[2]=ylim[2]+0.5

## plot(as.numeric(mean_dist["Original",]), as.numeric(mean_nb_celltypes_dist["Original",]), type="l", col=dataset.colors["Original"], ylim=ylim, xlab="", ylab="", axes=F)
## lines(as.numeric(mean_dist["Simulated",]), as.numeric(mean_nb_celltypes_dist["Simulated",]), col=dataset.colors["Simulated"], lwd=1.5)

plot(as.numeric(mean_dist["Original",]), as.numeric(mean_nb_celltypes_dist["Original",]), pch=20, col=dataset.colors["Original"], ylim=ylim, xlab="", ylab="", axes=F)
points(as.numeric(mean_dist["Simulated",]), as.numeric(mean_nb_celltypes_dist["Simulated",]), col=dataset.colors["Simulated"], lwd=1.5, pch=20)

## X axis
xax=pretty(range(as.numeric(mean_dist)))
labels=xax/1e6

axis(side=1, at=xax, labels=labels, mgp=c(3, 0.65, 0), cex.axis=1)
axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1, las=2)

## axis labels
mtext("distance between interacting fragments (Mb)", side=1, line=2.15, cex=0.8)
mtext("mean number of cell types", side=2, line=2.5, cex=0.8)

## confidence intervals

for(dataset in rownames(mean_dist)){
  segments(as.numeric(mean_dist[dataset,]), as.numeric(dist_conf_low_celltypes[dataset,]),  as.numeric(mean_dist[dataset,]), as.numeric(dist_conf_high_celltypes[dataset,]), col=dataset.colors[dataset])
}

## legend 

legend("topright", legend=c("PCHi-C data", "simulated data"), col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1.1, inset=c(0.05, -0.2), xpd=NA, title=ref_sp)

## plot label

mtext("B", side=3, line=1.5, at=-4.2e5, font=2, cex=1.1)

################################################################################################

#################### C - Cumulative number of interactions #######

labels=c("C", "D")

names(labels)=c("human", "mouse")

for (sp in c("human", "mouse")){

  nbmax=dim(cumul_int[[sp]][["observed"]])[1]
  xlim=c(0.5, nbmax+0.5)
  ylim=c(0, 100)
  
  plot(100*apply(cumul_int[[sp]][["simulated"]], 1, mean)/max(cumul_int[[sp]][["simulated"]]), pch=19, col=dataset.colors["Simulated"],
       xlab="", ylab="", main="", cex=0.5, mgp=c(2,1,0), axes=F, xlim=xlim, ylim=ylim)

  mtext(sp, side=3, cex=0.95)
  
  points(100*apply(cumul_int[[sp]][["observed"]], 1, mean)/max(cumul_int[[sp]][["observed"]]),  pch=19, col=dataset.colors["Original"], cex=0.5)
  
  axis(side=1, mgp=c(3, 0.65, 0), cex.axis=1)
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1, las=2)
  
  mtext("number of samples", side=1, line=2.25, cex=0.8)
  
  if(sp=="human"){
    legend("bottomright", legend=c("PCHi-C data", "simulated data"), col=dataset.colors[c("Original", "Simulated")], pch=20, bty='n', xpd=NA)
  }
  
  mtext("cumulative % of interactions", side=2, line=2.5, cex=0.8)

  pos.xlab=xlim[1]-diff(xlim)/4
  
  mtext(labels[sp], side=3, line=1, at=pos.xlab, font=2, cex=1.1)
}

###########################################################################################

dev.off()

###########################################################################################
