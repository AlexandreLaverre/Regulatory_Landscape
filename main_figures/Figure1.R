###########################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("parameters.R") ## paths are defined based on the user name

  pathEnhancers="../../../RegulatoryLandscapesManuscript/SupplementaryDataset4/"

  library(plotrix)

  set.seed(19)
}

###########################################################################################

## load all necessary data, scripts and libraries for the figure

if(load){
 
  library(ape)
    
  sp="human"
  
  ## functions for genome browser plots

  source(paste(pathScripts, "/genome_browser/plot.annotations.R", sep=""))
  source(paste(pathScripts, "/genome_browser/plot.interactions.R", sep=""))
  
  load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.Shh.figure.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.clustering.RData", sep=""))

  ## observed and simulated contacts - already bait-other, in the right distance range
  load(paste(pathFigures, "RData/data.bootstrap.nb.cell.types.",sp,".RData",sep=""))

  ## data for simulations
  load(paste(pathFigures, "RData/data.bait.annotation.RData",sep=""))

  baits=bait.info[[sp]]

  load(paste(pathFigures, "RData/data.restriction.map.RData", sep=""))

  map=restriction.map[[sp]]

  ## enhancer coordinates, ENCODE
  enhancers=list()
  
  for(enh in enhancer.datasets[[sp]]){
    enhancers[[sp]]=read.table(paste(pathEnhancers, sp, "/",enh,"/enhancer_coordinates.bed", sep=""), h=T, stringsAsFactors=F)
  }

  load=FALSE ## we only do this once
}

############################################################################################

## prepare data for figure

if(prepare){
  
  ## dendrogram based on the % of shared interactions between samples, observed-simulated 

  hcl=sample.clustering[[sp]][["hclust.alldist"]]
  sample.order=sample.clustering[[sp]][["sample.order.alldist"]]

  ## statistics for the test

  print(paste("interactions shared across cell types"))
  print(paste("observed", round(100*length(which(obs$nb_celltypes>=2))/nrow(obs))))
  print(paste("simulated", round(100*length(which(sim$nb_celltypes>=2))/nrow(sim))))

  print(paste("long range (>500kb) shared across cell types"))
  print(paste("observed", round(100*length(which(obs$nb_celltypes>=2 & obs$distance>=5e5))/length(which(obs$distance>=5e5)))))
  print(paste("simulated", round(100*length(which(sim$nb_celltypes>=2 & sim$distance>=5e5))/length(which(sim$distance>=5e5)))))

  ## data for simulations
  bait="chr11:66263828:66266552"
  cell="CD34"
  
  this.obs=obs[which(obs$id_bait==bait & !is.na(obs[,cell])),]
  this.sim=sim[which(sim$id_bait==bait & !is.na(sim[,cell])),]

  chr=this.obs$chr_bait[1]

  xrange=range(c(this.obs$start, this.sim$start, this.obs$end, this.sim$end))

  this.map=map[which(map$chr==chr & map$start>=xrange[1] & map$end<=xrange[2]),]

  this.baits=baits[which(baits$chr==chr & baits$start>=xrange[1] & baits$end<=xrange[2]),]

  this.enhancers=list()
  
  for(enh in enhancer.datasets[[sp]]){
    all.enhancers=enhancers[[sp]]
    this.enhancers[[enh]]=all.enhancers[which(all.enhancers$chr==chr & all.enhancers$start>=xrange[1] & all.enhancers$end<=xrange[2]),]
}

  ## all obs, only this sample

  all.obs=obs[which(!is.na(obs[,cell])),]
  all.sim=sim[which(!is.na(sim[,cell])),]
  
  ## we finish preparing the data
  prepare=FALSE
}

#############################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "GenomeResearch_Figures/Figure1.pdf", sep=""), width=6.85, height=11)

## layout

m=matrix(rep(NA, 50*10), nrow=50)

for(i in 1){
  m[i,]=c(rep(1, 6), rep(6, 4))
}

for(i in 2:5){
  m[i,]=c(rep(2, 6), rep(6, 4))
}

for(i in 6:7){
 m[i,]=c(rep(3, 6), rep(6, 4))
}

for(i in 8:12){
 m[i,]=c(rep(4, 6), rep(7, 4))
}

for(i in 13){
 m[i,]=c(rep(5, 6), rep(7, 4))
}

## SHH example

for(i in 14){
  m[i,]=c(rep(17,10)) ## empty plot
}

for(i in 15:17){
  m[i,]=c(rep(8,1), rep(9,9)) ## annotations
}

for(i in c(18)){
  m[i,]=c(rep(10, 1), rep(11,9))
}

for(i in c(19:40)){
  m[i,]=c(rep(12, 1), rep(13,9))
}

## ## enhancers
## for(i in 41:44){
##   m[i,]=c(rep(14, 1), rep(15,9))
## }
for(i in 41){
  m[i,]=c(rep(14, 10))
}

for(i in 42:50){
  m[i,]=c(rep(15, 5), rep(16, 5))
}


layout(m)

############################################################################################
## margin

par(mar=c(0.2, 5.2, 0.2, 1.2))

##########################################################################

## 1. legend plot

plot(1, type="n", xlim=c(0,1.15), ylim=c(0,1), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")

legend("topright", c("PCHi-C data", "simulated data"), col=dataset.colors, lty=1, inset=0.01, bty="n", xpd=NA, cex=0.95, seg.len=1)

smally=-0.72

rect(0.792, -0.85+smally, 0.80, -0.45+smally, col="gray80", border="gray80", xpd=NA)
rect(0.81, -0.85+smally, 0.825, -0.45+smally, col="gray20", border="gray20", xpd=NA)
rect(0.835, -0.85+smally, 0.845, -0.45+smally, col="gray80", border="gray80", xpd=NA)

text(x=0.865, y=-0.6+smally, labels="restriction fragments", cex=0.95, xpd=NA, adj=c(0,0.5))

rect(0.80, -1.55+smally, 0.802, -1.15+smally, col="gray20", border="gray20", xpd=NA)
rect(0.84, -1.55+smally, 0.842, -1.15+smally, col="gray20", border="gray20", xpd=NA)
rect(0.82, -1.55+smally, 0.822, -1.15+smally, col="gray20", border="gray20", xpd=NA)

text(x=0.865, y=-1.325+smally, labels="enhancers", cex=0.95, xpd=NA, adj=c(0,0.5))

rect(0.82, -2.2+smally, 0.835, -1.85+smally, col="gray80", border="red", xpd=NA)

text(x=0.865, y=-2+smally, labels="bait", cex=0.95, xpd=NA, adj=c(0,0.5))

mtext("example bait, one sample", side=3, line=-1, at=0.15, cex=0.65, font=1)

mtext("A", side=3, line=-1, at=-0.19, cex=1, font=2, xpd=NA)

##########################################################################

## 2. observed interactions

yrange=c(0, 1)

plot(1, type="n", xlim=xrange, ylim=c(0,0.6), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")

for(i in 1:nrow(this.obs)){
  pos.bait=(this.obs$start_bait[i]+this.obs$end_bait[i])/2
  pos.frag=(this.obs$start[i]+this.obs$end[i])/2
  x.center=(pos.bait+pos.frag)/2
  R.x=abs(x.center-pos.bait)
  R.y=0.5
  
  draw.ellipse(x=x.center, y=0, a=R.x, b=0.5, col=NA, border=dataset.colors["Original"])
}

##########################################################################

## 3. restriction fragments

plot(1, type="n", xlim=xrange, ylim=c(-1,1), axes=F, xlab="", ylab="", xaxs="i")

for(i in 1:dim(this.map)[1]){
  if(i%%2==0){
    rect(this.map$start[i], 0.25, this.map$end[i], 0.75, col="gray80", border=NA)
  } else{
    rect(this.map$start[i], 0.25, this.map$end[i], 0.75, col="gray20", border=NA)
  }
}

mtext("fragments", side=2, las=2, at=0.55, cex=0.65,line=1)

## coordinates for the bait

rect(this.obs$start_bait[1], 0.25, this.obs$end_bait[1], 0.75, col="gray80", border="red")
abline(h=-0.5, lty=1, col="gray80")

this.enh=this.enhancers[["ENCODE"]]

for(i in 1:dim(this.enh)[1]){
  rect(this.enh$start[i], -0.25, this.enh$end[i], -0.75, col="gray20", border=NA)
}

mtext("enhancers", side=2, las=2, at=-0.45, cex=0.65,line=1)

##########################################################################

## 4. simulated interactions

yrange=c(0, 1)

plot(1, type="n", xlim=xrange, ylim=c(-0.6,0), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")

for(i in 1:nrow(this.sim)){
  pos.bait=(this.sim$start_bait[i]+this.sim$end_bait[i])/2
  pos.frag=(this.sim$start[i]+this.sim$end[i])/2
  x.center=(pos.bait+pos.frag)/2
  R.x=abs(x.center-pos.bait)
  R.y=0.5
  
  draw.ellipse(x=x.center, y=0, a=R.x, b=-0.5, col=NA, border=dataset.colors["Simulated"])
  
}

##########################################################################
## 5. axis plot

plot(1, type="n", xlim=xrange, ylim=c(-0.6,0), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")

xaxs=pretty(xrange/1e6, n=10)
axis(side=1, mgp=c(3, 0.5, 0), cex.axis=0.9, line=-1.5, at=xaxs*1e6, labels=paste0(xaxs, "Mb"))

mtext(chr, side=1, line=-2, at=xrange[1]-diff(xrange)/20, cex=0.55)

##########################################################################

## 6. density plot, distance, observed interactions

d.obs=density(all.obs$distance, bw=0.5)
d.sim=density(all.sim$distance, bw=0.5)

xlim=range(c(0, maxDistance))
ylim=range(c(d.obs$y, d.sim$y))

xaxs=c(0, 500e3, 1e6, 1.5e6, 2e6, 2.5e6)
xaxslabels=c("0", "0.5", "1", "1.5", "2", "2.5")  

par(mar=c(2.8, 7.1, 1.1, 0.25))
plot(d.obs$x, d.obs$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Original"])
axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.9)
mtext("distance bait-fragment (Mb)", side=1, line=1.5, cex=0.65)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.9)

mtext("density, nb. contacts", side=2, line=2.1, cex=0.65)

mtext("all baits, one sample", side=3, line=0, at=xlim[2], adj=1, cex=0.65, font=1)

mtext("B", side=3, line=0, at=-0.75e6, cex=1, font=2, xpd=NA)

##########################################################################

## 7. density plot, distance, simulated interactions

plot(d.sim$x, d.sim$y, type="l", xlab="", ylab="", xlim=xlim, ylim=ylim, axes=F, col=dataset.colors["Simulated"])
axis(side=1, mgp=c(3, 0.5, 0), at=xaxs, labels=xaxslabels, cex.axis=0.9)
mtext("distance bait-fragment (Mb)", side=1, line=1.5, cex=0.65)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.9)

mtext("density, nb. contacts", side=2, line=2.1, cex=0.65)

mtext("C", side=3, line=0, at=-0.75e6, cex=1, font=2, xpd=NA)

##########################################################################

## 8. label plot
print("label plot")
## empty plot
par(mar=c(0, 0.1, 0, 0.3)) 
plot.new()

############################################################################################

## 9. annotations
print("annotations")
## annotations in the Shh region

par(mar=c(0.25, 0.5, 2.0, 9.75))
## plot(1, type="n", xlab="", ylab="", axes=F, xlim=shhxlim, ylim=c(0,1), xaxs="i", yaxs="i")

## quick fix for gene names, will be removed later
shhgenecoords$name=rep(NA, dim(shhgenecoords)[1])
shhgenecoords$name[which(shhgenecoords$id==shhid)]="SHH"
shhgenecoords$name[which(shhgenecoords$id=="ENSG00000105983")]="LMBR1"
shhgenecoords$name[which(shhgenecoords$id=="ENSG00000182648")]="LINC01006" 

plot.annotations.genes(gene.coords=shhgenecoords, focus.gene=shhid, gene.biotypes=c("protein_coding", "lincRNA"), xlim=shhxlim, col.focus=col.Shh, col.other="gray60", axis=F, axisunit="Mb", axisside=3, cex.name=0.8, name.position="top", show.arrows=T, highlighted.genes=c("ENSG00000105983", "ENSG00000182648"))

## axis label

mtext("genes", side=2, las=2, cex=0.7, line=1.75)


## plot label

mtext("D", side=3, line=0.5, at=shhxlim[1]-diff(shhxlim)/7.5, font=2, cex=1.1)

#############################################################################################

## empty plot
par(mar=c(0, 0.1, 0, 0.3)) 
plot.new()

#############################################################################################

## baits in the Shh region

par(mar=c(0.65, 0.5, 0, 10.5)) ## left and right margin should be the same as above
plot(1, type="n", xlab="", ylab="", axes=F, xlim=shhxlim, ylim=c(0,1), xaxs="i", yaxs="i")

segments(shhxlim[1], 0.75, shhxlim[2], 0.75, lwd=0.5, lty=3, col="gray40")

rect(allshhbaits$start, 0.5, allshhbaits$end, 1, col="gray40", border="gray40")

xax=pretty(shhxlim/1e6)
axis(side=1, at=xax*1e6, labels=paste(xax, "Mb", sep=""), cex.axis=0.85, line=-0.25, mgp=c(3, 0.35,0))

## shhchr
mtext("chr7", at=shhxlim[2]+diff(shhxlim)/20, line=-1.5, side=3, cex=0.7)

mtext("baits", side=2, las=2, cex=0.7, line=1.75, at=0.75)

##############################################################################################

## Dendrogram of samples
par(mar=c(0.5, 0.15, 0.105, 0.3)) 
plot(as.phylo(hcl), direction="rightwards", show.tip.label=FALSE)

##############################################################################################

## interactions in the Shh region
par(mar=c(0.5, 0.5, 0.1, 10.5)) ## left and right margin should be the same as above

ylim=c(0, length(samples)+1)
height=0.25
ypos=1:length(samples)
names(ypos)=sample.order # re-order according to AFC

plot(1, type="n", xlab="", ylab="", axes=F, xlim=shhxlim, ylim=ylim, xaxs="i", yaxs="i")

for(sample in samples){
  abline(h=ypos[sample], col="gray40",lwd=0.5, lty=3)
  
  this.int=shhinteractions[which(!is.na(shhinteractions[,sample])),]
  if(dim(this.int)[1]>0){
    rect(this.int$start, ypos[sample]-height, this.int$end, ypos[sample]+height, col=col.Shh, border=NA)
  }
}

## labels for the cell types

ywidth=diff(ypos)[1]
segx=shhxlim[2]+diff(shhxlim)/100

for(c in unique(celltypes)){
  all.ypos=ypos[samples[which(celltypes[samples]==c)]]
  segments(segx, min(all.ypos)-ywidth/3, segx, max(all.ypos)+ywidth/3, xpd=NA)
  
  mtext(syn.celltypes[c], side=4, line=0.75, las=2, cex=0.65, at=mean(all.ypos))
}

##############################################################################################

par(mar=c(0, 0.1, 0, 0.3)) 
plot.new()

###############################################################################################

par(mar = c(3.25, 3.75, 0.25, 1)) # external margins

#################### Histogram with number of samples in which an interaction is observed #####################

b=barplot(as.matrix(pc_nb_celltypes_matrix), beside=T, xlab='',
          names=rep("", dim(pc_nb_celltypes_matrix)[2]), ylim=c(0,80), space=c(0.4,1),
          ylab="", border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
          lwd=1.5,  mgp=c(3, 0.75, 0), cex.axis=1.1, las=2)

mtext(colnames(nb_celltypes_matrix), at=apply(b, 2, mean), side=1, line=0.5, cex=0.75)

## axis labels
mtext("number of cell types", side=1, line=2, cex=0.75)
mtext("% interactions", side=2, line=2.5, cex=0.75)

## legend & plot label
##legend("topright", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],  fill=dataset.colors[c("Original", "Simulated")], bty='n', cex=1.1, inset=c(0.05, -0.1), xpd=NA)
mtext("E", side=3, line=1, at=-3.9, font=2, cex=1.1)

################################################################################################

#################### Distribution of number of cell types according to distance #####################

ylim=c(0, max(c(as.numeric(mean_nb_celltypes_dist["Original",]), as.numeric(mean_nb_celltypes_dist["Simulated",]))))
ylim[2]=ylim[2]+1

plot(as.numeric(mean_dist["Original",]), as.numeric(mean_nb_celltypes_dist["Original",]), pch=20, col=dataset.colors["Original"], ylim=ylim, xlab="", ylab="", axes=F)
points(as.numeric(mean_dist["Simulated",]), as.numeric(mean_nb_celltypes_dist["Simulated",]), pch=20, col=dataset.colors["Simulated"], lwd=1.5)


## X axis
xax=pretty(range(as.numeric(mean_dist)))
labels=xax/1e6
axis(side=1, at=xax, labels=labels, mgp=c(3, 0.65, 0), cex.axis=1.1)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1, las=2)

## axis labels
mtext("distance between interacting fragments (Mb)", side=1, line=2, cex=0.75)
mtext("mean number of cell types", side=2, line=2.5, cex=0.75)

## confidence intervals

for(dataset in rownames(mean_dist)){
  segments(as.numeric(mean_dist[dataset,]), as.numeric(dist_conf_low_celltypes[dataset,]),  as.numeric(mean_dist[dataset,]), as.numeric(dist_conf_high_celltypes[dataset,]), col=dataset.colors[dataset])
}

## legend & plot label

## legend("topright", legend=c("PCHi-C data", "simulated data"), col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1.1, inset=c(0.05, -0.1), xpd=NA)

mtext("F", side=3, line=1, at=-3.15e5, font=2, cex=1.1)

###########################################################################################

############################################################################################
dev.off()

###########################################################################################
