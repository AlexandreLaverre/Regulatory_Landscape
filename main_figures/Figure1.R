###########################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("parameters.R") ## paths are defined based on the user name
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

  load=FALSE ## we only do this once
}

############################################################################################

## prepare data for figure

if(prepare){
  ## observed and simulated contacts - already bait-other, in the right distance range
  obs=observed.contacts[[sp]]
  sim=simulated.contacts[[sp]]

  info=sampleinfo[[sp]]
  rownames(info)=info$Sample.ID
  
  samples=info$Sample.ID 
  celltypes=info$Broad.cell.type.or.tissue
  names(celltypes)=samples
  
  print(paste("there are", length(samples), "samples"))
  
  breaks_samples=c(0, 1, 5, 10, 15, 20, length(samples))
  nb_samples_names=c("1", "2-5", "6-10", "11-15", "16-20", paste("21-", length(samples), sep=""))
    
  obs$nb_samples <- apply(obs[,samples], 1, function(x) sum(!is.na(x)))
  obs$sample_class <- cut(obs$nb_samples, breaks=breaks_samples, include.lowest = T)
  obs$dist_class <- cut(obs$distance, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)

  obs$nb_celltypes <- apply(obs[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
  obs$celltype_class<- cut(obs$nb_celltypes, breaks=c(0:7, max(obs$nb_celltypes)), include.lowest=T)
  levels(obs$celltype_class)=c(as.character(1:7), ">7")
  
  sim$nb_samples <- apply(sim[,samples], 1, function(x) sum(!is.na(x)))
  sim$sample_class <- cut(sim$nb_samples, breaks=breaks_samples, include.lowest = T)
  sim$dist_class <- cut(sim$distance, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)

  sim$nb_celltypes <- apply(sim[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
  sim$celltype_class<- cut(sim$nb_celltypes, breaks=c(0:7, max(sim$nb_celltypes)), include.lowest=T)
  levels(sim$celltype_class)=c(as.character(1:7), ">7")

  filtered_data <- list("Original"=obs, "Simulated"=sim) ## data is already unbaited, in cis, in the right distance range
 
  ## compute number of interactions in each nb samples class
  
  nb_samples_matrix <- sapply(filtered_data, function(x) as.numeric(table(x$sample_class)))
  rownames(nb_samples_matrix) = nb_samples_names

  nb_samples_matrix=t(nb_samples_matrix)
  pc_nb_samples_matrix=100*nb_samples_matrix/apply(nb_samples_matrix,1, sum)

  ## same for cell types
  
  nb_celltypes_matrix <- sapply(filtered_data, function(x) as.numeric(table(x$celltype_class)))
  rownames(nb_celltypes_matrix) = levels(obs$celltype_class)

  nb_celltypes_matrix=t(nb_celltypes_matrix)
  pc_nb_celltypes_matrix=100*nb_celltypes_matrix/apply(nb_celltypes_matrix,1, sum)

  ## divide interactions based on distance, compute mean number of samples by distance class

  mean_nb_celltypes_dist <- t(sapply(filtered_data, function(x)   tapply(x$nb_celltypes, as.factor(x$dist_class), mean, na.rm=T)))
  mean_dist <- t(sapply(filtered_data, function(x)   tapply(x$distance, as.factor(x$dist_class), mean, na.rm=T)))
  
  dist_conf_low_celltypes <- t(sapply(filtered_data, function(x) tapply(x$nb_celltypes, as.factor(x$dist_class), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
  dist_conf_high_celltypes <- t(sapply(filtered_data, function(x) tapply(x$nb_celltypes, as.factor(x$dist_class), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

  # dendrogram based on the % of shared interactions between samples, observed-simulated 

  hcl=sample.clustering[[sp]][["hclust.alldist"]]
  sample.order=sample.clustering[[sp]][["sample.order.alldist"]]
  
  ## we finish preparing the data
  prepare=FALSE
}

#############################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "Figure1.pdf", sep=""), width=6.85, height=7.5)

## layout

m=matrix(rep(NA, 51*10), nrow=51)

for(i in 1:6){
  m[i,]=c(rep(1,1), rep(2,9))
}

for(i in c(7)){
  m[i,]=c(rep(3, 1), rep(4,9))
}

for(i in c(8:29)){
  m[i,]=c(rep(5, 1), rep(6,9))
}


for(i in 30:32){
  m[i,]=c(rep(7, 1), rep(8,9))
}

for(i in 33:51){
  m[i,]=c(rep(9, 5), rep(10, 5))
}

layout(m)

############################################################################################

## empty plot
par(mar=c(0, 0.1, 0, 0.3)) 
plot.new()

############################################################################################

## annotations in the Shh region

par(mar=c(0.25, 0.5, 2.0, 9.75))
## plot(1, type="n", xlab="", ylab="", axes=F, xlim=shhxlim, ylim=c(0,1), xaxs="i", yaxs="i")

## quick fix for gene names, will be removed later
shhgenecoords$name=rep(NA, dim(shhgenecoords)[1])
shhgenecoords$name[which(shhgenecoords$id==shhid)]="SHH"
shhgenecoords$name[which(shhgenecoords$id=="ENSG00000105983")]="LMBR1"
shhgenecoords$name[which(shhgenecoords$id=="ENSG00000182648")]="LINC01006" 

plot.annotations.genes(gene.coords=shhgenecoords, focus.gene=shhid, gene.biotypes=c("protein_coding", "lincRNA"), xlim=shhxlim, col.focus=col.Shh, col.other="gray60", axis=T, axisunit="Mb", axisside=3, cex.name=0.9, name.position="top", show.arrows=T, highlighted.genes=c("ENSG00000105983", "ENSG00000182648"))

## axis label

mtext("genes", side=2, las=2, cex=0.7, line=1.75)

## shhchr

mtext("chr7", at=shhxlim[2]+diff(shhxlim)/20, line=0.5, side=3, cex=0.75)

## plot label

mtext("a", side=3, line=0.5, at=shhxlim[1]-diff(shhxlim)/7.5, font=2, cex=1.2)

#############################################################################################

## empty plot
par(mar=c(0, 0.1, 0, 0.3)) 
plot.new()

#############################################################################################

## baits in the Shh region

par(mar=c(0.5, 0.5, 0.1, 10.5)) ## left and right margin should be the same as above
plot(1, type="n", xlab="", ylab="", axes=F, xlim=shhxlim, ylim=c(0,1), xaxs="i", yaxs="i")

segments(shhxlim[1], 0.5, shhxlim[2], 0.5, lwd=0.5, lty=3, col="gray40")

rect(allshhbaits$start, 0.15, allshhbaits$end, 0.85, col="gray40", border="gray40")

mtext("baits", side=2, las=2, cex=0.7, line=1.75)

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

##############################################################################################

## enhancers in the Shh region

par(mar=c(0.15, 0.5, 0.1, 10.5)) ## left and right margin should be the same as above
plot(1, type="n", xlab="", ylab="", axes=F, xlim=shhxlim, ylim=c(0,1), xaxs="i", yaxs="i")

ypos=seq(from=0.1, to=0.9, length=length(shhenhancers))
names(ypos)=names(shhenhancers)

height=0.05

for(ed in names(shhenhancers)){
  this.enhancers=shhenhancers[[ed]]
  segments(shhxlim[1], ypos[ed], shhxlim[2], ypos[ed], lwd=0.5, lty=3, col="gray40")
  rect(this.enhancers$start, ypos[ed]-height, this.enhancers$end, ypos[ed]+height, col="gray40", border=NA)

  mtext(enh.syn[ed], side=4, at=ypos[ed], line=0.35, las=2, cex=0.65)
}

mtext("predicted", side=2, las=2, cex=0.7, line=0.75, at=0.75)
mtext("enhancers", side=2, las=2, cex=0.7, line=0.75, at=0.4)

## Zrs enhancer

## ZRS coordinates from Vista database: chr7:156,791,088-156,791,875 (after liftOver transformation to hg38)

zrspos=(156791088+156791875)/2

segments(zrspos, 0,  zrspos, 1,col="gray40")

mtext("ZRS", font=1, cex=0.7, side=3, at=zrspos, line=0.15)

###############################################################################################

par(mai = c(1, 0.8, 0.5, 0.1)) # internal margins
par(mar = c(3.5, 3.75, 3.1, 1)) # external margins

#################### Fig 1.B - Histogram with number of samples in which an interaction is observed #####################

b=barplot(as.matrix(pc_nb_celltypes_matrix), beside=T, xlab='',
          names=rep("", dim(pc_nb_celltypes_matrix)[2]), ylim=c(0,80), space=c(0.4,1),
          ylab="", border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
          lwd=1.5,  mgp=c(3, 0.75, 0), cex.axis=1.1, las=2)

mtext(colnames(nb_celltypes_matrix), at=apply(b, 2, mean), side=1, line=0.5, cex=0.75)

## axis labels
mtext("number of cell types", side=1, line=2.25, cex=0.8)
mtext("% of interactions", side=2, line=2.5, cex=0.8)

## legend & plot label
legend("topright", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n', density=dataset.density, angle=dataset.angle,
       cex=1.1, inset=c(0.05, -0.1), xpd=NA)
mtext("b", side=3, line=1, at=-3.9, font=2, cex=1.2)

################################################################################################

#################### Fig 1.C - Distribution of number of cell types according to distance #####################

ylim=c(0, max(c(as.numeric(mean_nb_celltypes_dist["Original",]), as.numeric(mean_nb_celltypes_dist["Simulated",]))))
ylim[2]=ylim[2]+1

plot(as.numeric(mean_dist["Original",]), as.numeric(mean_nb_celltypes_dist["Original",]), type="l", col=dataset.colors["Original"], ylim=ylim, xlab="", ylab="", axes=F)
lines(as.numeric(mean_dist["Simulated",]), as.numeric(mean_nb_celltypes_dist["Simulated",]), col=dataset.colors["Simulated"], lwd=1.5)

## X axis
xax=pretty(range(as.numeric(mean_dist)))
labels=xax/1e6
axis(side=1, at=xax, labels=labels, mgp=c(3, 0.65, 0), cex.axis=1.1)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1, las=2)

## axis labels
mtext("distance between interacting fragments (Mb)", side=1, line=2.25, cex=0.8)
mtext("mean number of cell types", side=2, line=2.5, cex=0.8)

## confidence intervals

for(dataset in rownames(mean_dist)){
  segments(as.numeric(mean_dist[dataset,]), as.numeric(dist_conf_low_celltypes[dataset,]),  as.numeric(mean_dist[dataset,]), as.numeric(dist_conf_high_celltypes[dataset,]), col=dataset.colors[dataset])
}

## legend & plot label

legend("topright", legend=c("PCHi-C data", "simulated data"), col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1.1, inset=c(0.05, -0.1), xpd=NA)
mtext("c", side=3, line=1, at=-3.15e5, font=2, cex=1.2)

###########################################################################################

dev.off()

###########################################################################################
