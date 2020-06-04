#########################################################################################################################

objects=ls()

if(! "path"%in%objects){
  ## basic parameters
  
  path <- "/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
  pathFigures=paste(path, "Figures/", sep="")

  ## function for genome browser plots

  source("../genome_browser/plot.annotations.R")
  source("../genome_browser/plot.interactions.R")
  
  sp = "human"

  ## colors for the datasets

  dataset.colors=c("forestgreen", "firebrick1")
  names(dataset.colors)=c("Original", "Simulated")

  ## minimum and maximum distance for considered interactions

  minDistance=25e3
  maxDistance=2.5e6
  
  load=T
  prepare=T
}

#########################################################################################################################

## load all data for the figure

if(load){
  obs <- read.table(paste(path, "SupplementaryDataset1/", sp, "/all_interactions.txt", sep=""), header=T)
  simul <- read.table(paste(path, "SupplementaryDataset2/", sp, "/simulated_all_interactions.txt", sep=""), header=T)

  load("RData/data.Shh.figure.RData")
  
  load=FALSE ## we only do it once
}

#########################################################################################################################

## prepare data for figure

if(prepare){
  samples=colnames(obs)[9:dim(obs)[2]] ## first 8 columns contain other info for interactions
  
  print(paste("there are", length(samples), "samples"))

  obs$nb_samples <- apply(obs[,samples], 1, function(x) sum(!is.na(x)))
  obs$sample_class <- cut(obs$nb_samples, breaks=c(0, 1, 5, 10, 15, 20, length(samples)), include.lowest = T)
  obs$dist_class <- cut(obs$distance, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)
  
  simul$nb_samples <- apply(simul[,samples], 1, function(x) sum(!is.na(x)))
  simul$sample_class <- cut(simul$nb_samples, breaks=c(0, 1,  5, 10, 15, 20, length(samples)), include.lowest = T)
  simul$dist_class <- cut(simul$distance, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)

  filtered_data <- list("Original"=obs, "Simulated"=simul)
  
  ## select only unbaited interactions, in cis
  filtered_data <- lapply(filtered_data, function(x) x[which(x$type == "unbaited" & as.character(x$chr_bait) == as.character(x$chr)),])
  
  ## select only interactions found within the selected distance range
  filtered_data <- lapply(filtered_data, function(x) x[which(x$distance <= maxDistance & x$distance >= minDistance),])
  
  ## compute number of interactions in each nb samples class
  
  nb_samples_matrix <- sapply(filtered_data, function(x) as.numeric(table(x$sample_class)))
  rownames(nb_samples_matrix)=c("1", "2-5", "6-10", "11-15", "16-20", paste("21-", length(samples), sep=""))

  nb_samples_matrix=t(nb_samples_matrix)
  pc_nb_samples_matrix=100*nb_samples_matrix/apply(nb_samples_matrix,1, sum)

  ## divide interactions based on distance, compute mean number of samples by distance class
  
  mean_nb_samples_dist <- t(sapply(filtered_data, function(x)   tapply(x$nb_samples, as.factor(x$dist_class), mean, na.rm=T)))
  mean_dist <- t(sapply(filtered_data, function(x)   tapply(x$distance, as.factor(x$dist_class), mean, na.rm=T)))

  dist_conf_low <- t(sapply(filtered_data, function(x) tapply(x$nb_samples, as.factor(x$dist_class), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
  dist_conf_high <- t(sapply(filtered_data, function(x) tapply(x$nb_samples, as.factor(x$dist_class), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))
                

  prepare=FALSE
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(paste(pathFigures, "Figure1.pdf", sep=""), width=6.85, height=7)


## layout

m=matrix(rep(NA, 50*10), nrow=50)

for(i in 1:2){
  m[i,]=rep(1,10)
}

for(i in 3:25){
  m[i,]=rep(2, 10)
}

for(i in 26:50){
  m[i,]=c(rep(3, 5), rep(4, 5))
}

layout(m)

#########################################################################################################################

## annotations in the Shh region

par(mar=c(0.5, 2.75, 0.1, 1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=shhxlim, ylim=c(0,1), xaxs="i", yaxs="i")

#########################################################################################################################

## interactions in the Shh region
par(mar=c(0.5, 2.75, 0.1, 1)) ## left and right margin should be the same as above
plot(1, type="n", xlab="", ylab="", axes=F, xlim=shhxlim, ylim=c(0,1), xaxs="i", yaxs="i")

#########################################################################################################################

par(mai = c(1, 0.8, 0.5, 0.1)) # internal margins
par(mar = c(3.5, 3.75, 2.1, 1)) # external margins

#################### Fig 1.B - Histogram with number of samples in which an interaction is observed #####################

b=barplot(as.matrix(pc_nb_samples_matrix), beside=T, xlab='', names=rep("", dim(pc_nb_samples_matrix)[2]), ylim=c(0,100), space=c(0.4,1), ylab="", border=dataset.colors[c("Original", "Simulated")], col="white", lwd=1.5,  mgp=c(3, 0.75, 0), cex.axis=1.1)

mtext(colnames(nb_samples_matrix), at=apply(b, 2, mean), side=1, line=0.5, cex=0.75)

## axis labels
mtext("number of samples", side=1, line=2.25, cex=0.9)
mtext("% of interactions", side=2, line=2.5, cex=0.9)

## legend

legend("topright", legend=c("original", "simulated"), border=dataset.colors[c("Original", "Simulated")],fill="white", bty='n', cex=1.1, inset=c(0.05, -0.1), xpd=NA)

## plot label
mtext("B", side=3, line=1, at=-2.9, font=2, cex=1.1)

#########################################################################################################################

#################### Fig 1.C - Distribution of number of samples according to distance #####################

plot(as.numeric(mean_dist["Original",]), as.numeric(mean_nb_samples_dist["Original",]), type="l", col=dataset.colors["Original"], ylim=c(0,6), xlab="", ylab="", mgp=c(3, 0.75, 0), xaxt="n", cex.axis=1.1)
lines(as.numeric(mean_dist["Simulated",]), as.numeric(mean_nb_samples_dist["Simulated",]), col=dataset.colors["Simulated"], lwd=1.5)

## X axis
xax=pretty(range(as.numeric(mean_dist)))
labels=xax/1e6
axis(side=1, at=xax, labels=labels, mgp=c(3, 0.65, 0), cex.axis=1.1)

## axis labels
mtext("distance between interacting fragments (Mb)", side=1, line=2.25, cex=0.9)
mtext("mean number of samples", side=2, line=2.5, cex=0.9)

## confidence intervals

for(dataset in rownames(mean_dist)){
  segments(as.numeric(mean_dist[dataset,]), as.numeric(dist_conf_low[dataset,]),  as.numeric(mean_dist[dataset,]), as.numeric(dist_conf_high[dataset,]), col=dataset.colors[dataset])
}

## plot label
mtext("C", side=3, line=1, at=-3.9e5, font=2, cex=1.2)

#########################################################################################################################

dev.off()

#########################################################################################################################
