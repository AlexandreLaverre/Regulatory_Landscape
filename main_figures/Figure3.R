#########################################################################################################################

objects=ls()

if(!"pathFigures"%in%objects){

 source("parameters.R")

 library(ape)
 library(vioplot)

 load=T
 prepare=T
}

#########################################################################################################################

if(load){
 ref_sp = "human"
 close_sp= "macaque"
 target_sp = "mouse"

 enhancers = enhancer.datasets[[ref_sp]]

 selenh="ENCODE"

 load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", ref_sp, ".RData", sep=""))


 ## enhancer statistics
 load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
 
 enh.stats.obs=enhancer.statistics[[ref_sp]][[selenh]][["original"]]
 enh.stats.sim=enhancer.statistics[[ref_sp]][[selenh]][["simulated"]]

 ## enhancer alignment
 align_enhancer_obs=list_align_enh[[selenh]][["enh_align_obs"]]
 align_enhancer_sim=list_align_enh[[selenh]][["enh_align_simul"]]
 
 ## fragment statistics
 load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))
 
 frag.stats.obs=fragment.statistics[[ref_sp]][["original"]]
 frag.stats.sim=fragment.statistics[[ref_sp]][["simulated"]]
 
 load=F
}

#########################################################################################################################

if(prepare){
 
  ## add nb genes 

  frag_align_obs$nb_genes_500kb=frag.stats.obs[frag_align_obs$ID, "nb_genes_500kb"]
  frag_align_simul$nb_genes_500kb=frag.stats.sim[frag_align_simul$ID, "nb_genes_500kb"]
  
  align_enhancer_obs$nb_genes_500kb=enh.stats.obs[align_enhancer_obs$ID, "nb_genes_500kb"]
  align_enhancer_sim$nb_genes_500kb=enh.stats.sim[align_enhancer_sim$ID, "nb_genes_500kb"]

  ## add nb genes class

  frag_align_obs$class_genes_500kb=cut(frag_align_obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag_align_obs$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))

  frag_align_simul$class_genes_500kb=cut(frag_align_simul$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(frag_align_simul$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))

  align_enhancer_obs$class_genes_500kb=cut(align_enhancer_obs$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(align_enhancer_obs$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))
  
  align_enhancer_sim$class_genes_500kb=cut(align_enhancer_sim$nb_genes_500kb, breaks=c(seq(from=0, to=30, by=5), max(align_enhancer_sim$nb_genes_500kb)), include.lowest=T, labels=c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", ">30"))

  ## add repeats
  
  frag.stats.obs$pcrepeat=100*frag.stats.obs$repeat_bp/frag.stats.obs$length
  frag.stats.sim$pcrepeat=100*frag.stats.sim$repeat_bp/frag.stats.sim$length

  frag_align_obs$pcrepeat=frag.stats.obs[frag_align_obs$ID, "pcrepeat"]
  frag_align_simul$pcrepeat=frag.stats.sim[frag_align_simul$ID, "pcrepeat"]
  
  enh.stats.obs$pcrepeat=100*enh.stats.obs$repeat_bp/enh.stats.obs$length
  enh.stats.sim$pcrepeat=100*enh.stats.sim$repeat_bp/enh.stats.sim$length
  
  align_enhancer_obs$pcrepeat=enh.stats.obs[align_enhancer_obs$ID, "pcrepeat"]
  align_enhancer_sim$pcrepeat=enh.stats.sim[align_enhancer_sim$ID, "pcrepeat"]

  prepare=FALSE
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#########################################################################################################################

pdf(paste(pathFigures, "Figure3.pdf", sep=""), width=6.85, height=9)

par(mai = c(0.5, 0.1, 0.3, 0.1)) #bottom, left, top and right

m=matrix(rep(NA, 7*6), nrow=7)

for(i in 1:3){
 m[i,]=c(1,1,2,2,3,3)
}

for(i in 4:5){
 m[i,]=c(4,4,4,5,5,5)
}

for(i in 6:7){
 m[i,]=c(6,6,6,7,7,7)
}

layout(m)

################################## a - Phylogenetic tree ################################################################

tree <- read.tree(paste(pathFigures, "RData/Ensembl_species_tree", sep=""))

tree <- keep.tip(tree, c("Mus_musculus", "Homo_sapiens", "Rattus_norvegicus", "Macaca_mulatta", "Oryctolagus_cuniculus", "Canis_lupus_familiaris", "Bos_taurus", "Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus"))
species <- c("macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
species_names <- c("human", species)

par(mar=c(4.1,1.1, 2.1, 1))
plot(tree, cex=1.2, y.lim=c(0.3,10.5), x.lim=c(0,1.07), label.offset = 0.01, show.tip.label = F, main="")
tiplabels(species_names, bg = NA, adj = -0.1, frame="none", cex=1.3)

 # legend for the plot

legend("bottomleft", fill=dataset.colors, border=dataset.colors, legend = c("PCHi-C data", "simulated data"), bty='n', cex=1.3, xpd=T, inset=c(-0.01, -0.1), horiz=FALSE)

# label
mtext("a", side=3, line=0.5, at=-0.05, font=2, cex=1.2)

######################## b - Restriction fragments sequence conservation ########################

ylim=c(-2, 38.5)
xlim=c(0, 100)

plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=xlim, main="", bty="n")

# Simulated
ypos.sim=c(4,8,12,16,20,24,28,32,36)-0.27

par(bty="n")

vioplot(100*frag_align_simul[,species], at=ypos.sim, add=T, border=dataset.colors["Simulated"], col=rgb(t(col2rgb(dataset.colors["Simulated"])/255), alpha = 0.6), plotCentre="line", axes=F, xaxt="n", yaxt="n", horizontal = T, las=1, cex.main = 1.2, main="", bty="n")

# Original
ypos.obs=c(5,9,13,17,21,25,29,33,37)+0.27

vioplot(100*frag_align_obs[,species], at=ypos.obs, add=T, axes=F, xaxt="n", yaxt="n", horizontal = T, border=dataset.colors["Original"], col=rgb(t(col2rgb(dataset.colors["Original"])/255), alpha = 0.6), plotCentre="line")

## add mean point

points(x=100*apply(frag_align_simul[,species], 2, mean, na.rm=T), y = ypos.sim, col = "white", pch=20, cex=0.8)
points(x=100*apply(frag_align_obs[,species], 2, mean, na.rm=T), y = ypos.obs, col = "white", pch=20, cex=0.8)

## axis and legend

axis(1, pos=0.7, at=seq(0,100,20), labels=c("0", "20", "40", "60", "80", "100"), cex.axis=1.2)
mtext("% aligned sequence", side=1, xpd = TRUE, cex=0.8, line=0.5)

mtext("restriction fragments", side=3, line=-1, cex=0.8)

mtext("b", side=3, line=0.5, at=-8, font=2, cex=1.2)

########################## c - ENCODE enhancers sequence conservation ########################


ylim=c(-2, 38.5)
xlim=c(0, 100)

plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=xlim, main="", bty="n")

# simulated enhancer
vioplot(100*align_enhancer_sim[,species], at=ypos.sim, col=rgb(t(col2rgb(dataset.colors["Simulated"])/255), alpha = 0.6), border=dataset.colors["Simulated"], add=T, axes=F, xaxt="n", yaxt="n", horizontal = T, las=1, cex.main = 1.2, main="", plotCentre="line")

# original enhancer
vioplot(100*align_enhancer_obs[,species],at=ypos.obs, col=rgb(t(col2rgb(dataset.colors["Original"])/255), alpha = 0.6), border=dataset.colors["Original"], add=T, axes=F, xaxt="n", yaxt="n", horizontal = T, plotCentre="line")

# add mean point
points(x = apply(100*align_enhancer_sim[,species], 2, mean, na.rm=T), y=ypos.sim, col = "white", pch=20, cex=0.8)
points(x = apply(100*align_enhancer_obs[,species], 2, mean, na.rm=T), y = ypos.obs, col = "white", pch=20, cex=0.8)

# axis and legend
axis(1, pos=0.7, at=seq(0,100,20), labels=c("0", "20", "40", "60", "80", "100"), cex.axis=1.2)

mtext("% aligned sequence", side=1, xpd = TRUE, cex=0.8, line=0.5)

mtext("enhancers", side=3, line=-1, cex=0.8)

mtext("c", side=3, line=0.5, at=-8, font=2, cex=1.2)

######################## Conserved sequence human to  mouse vs distance to promoters  ########################

par(mai = c(0.8, 0.6, 0.2, 0.2)) #bottom, left, top and right
par(mar=c(4.1, 4.5, 2, 1.5))

nbclasses=length(levels( frag_align_obs$dist_class))
xpos=1:nbclasses

xlim=c(-0.5, max(xpos)+1)

## axis position
class_leg <- c("0", "0.5", "1", "1.5", "2")
xax=seq(from=0, to=max(xpos)+1, by=10)

labels=c("d", "e")
names(labels)=c("restriction fragments", "enhancers")

other_sp=target_sp

for(type in c("restriction fragments", "enhancers")){
  
  if(type=="restriction fragments"){
    data.obs=frag_align_obs
    data.sim=frag_align_simul
  }

  if(type=="enhancers"){
    data.obs=align_enhancer_obs
    data.sim=align_enhancer_sim
  }
  
  mean.val.obs=tapply(100*data.obs[, other_sp], data.obs$dist_class, function(x) mean(x, na.rm=T))
  ci.low.obs=tapply(100*data.obs[, other_sp], data.obs$dist_class, function(x) t.test(x)[["conf.int"]][1])
  ci.high.obs=tapply(100*data.obs[, other_sp], data.obs$dist_class, function(x) t.test(x)[["conf.int"]][2])
  
  mean.val.sim=tapply(100*data.sim[, other_sp], data.sim$dist_class, function(x) mean(x, na.rm=T))
  ci.low.sim=tapply(100*data.sim[, other_sp], data.sim$dist_class, function(x) t.test(x)[["conf.int"]][1])
  ci.high.sim=tapply(100*data.sim[, other_sp], data.sim$dist_class, function(x) t.test(x)[["conf.int"]][2])

 ylim=range(c(ci.low.obs, ci.high.obs, ci.low.sim, ci.high.sim))

 dy=diff(ylim)/20
 ylim=ylim+c(-dy, dy)

 plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")

 points(xpos, mean.val.obs, col=dataset.colors["Original"], pch=20)
 segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])

 points(xpos, mean.val.sim, col=dataset.colors["Simulated"], pch=20)
 segments(xpos, ci.low.sim, xpos, ci.high.sim, col=dataset.colors["Simulated"])

  axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)

  if(type=="enhancers"){
    mtext("distance to promoters (Mb)", side=1, line=2.2, cex=0.8)
  }

  if(type=="restriction fragments"){
    mtext("distance to baits (Mb)", side=1, line=2.2, cex=0.8)
  }
  
 axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
 mtext("% aligned sequence", side=2, line=3, cex=0.8)

 mtext(paste(ref_sp, " vs. ", other_sp, ", ", type,sep=""), side=3, cex=0.8, line=1)

 mtext(labels[type], side=3, line=1, at=-7.75, font=2, cex=1.2)
}

#######################################################################################################

## sequence conservation for gene classes

par(mai = c(0.8, 0.6, 0.2, 0.2)) #bottom, left, top and right

par(mar=c(5.1, 4.5, 1, 1.5))

nbclasses=length(levels(frag_align_obs$class_genes_500kb))
xpos=1:nbclasses

xlim=c(0.5, max(xpos)+0.5)

smallx=c(-0.1, 0.1)
names(smallx)=c("obs", "sim")

## axis position

xax=xpos
class_leg=levels(frag_align_obs$class_genes_500kb)

labels=c("f", "g")
names(labels)=c("restriction fragments", "enhancers")

other_sp=target_sp

for(type in c("restriction fragments", "enhancers")){
  
  if(type=="restriction fragments"){
    data.obs=frag_align_obs
    data.sim=frag_align_simul
  }

  if(type=="enhancers"){
    data.obs=align_enhancer_obs
    data.sim=align_enhancer_sim
  }

  ## all 
  mean.val.obs=tapply(100*data.obs[, other_sp], data.obs$class_genes_500kb, function(x) mean(x, na.rm=T))
  ci.low.obs=tapply(100*data.obs[, other_sp], data.obs$class_genes_500kb, function(x) t.test(x)[["conf.int"]][1])
  ci.high.obs=tapply(100*data.obs[, other_sp], data.obs$class_genes_500kb, function(x) t.test(x)[["conf.int"]][2])

  mean.val.sim=tapply(100*data.sim[, other_sp], data.sim$class_genes_500kb, function(x) mean(x, na.rm=T))
  ci.low.sim=tapply(100*data.sim[, other_sp], data.sim$class_genes_500kb, function(x) t.test(x)[["conf.int"]][1])
  ci.high.sim=tapply(100*data.sim[, other_sp], data.sim$class_genes_500kb, function(x) t.test(x)[["conf.int"]][2])

  ## no repeats
  
  mean.val.obs.norep=tapply(100*data.obs[which(data.obs$pcrepeat==0), other_sp], data.obs[which(data.obs$pcrepeat==0), "class_genes_500kb"], function(x) mean(x, na.rm=T))
  ci.low.obs.norep=tapply(100*data.obs[which(data.obs$pcrepeat==0), other_sp], data.obs[which(data.obs$pcrepeat==0), "class_genes_500kb"], function(x) t.test(x)[["conf.int"]][1])
  ci.high.obs.norep=tapply(100*data.obs[which(data.obs$pcrepeat==0), other_sp], data.obs[which(data.obs$pcrepeat==0), "class_genes_500kb"], function(x) t.test(x)[["conf.int"]][2])

  mean.val.sim.norep=tapply(100*data.sim[which(data.sim$pcrepeat==0), other_sp], data.sim[which(data.sim$pcrepeat==0), "class_genes_500kb"], function(x) mean(x, na.rm=T))
  ci.low.sim.norep=tapply(100*data.sim[which(data.sim$pcrepeat==0), other_sp], data.sim[which(data.sim$pcrepeat==0), "class_genes_500kb"], function(x) t.test(x)[["conf.int"]][1])
  ci.high.sim.norep=tapply(100*data.sim[which(data.sim$pcrepeat==0), other_sp], data.sim[which(data.sim$pcrepeat==0), "class_genes_500kb"], function(x) t.test(x)[["conf.int"]][2])

  ## plot
  
  ylim=range(c(ci.low.obs, ci.high.obs, ci.low.sim, ci.high.sim, ci.high.sim.norep, ci.high.obs.norep))
  
  dy=diff(ylim)/20
  ylim=ylim+c(-dy, dy)
  
  plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")
  
  points(xpos+smallx["obs"], mean.val.obs, col=dataset.colors["Original"], pch=20, cex=1.35)
  segments(xpos+smallx["obs"], ci.low.obs, xpos+smallx["obs"], ci.high.obs, col=dataset.colors["Original"])

  points(xpos+smallx["sim"], mean.val.sim, col=dataset.colors["Simulated"], pch=20, cex=1.35)
  segments(xpos+smallx["sim"], ci.low.sim, xpos+smallx["sim"], ci.high.sim, col=dataset.colors["Simulated"])

  points(xpos+smallx["obs"], mean.val.obs.norep, bg="white", col=dataset.colors["Original"], pch=21)
  segments(xpos+smallx["obs"], ci.low.obs.norep, xpos+smallx["obs"], ci.high.obs.norep, col=dataset.colors["Original"])
  
  points(xpos+smallx["sim"], mean.val.sim.norep, bg="white", col=dataset.colors["Simulated"], pch=21)
  segments(xpos+smallx["sim"], ci.low.sim.norep, xpos+smallx["sim"], ci.high.sim.norep, col=dataset.colors["Simulated"])

  ## axes
  axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1, las=2)
  mtext("number of genes within 500kb", side=1, line=3.85, cex=0.8)
  
  axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
  mtext("% aligned sequence", side=2, line=3, cex=0.8)
  
  abline(v=xpos[-length(xpos)]+0.5, lty=2, col="gray40")
  
  mtext(labels[type], side=3, line=1, at=-0.75, font=2, cex=1.2)

  if(type=="enhancers"){
    legend("topright", box.col="white", bg="white", pch=21, pt.bg=c("black", "white"), legend=c("all data", "without repeats"),xpd=NA, inset=c(0.01, -0.05), cex=1.1)
  }
    
}


#######################################################################################################

dev.off()

#######################################################################################################
