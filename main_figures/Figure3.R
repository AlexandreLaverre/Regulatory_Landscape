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

 load(paste(pathFigures, "RData/data.sequence.conservation.pcungapped.", ref_sp, ".Rdata", sep=""))

 load=F
}

#########################################################################################################################

if(prepare){
  align_enhancer_obs=list_align_enh[[selenh]][["enh_align_obs"]]
  align_enhancer_simul=list_align_enh[[selenh]][["enh_align_simul"]]

  prepare=F
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
 m[i,]=c(4,4,4,6,6,6)
}

for(i in 6:7){
 m[i,]=c(5,5,5,7,7,7)
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

legend("bottomleft", fill=dataset.colors, border=dataset.colors, legend = c("PCHi-C data", "simulated data"), bty='n', cex=1.3, xpd=T, inset=c(-0.01, -0.05), horiz=FALSE)

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

points(x=100*apply(frag_align_simul[,species], 2, mean), y = ypos.sim, col = "white", pch=20, cex=0.8)
points(x=100*apply(frag_align_obs[,species], 2, mean), y = ypos.obs, col = "white", pch=20, cex=0.8)

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
vioplot(100*align_enhancer_simul[,species], at=ypos.sim, col=rgb(t(col2rgb(dataset.colors["Simulated"])/255), alpha = 0.6), border=dataset.colors["Simulated"], add=T, axes=F, xaxt="n", yaxt="n", horizontal = T, las=1, cex.main = 1.2, main="", plotCentre="line")

# original enhancer
vioplot(100*align_enhancer_obs[,species],at=ypos.obs, col=rgb(t(col2rgb(dataset.colors["Original"])/255), alpha = 0.6), border=dataset.colors["Original"], add=T, axes=F, xaxt="n", yaxt="n", horizontal = T, plotCentre="line")

# add mean point
points(x = apply(100*align_enhancer_simul[,species], 2, mean, na.rm=T), y=ypos.sim, col = "white", pch=20, cex=0.8)
points(x = apply(100*align_enhancer_obs[,species], 2, mean, na.rm=T), y = ypos.obs, col = "white", pch=20, cex=0.8)

# axis and legend
axis(1, pos=0.7, at=seq(0,100,20), labels=c("0", "20", "40", "60", "80", "100"), cex.axis=1.2)

mtext("% aligned sequence", side=1, xpd = TRUE, cex=0.8, line=0.5)

mtext("enhancers", side=3, line=-1, cex=0.8)

mtext("c", side=3, line=0.5, at=-8, font=2, cex=1.2)

######################## d - Conserved sequence human to macaque & mouse vs distance to promoters, restriction fragments ########################

par(mai = c(0.8, 0.6, 0.2, 0.2)) #bottom, left, top and right
par(mar=c(4.1, 4.5, 2, 1.5))

nbclasses=length(levels( frag_align_obs$dist_class))
xpos=1:nbclasses

xlim=c(-0.5, max(xpos)+1)

## axis position
class_leg <- c("0", "0.5", "1", "1.5", "2")
xax=seq(from=0, to=max(xpos)+1, by=10)

labels=c("d", "f")
names(labels)=c(close_sp, target_sp)

for(other_sp in c(close_sp, target_sp)){
 mean.val.obs=tapply(100*frag_align_obs[, other_sp], frag_align_obs$dist_class, function(x) mean(x, na.rm=T))
 ci.low.obs=tapply(100*frag_align_obs[, other_sp], frag_align_obs$dist_class, function(x) t.test(x)[["conf.int"]][1])
 ci.high.obs=tapply(100*frag_align_obs[, other_sp], frag_align_obs$dist_class, function(x) t.test(x)[["conf.int"]][2])

 mean.val.simul=tapply(100*frag_align_simul[, other_sp], frag_align_simul$dist_class, function(x) mean(x, na.rm=T))
 ci.low.simul=tapply(100*frag_align_simul[, other_sp], frag_align_simul$dist_class, function(x) t.test(x)[["conf.int"]][1])
 ci.high.simul=tapply(100*frag_align_simul[, other_sp], frag_align_simul$dist_class, function(x) t.test(x)[["conf.int"]][2])

 ylim=range(c(ci.low.obs, ci.high.obs, ci.low.simul, ci.high.simul))

 dy=diff(ylim)/20
 ylim=ylim+c(-dy, dy)

 plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")

 lines(xpos, mean.val.obs, col=dataset.colors["Original"])
 segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])

 lines(xpos, mean.val.simul, col=dataset.colors["Simulated"])
 segments(xpos, ci.low.simul, xpos, ci.high.simul, col=dataset.colors["Simulated"])

 axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
 mtext("distance to promoters (Mb)", side=1, line=2.2, cex=0.8)

 axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
 mtext("% aligned sequence", side=2, line=3, cex=0.8)

 mtext(paste(ref_sp, " vs. ", other_sp, ", restriction fragments",sep=""), side=3, cex=0.8)

 mtext(labels[other_sp], side=3, line=1, at=-7.75, font=2, cex=1.2)
}

#######################################################################################################

## same, for enhancers

par(mai = c(0.8, 0.6, 0.2, 0.2)) #bottom, left, top and right

par(mar=c(4.1, 4.5, 2, 1.5))

nbclasses=length(levels( frag_align_obs$dist_class))
xpos=1:nbclasses

xlim=c(-0.5, max(xpos)+1)

## axis position
class_leg <- c("0", "0.5", "1", "1.5", "2")
xax=seq(from=0, to=max(xpos)+1, by=10)

labels=c("e", "g")
names(labels)=c(close_sp, target_sp)

for(other_sp in c(close_sp, target_sp)){
 mean.val.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) mean(x, na.rm=T))
 ci.low.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) t.test(x)[["conf.int"]][1])
 ci.high.obs=tapply(100*align_enhancer_obs[, other_sp], align_enhancer_obs$dist_class, function(x) t.test(x)[["conf.int"]][2])

 mean.val.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) mean(x, na.rm=T))
 ci.low.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) t.test(x)[["conf.int"]][1])
 ci.high.simul=tapply(100*align_enhancer_simul[, other_sp], align_enhancer_simul$dist_class, function(x) t.test(x)[["conf.int"]][2])

 ylim=range(c(ci.low.obs, ci.high.obs, ci.low.simul, ci.high.simul))

 dy=diff(ylim)/20
 ylim=ylim+c(-dy, dy)

 plot(1, type="n", xlab="", ylab="", axes=F, main="", xlim=xlim, ylim=ylim, xaxs="i")

 lines(xpos, mean.val.obs, col=dataset.colors["Original"])
 segments(xpos, ci.low.obs, xpos, ci.high.obs, col=dataset.colors["Original"])

 lines(xpos, mean.val.simul, col=dataset.colors["Simulated"])
 segments(xpos, ci.low.simul, xpos, ci.high.simul, col=dataset.colors["Simulated"])

 axis(side=1, at=xax, mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)
 mtext("distance to promoters (Mb)", side=1, line=2.2, cex=0.8)

 axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
 mtext("% aligned sequence", side=2, line=3, cex=0.8)

 mtext(paste(ref_sp, " vs. ", other_sp, ", enhancers", sep=""), side=3, cex=0.8)

 mtext(labels[other_sp], side=3, line=1, at=-7.75, font=2, cex=1.2)
}

#######################################################################################################

dev.off()

#######################################################################################################
