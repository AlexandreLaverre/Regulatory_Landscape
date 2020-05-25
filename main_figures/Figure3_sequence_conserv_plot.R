library(ape)
library(vioplot)
setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")

path <- "/home/laverre/Data/Regulatory_landscape/result/Main_figures/"
pdf(paste(path, "/Figure3_sequence_conservation.pdf", sep=""), width=7, height=10)
par(mai = c(0.5, 0.1, 0.3, 0.1)) #bottom, left, top and right 
a <- matrix(c(1,1,2,2,3,3),ncol=6, nrow = 2, byrow=T)
b <- matrix(c(4,4,4,5,5,5), nrow=1, byrow=F)
c <- matrix(c(6,6,6,7,7,7), nrow=1, byrow=F)
d <- rbind(a,b,c)
layout(d)

######################## A - Phylogenetic tree ######################## 
tree <- read.tree("/home/laverre/Documents/Regulatory_Landscape/data/ensembl_tree")
tree <- keep.tip(tree, c("Mus_musculus", "Homo_sapiens", "Rattus_norvegicus", "Macaca_mulatta", "Oryctolagus_cuniculus", "Canis_lupus_familiaris", "Bos_taurus", "Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus"))


plot(tree, cex=1.2, y.lim=c(0.3,10.3), x.lim=c(0,1.07), label.offset = 0.01, show.tip.label = F, main="A.")
tiplabels(c("human", "macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken"), bg = NA, adj = -0.1, frame="none", cex=1.3)
tiplabels(c("human", "macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken"), bg = NA, adj = -0.1, frame="none", cex=1.3)

par(xpd=TRUE)
legend(y=0.5, x=0, fill=c("firebrick1","dodgerblue", "forestgreen"), 
       legend = c("Simulated", "Original with enhancer", "Original"),
       text.width=c(0.4), ncol=2, bty='n', cex=1.3)

######################## B - Restriction fragments sequence conservation ######################## 
load(paste(path, "Fig3_seq_conserv_human.Rdata", sep=""))
species <- c("macaque", "mouse", "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")


vioplot(c(0,align_simul[,species]), at=c(0,4,8,12,16,20,24,28,32,36), border="firebrick1", col=rgb(t(col2rgb('firebrick1')/255), alpha = 0.6),
        plotCentre="line", axes=F, yaxt='n', horizontal = T, las=1, cex.main = 1.2, main="B. Contacted sequences")

vioplot(c(0,align_obs[,species]), at=c(1,5,9,13,17,21,25,29,33,37), add=T, axes=F, horizontal = T, 
        border="forestgreen", col=rgb(t(col2rgb('forestgreen')/255), alpha = 0.6), plotCentre="line")
vioplot(c(0,align_obs_enh[,species]), at=c(2,6,10,14,18,22,26,30,34,38), add=T, axes=F, horizontal = T,
        border="dodgerblue", col=rgb(t(col2rgb('dodgerblue')/255), alpha = 0.6), plotCentre="line")

axis(1, pos=0.7, at=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))
mtext("Alignment score", side=1, xpd = TRUE, cex=0.8)

######################## C - ENCODE enhancers sequence conservation ######################## 

vioplot(c(0, align_enhancers_simul[,species]), at=c(0,3,6,9,12,15,18,21,24,27), col=rgb(t(col2rgb('firebrick1')/255), alpha = 0.6), border="firebrick1",
        axes=F, yaxt='n', horizontal = T, las=1, cex.main = 1.2, main="C. Contacted ENCODE enhancers", plotCentre="line")

vioplot(c(0, align_enhancers_obs[,species]),
        at=c(1,4,7,10,13,16,19,22,25,28), col=rgb(t(col2rgb('forestgreen')/255), alpha = 0.6), border="forestgreen", add=T, axes=F, horizontal = T, plotCentre="line")

axis(1, pos=0.7, at=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))
mtext("Alignment score", side=1, xpd = TRUE, cex=0.8)

######################## D - Conserved sequence Human to mouse vs distance to promoters ########################
par(mai = c(0.8, 0.6, 0.2, 0.2)) #bottom, left, top and right 
par(xpd=FALSE)
class_leg <- c("25Kb", "500Kb", "1Mb", "1.5Mb", "2Mb")

xmin=0.25
xmax=0.45
CEX=1.2
CEX_lines=1

plot(obs_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="D. Human To Mouse",
     xlab="", ylab="Alignment score", xaxt = "n", ylim=c(xmin,xmax), cex.lab=CEX, cex.axis=CEX)

for (row in 1:nrow(obs_dist)){
  segments(x0=row,y0=obs_dist[row,"int_start"],x1=row,y1=obs_dist[row,"int_end"], col='forestgreen', lwd=0.3)}

lines(simul_dist[,"inter"], type="l", col="firebrick1", cex=CEX_lines)
for (row in 1:nrow(simul_dist)){
  segments(x0=row,y0=simul_dist[row,"int_start"],x1=row,y1=simul_dist[row,"int_end"], col='firebrick1', lwd=0.3)}

lines(obs_enh_dist[,"inter"], type="l", col="dodgerblue", cex=CEX_lines)
for (row in 1:nrow(obs_enh_dist)){
  segments(x0=row,y0=obs_enh_dist[row,"int_start"],x1=row,y1=obs_enh_dist[row,"int_end"], col='dodgerblue', lwd=0.3)}

axis(1, at=seq(1,nrow(obs_dist)+1,10), labels=F)
text(seq(1,nrow(obs_dist)+1,10),par("usr")[3]-0.02, class_leg, xpd = TRUE, cex=CEX)
mtext("Distance to promoters (bp)", side=1, line=2, cex=0.8)
 
######################## E - Conserved sequence Human to Macaque vs distance to promoters ######################## 
xmin=0.78
xmax=0.95

plot(obs_dist_mac[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="E. Human To Macaque",
     xlab="", ylab="Alignment score", xaxt = "n", ylim=c(xmin,xmax), cex.lab=CEX, cex.axis=CEX)

for (row in 1:nrow(obs_dist_mac)){
  segments(x0=row,y0=obs_dist_mac[row,"int_start"],x1=row,y1=obs_dist_mac[row,"int_end"], col='forestgreen', lwd=0.3)}

lines(simul_dist_mac[,"inter"], type="l", col="firebrick1", cex=CEX_lines)
for (row in 1:nrow(simul_dist_mac)){
  segments(x0=row,y0=simul_dist_mac[row,"int_start"],x1=row,y1=simul_dist_mac[row,"int_end"], col='firebrick1', lwd=0.3)}

lines(obs_enh_dist_mac[,"inter"], type="l", col="dodgerblue", cex=CEX_lines)
for (row in 1:nrow(obs_enh_dist_mac)){
  segments(x0=row,y0=obs_enh_dist_mac[row,"int_start"],x1=row,y1=obs_enh_dist_mac[row,"int_end"], col='dodgerblue', lwd=0.3)}

axis(1, at=seq(1,nrow(obs_dist)+1,10), labels=F)
text(seq(1,nrow(obs_dist)+1,10),par("usr")[3]-0.02, class_leg, xpd = TRUE, cex=CEX)
mtext("Distance to promoters (bp)", side=1, line=2.5, cex=0.8)

######################## F - Repeat proportion vs distance to promoters ######################## 
ymin=32
ymax=55

plot(obs_repet_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="F.",
     xlab="", ylab="Repeat proportion (%)", xaxt = "n", ylim=c(ymin,ymax), cex.lab=CEX, cex.axis=CEX, las=2)

for (row in 1:nrow(obs_repet_dist)){
  segments(x0=row,y0=obs_repet_dist[row,"int_start"],x1=row,y1=obs_repet_dist[row,"int_end"], col='forestgreen', lwd=0.3)}

lines(simul_repet_dist[,"inter"], type="l", col="firebrick1", cex=CEX_lines)
for (row in 1:nrow(simul_repet_dist)){
  segments(x0=row,y0=simul_repet_dist[row,"int_start"],x1=row,y1=simul_repet_dist[row,"int_end"], col='firebrick1', lwd=0.3)}

lines(obs_enh_repet_dist[,"inter"], type="l", col="dodgerblue", cex=CEX_lines)
for (row in 1:nrow(obs_enh_repet_dist)){
  segments(x0=row,y0=obs_enh_repet_dist[row,"int_start"],x1=row,y1=obs_enh_repet_dist[row,"int_end"], col='dodgerblue', lwd=0.3)}

axis(1, at=seq(1,nrow(obs_repet_dist)+1,10), labels=F)
text(seq(1,nrow(obs_repet_dist)+1,10),par("usr")[3]-3, class_leg, xpd = TRUE, cex=CEX)
mtext("Distance to promoters (bp)", side=1, line=2.5, cex=0.8)

######################## G - Exonic proportion vs distance to promoters ######################## 
xmin=0
xmax=12

plot(obs_exon_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="G.",
     xlab="", ylab="Exonic proportion (%)", xaxt = "n", ylim=c(xmin,xmax), cex.lab=CEX, cex.axis=CEX, las=2)

for (row in 1:nrow(obs_exon_dist)){
  segments(x0=row,y0=obs_exon_dist[row,"int_start"],x1=row,y1=obs_exon_dist[row,"int_end"], col='forestgreen', lwd=0.3)}

lines(simul_exon_dist[,"inter"], type="l", col="firebrick1", cex=CEX_lines)
for (row in 1:nrow(simul_exon_dist)){
  segments(x0=row,y0=simul_exon_dist[row,"int_start"],x1=row,y1=simul_exon_dist[row,"int_end"], col='firebrick1', lwd=0.3)}

lines(obs_enh_exon_dist[,"inter"], type="l", col="dodgerblue", cex=CEX_lines)
for (row in 1:nrow(obs_enh_exon_dist)){
  segments(x0=row,y0=obs_enh_exon_dist[row,"int_start"],x1=row,y1=obs_enh_exon_dist[row,"int_end"], col='dodgerblue', lwd=0.3)}

axis(1, at=seq(1,nrow(obs_exon_dist)+1,10), labels=F)
text(seq(1,nrow(obs_exon_dist)+1,10),par("usr")[3]-1.5, class_leg, xpd = TRUE, cex=CEX)
mtext("Distance to promoters (bp)", side=1, line=2.5, cex=0.8)

dev.off()
