library(ape)
library(vioplot)
setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")

path <- "/home/laverre/Data/Regulatory_landscape/result/Main_figures/"
pdf(paste(path, "/Figure3_A_sequence_conservation.pdf", sep=""), width=7, height=5)
par(mai = c(0.3, 0.2, 0.3, 0.3))
a <- matrix(c(1,1,2,2,3,3),ncol=6, nrow = 2, byrow=T)
b <- matrix(c(3,3,3,4,4,4), nrow=1, byrow=F)
c <- rbind(a,b)
layout(c)

tree <- read.tree("/home/laverre/Documents/Regulatory_Landscape/data/ensembl_tree")
tree <- keep.tip(tree, c("Mus_musculus", "Homo_sapiens", "Rattus_norvegicus", "Macaca_mulatta", "Oryctolagus_cuniculus", "Canis_lupus_familiaris", "Bos_taurus", "Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus"))

plot(tree, cex=1.2, y.lim=c(0.4,10.3), x.lim=c(0,1.07), label.offset = 0.01, show.tip.label = F, main="Fig3.A")
tiplabels(c("human", "macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken"), bg = NA, adj = -0.1, frame="none", cex=1.3)
tiplabels(c("human", "macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken"), bg = NA, adj = -0.1, frame="none", cex=1.3)

legend(y=0.7, x=0, fill=c("forestgreen","dodgerblue", "firebrick1"), 
       legend = c("Original", "Original with enh", "Simulated"),
       text.width=c(0.3), ncol=2, bty='n', cex=1.1)

load("Fig4_human_mean.Rdata")

vioplot(c(0,conserv_simul), at=c(0,4,8,12,16,20,24,28,32,36), names=c("human", colnames(conserv)), border="firebrick1",
        col='firebrick1', axes=F, yaxt='n', horizontal = T, las=1, cex.names = 0.8, main="Human conservation")

vioplot(c(0,conserv), at=c(1,5,9,13,17,21,25,29,33,37), col='forestgreen', add=T, axes=F, horizontal = T, border="forestgreen")
vioplot(c(0,conserv_enh), at=c(2,6,10,14,18,22,26,30,34,38), col='dodgerblue', add=T, axes=F, horizontal = T, border="dodgerblue")
axis(1, at=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))

load("Fig4_mouse_mean.Rdata")

vioplot(c(conserv_simul[,c(3,4)],0,conserv_simul[,c(1,2,5,6,7,8,9)]), at=c(0,4,8,12,16,20,24,28,32,36), col='firebrick1', border="firebrick1", axes=F, yaxt='n', horizontal = T, las=1, cex.names = 0.8, main="Mouse conservation")
vioplot(c(conserv[,c(3,4)],0,conserv[,c(1,2,5,6,7,8,9)]), at=c(1,5,9,13,17,21,25,29,33,37), col='forestgreen', border="forestgreen", add=T, axes=F, horizontal = T)
vioplot(c(conserv_enh[,c(3,4)],0,conserv_enh[,c(1,2,5,6,7,8,9)]), at=c(2,6,10,14,18,22,26,30,34,38), col='dodgerblue', border="dodgerblue", add=T, axes=F, horizontal = T)

axis(1, at=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))

#### Fig 4.B Human to mouse vs distance to promoters
plot(conserv_dist_all_sp[["mouse"]][["obs"]][1:50,"inter"], type="l", col="forestgreen", cex=CEX_lines, main=paste(sp_origin, " to mouse", sep=""),
     xlab="", ylab="Ungapped Non-exonic Score", xaxt = "n", ylim=c(xmin,xmax), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)

for (row in 1:nrow(conserv_dist_all_sp[["mouse"]][["obs"]])){
  segments(x0=row,y0=conserv_dist_all_sp[["mouse"]][["obs"]][row,"int_start"],x1=row,y1=conserv_dist_all_sp[["mouse"]][["obs"]][row,"int_end"], col='forestgreen', lwd=0.3)}

lines(conserv_dist_all_sp[["mouse"]][["simul"]][1:50,"inter"], type="l", col="firebrick1", cex=CEX_lines)
for (row in 1:nrow(conserv_dist_all_sp[["mouse"]][["simul"]])){
  segments(x0=row,y0=conserv_dist_all_sp[["mouse"]][["simul"]][row,"int_start"],x1=row,y1=conserv_dist_all_sp[["mouse"]][["simul"]][row,"int_end"], col='firebrick1', lwd=0.3)}

lines(conserv_dist_all_sp[["mouse"]][["obs_enh"]][1:50,"inter"], type="l", col="dodgerblue", cex=CEX_lines)
for (row in 1:nrow(conserv_dist_all_sp[["mouse"]][["obs_enh"]])){
  segments(x0=row,y0=conserv_dist_all_sp[["mouse"]][["obs_enh"]][row,"int_start"],x1=row,y1=conserv_dist_all_sp[["mouse"]][["obs_enh"]][row,"int_end"], col='dodgerblue', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.02, class_leg, xpd = TRUE, cex=CEX)
dev.off()
