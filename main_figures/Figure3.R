#########################################################################################################################
source("parameters.R") ## pathFiguress are defined based on the user name

library(ape)
library(vioplot)

ref_sp = "human"
target_sp = "mouse"

load(paste(pathFigures, "Fig3_", ref_sp, "_test.Rdata", sep=""))

enhancers = c("FANTOM5", "ENCODE")
if(ref_sp == "human"){enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")}
if(ref_sp == "human"){pdf_name="Figure3_test.pdf"}else{pdf_name="Sup_Figure11.pdf"}

col <- c("red", "navy", "forestgreen", "orange")
#########################################################################################################################

pdf(paste(pathFigures, pdf_name, sep=""), width=7, height=10)
par(mai = c(0.5, 0.1, 0.3, 0.1)) #bottom, left, top and right 
a <- matrix(c(1,1,2,2,3,3),ncol=6, nrow = 2, byrow=T)
b <- matrix(c(4,4,4,5,5,5), nrow=1, byrow=F)
c <- matrix(c(6,6,6,7,7,7), nrow=1, byrow=F)
d <- rbind(a,b,c)
layout(d)

######################## A - Phylogenetic tree ######################## 
tree <- read.tree(paste(pathFigures, "Ensembl_species_tree", sep=""))
tree <- keep.tip(tree, c("Mus_musculus", "Homo_sapiens", "Rattus_norvegicus", "Macaca_mulatta", "Oryctolagus_cuniculus", "Canis_lupus_familiaris", "Bos_taurus", "Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus"))
species <- c("macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
species_names <- c("human", species)

if(ref_sp == "mouse"){
  tree <- rotate(tree, c(1,5))
  species <- c( "rat", "rabbit", "human", "macaque", "cow", "dog", "elephant", "opossum", "chicken")}

plot(tree, cex=1.2, y.lim=c(0.3,10.3), x.lim=c(0,1.07), label.offset = 0.01, show.tip.label = F, main="")
tiplabels(species_names, bg = NA, adj = -0.1, frame="none", cex=1.3)
tiplabels(species_names, bg = NA, adj = -0.1, frame="none", cex=1.3)

mtext("A", side=3, line=1, at=0, font=2, cex=1.2)

par(xpd=TRUE)
legend(y=0.5, x=0, fill=c("firebrick1","dodgerblue", "forestgreen"), 
       legend = c("Simulated", "Original with enhancer", "Original"),
       text.width=c(0.4), ncol=2, bty='n', cex=1.3)

######################## B - Restriction fragments sequence conservation ######################## 
vioplot(c(0,align_simul[,species]), at=c(0,4,8,12,16,20,24,28,32,36), border="firebrick1", col=rgb(t(col2rgb('firebrick1')/255), alpha = 0.6),
        plotCentre="line", axes=F, yaxt='n', horizontal = T, las=1, cex.main = 1.2, main="")

vioplot(c(0,align_obs[,species]), at=c(1,5,9,13,17,21,25,29,33,37), add=T, axes=F, horizontal = T, 
        border="forestgreen", col=rgb(t(col2rgb('forestgreen')/255), alpha = 0.6), plotCentre="line") 
vioplot(c(0,align_obs_enh[,species]), at=c(2,6,10,14,18,22,26,30,34,38), add=T, axes=F, horizontal = T,
        border="dodgerblue", col=rgb(t(col2rgb('dodgerblue')/255), alpha = 0.6), plotCentre="line") 

points(x=apply(align_simul[,species], 2, mean), y=c(4,8,12,16,20,24,28,32,36), col = "white", pch=20, cex=0.8)

points(x = apply(align_obs[,species], 2, mean), y = c(5,9,13,17,21,25,29,33,37), col = "white", pch=20, cex=0.8)

points(x = apply(align_obs_enh[,species], 2, mean), y = c(6,10,14,18,22,26,30,34,38), col = "white", pch=20, cex=0.8)

axis(1, pos=0.7, at=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.lab=1.2)
mtext("Alignment score", side=1, xpd = TRUE, cex=0.8)
mtext("B", side=3, line=1, at=0.1, font=2, cex=1.2)

######################## C - ENCODE enhancers sequence conservation ######################## 

vioplot(c(0, align_enhancers_simul[,species]), at=c(0,3,6,9,12,15,18,21,24,27), col=rgb(t(col2rgb('firebrick1')/255), alpha = 0.6), border="firebrick1",
        axes=F, yaxt='n', horizontal = T, las=1, cex.main = 1.2, main="", plotCentre="line")

vioplot(c(0, align_enhancers_obs[,species]),at=c(1,4,7,10,13,16,19,22,25,28), col=rgb(t(col2rgb('forestgreen')/255), alpha = 0.6), border="forestgreen",
        add=T, axes=F, horizontal = T, plotCentre="line")

points(x=apply(align_enhancers_simul[,species], 2, mean), y=c(3,6,9,12,15,18,21,24,27), col = "white", pch=20, cex=0.8)

points(x = apply(align_enhancers_obs[,species], 2, mean), y = c(4,7,10,13,16,19,22,25,28), col = "white", pch=20, cex=0.8)

axis(1, pos=0.7, at=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.lab=1.2)
mtext("Alignment score", side=1, xpd = TRUE, cex=0.8)
mtext("C", side=3, line=1, at=0.1, font=2, cex=1.2)

######################## D - Conserved sequence Human to mouse vs distance to promoters ########################
par(mai = c(0.8, 0.6, 0.2, 0.2)) #bottom, left, top and right 
par(xpd=FALSE)
class_leg <- c("0", "0.5", "1", "1.5", "2")

xmin=0.2
xmax=0.4
CEX=1.2
CEX_lines=1

plot(obs_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="",
     xlab="", ylab="Alignment score", xaxt = "n", ylim=c(xmin,xmax), cex.lab=CEX, cex.axis=CEX, las=2)

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
mtext("Distance to promoters (Mb)", side=1, line=2, cex=0.8)
mtext("D", side=3, line=1, at=-1.5, font=2, cex=1.2)

######################## E - Conserved sequence to closest specie vs distance to promoters ######################## 
if(ref_sp=="human"){ymin=0.78; ymax=0.95}else{ymin=0.55; ymax=0.75}

plot(obs_dist_mac[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="",
     xlab="", ylab="Alignment score", xaxt = "n", ylim=c(ymin,ymax), cex.lab=CEX, cex.axis=CEX, las=2)

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
mtext("Distance to promoters (Mb)", side=1, line=2.5, cex=0.8)
mtext("E", side=3, line=1, at=-1.5, font=2, cex=1.2)


######################## F - Repeat proportion vs distance to promoters ######################## 
if(ref_sp=="human"){ymin=20; ymax=50}else{ymin=25; ymax=50}

plot(obs_repet_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="",
     xlab="", ylab="No-exonic repeat proportion (%)", xaxt = "n", ylim=c(ymin,ymax), cex.lab=CEX, cex.axis=CEX, las=2)

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
mtext("Distance to promoters (Mb)", side=1, line=2.5, cex=0.8)
mtext("F", side=3, line=1, at=-1.5, font=2, cex=1.2)

 ######################## G - Conserv enhancers vs distance to promoters ######################## 
par(mai = c(0.8, 0.6, 0.2, 1)) #bottom, left, top and right 
if(ref_sp=="human"){ymin=0.25; ymax=0.55}else{ymin=0.5; ymax=0.7}

color_n = 1 # To change color between each enhancers dataset
 
#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=F)
for (enh in enhancers){
   if (enh == "CAGE"){
     plot(list_conserv_enh[[enh]], type="l", col=col[color_n], main="",
          xlab="", ylab="Alignment score", xaxt = "n", ylim=c(ymin,ymax), las=2)
     
   }else{lines(list_conserv_enh[[enh]], type="l", col=col[color_n])}
   
   for (row in 1:length(list_conserv_enh[[enh]])){
     segments(x0=row,y0=list_conserv_enh[[paste0(enh, "_start")]][row],
              x1=row,y1=list_conserv_enh[[paste0(enh, "_end")]][row], col=col[color_n], lty=3, lwd=0.6)
   }
   
   color_n = color_n + 1 
 }
 
 
axis(1, at=seq(1,length(list_conserv_enh[[enh]])+1, 10), labels=F)
text(seq(1,length(list_conserv_enh[[enh]])+1,10), par("usr")[3]-0.04, class_leg, xpd = TRUE, cex=CEX)
mtext("Distance to promoters (Mb)", side=1, line=2.5, cex=0.8)
mtext("G", side=3, line=1, at=-1.5, font=2, cex=1.2)
 
par(xpd=TRUE)
enhancers_name = c("FANTOM5", "ENCODE")
if (ref_sp == "human"){enhancers_name <- c(enhancers_name, "RoadMap\nEpigenomics", "GRO-seq")}
 
legend("right", inset=c(-0.55,0), col=col, legend = enhancers, bty='n', lty=1)

dev.off()

