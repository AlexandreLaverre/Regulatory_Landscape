library(ape)
library(vioplot)
setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")

ref_sp = "mouse"
target_sp = "human"

path <- "/home/laverre/Data/Regulatory_landscape/result/Figures/"

if(ref_sp == "human"){pdf_name="Figure3_bis.pdf"}else{pdf_name="Sup_Figure11.pdf"}

pdf(paste(path, pdf_name, sep=""), width=7, height=10)
par(mai = c(0.5, 0.1, 0.3, 0.1)) #bottom, left, top and right 
a <- matrix(c(1,1,2,2,3,3),ncol=6, nrow = 2, byrow=T)
b <- matrix(c(4,4,4,5,5,5), nrow=1, byrow=F)
c <- matrix(c(6,6,6,7,7,7), nrow=1, byrow=F)
d <- rbind(a,b,c)
layout(d)

class_leg <- c("0", "0.5", "1", "1.5", "2")

xmin=0.25
xmax=0.45
CEX=1.2
CEX_lines=1


######################## F - Repeat proportion vs distance to promoters ######################## 
if(ref_sp=="human"){ymin=20; ymax=50}else{ymin=25; ymax=50}

plot(obs_repet_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="",
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
mtext("Distance to promoters (Mb)", side=1, line=2.5, cex=0.8)
mtext("F", side=3, line=1, at=-1.5, font=2, cex=1.2)

# ######################## G - Exonic proportion vs distance to promoters ######################## 
xmin=0
xmax=10

plot(obs_exon_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main="",
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
mtext("Distance to promoters (Mb)", side=1, line=2.5, cex=0.8)
mtext("G", side=3, line=1, at=-1.5, font=2, cex=1.2)
