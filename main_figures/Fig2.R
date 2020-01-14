setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")

load("Fig2_cell.Rdata")

CEX = 1.1
LWD = 1

png("Fig2.png", width = 800, height = 600)
par(mfrow=c(2,2))
### A - % overlap with enhancers

barcenter <- barplot(proportion, border=rep(c("black", "black", "white"),4), col=rep(c("firebrick1", "dodgerblue3", "white"),4), 
                     ylim=c(0,0.7), ylab="Enhancers presence proportion", axisnames = F, cex.lab=CEX, cex.axis=CEX)
text(c(1.4,4.7,8.5,12.1), par("usr")[3]-0.005, labels = enh_dataset, pos = 1, xpd = TRUE, cex=CEX)
legend("topright", legend = c("Observed", "Simulated"), fill=c("firebrick1", "dodgerblue3"), bty='n', cex=CEX)

segments(barcenter, conf_up, barcenter, conf_low, lwd = 3, cex=CEX)
arrows(barcenter, conf_up, barcenter, conf_low, lwd = 1.5, angle = 90, code = 3, length = 0.05, cex=CEX)

for (x in seq(1,length(barcenter)-1, by=3)){
  segments(barcenter[x], proportion[x]+0.04, barcenter[x+1], proportion[x]+0.04) 
  text("***", x=(barcenter[x]+barcenter[x+1])/2, y=proportion[x]+0.06, cex=CEX)
}

# B - % overlap as a function of number of cell types

plot(human_enh_cell$CAGE, type="l", col="red", lwd=LWD, ylab="Enhancers presence proportion", xlab="Number of cell", xaxt = "n",
     cex.lab=CEX, ylim=c(0,1), cex.axis=CEX)
for (row in 1:nrow(human_enh_cell)){
  segments(x0=row,y0=human_enh_cell[row,]$CAGE_conflow,x1=row,y1=human_enh_cell[row,]$CAGE_confup, col='red', lwd=LWD)}

points(human_enh_cell$ENCODE, type="l", col="dodgerblue3", lwd=LWD)
for (row in 1:nrow(human_enh_cell)){
  segments(x0=row,y0=human_enh_cell[row,]$ENCODE_conflow,x1=row,y1=human_enh_cell[row,]$ENCODE_confup, col="dodgerblue3", lwd=LWD)}

points(human_enh_cell$RoadMap, type="l", col="forestgreen", lwd=LWD)
for (row in 1:nrow(human_enh_cell)){
  segments(x0=row,y0=human_enh_cell[row,]$RoadMap_conflow,x1=row,y1=human_enh_cell[row,]$RoadMap_confup, col="forestgreen", lwd=LWD)}

points(human_enh_cell$GRO_seq, type="l", col="orange", lwd=LWD)
for (row in 1:nrow(human_enh_cell)){
  segments(x0=row,y0=human_enh_cell[row,]$GRO_seq_conflow,x1=row,y1=human_enh_cell[row,]$GRO_seq_confup, col='orange', lwd=LWD)}

class_leg <- c("1", "3", "5", "7", "9", "11", "13")
axis(1, at=seq(1,51,4), labels=F)
text(seq(1,51,4), par("usr")[3]-0.02, labels = class_leg, pos = 1, xpd = TRUE)
legend("topleft", legend=c("RoadMap", "ENCODE", "GRO_seq","CAGE"), fill=c("forestgreen", "dodgerblue3", "orange","red"), bty='n', cex=CEX-0.1)

# C - % overlap as a function of the genomic distance
plot(human_enh_dist$CAGE[0:50], type="l", col="red", lwd=LWD, ylab="Enhancers presence proportion",
     xlab="Linear distance to promoters regions (pb)", xaxt = "n", cex.lab=CEX, ylim=c(0,0.8), cex.axis=CEX)
for (row in 1:nrow(human_enh_dist[0:50,])){
  segments(x0=row,y0=human_enh_dist[row,]$CAGE_conflow,x1=row,y1=human_enh_dist[row,]$CAGE_confup, col='red', lwd=LWD)}

points(human_enh_dist_simul$CAGE[0:50], type="l", lty=2, col="red", lwd=0.6)
for (row in 1:nrow(human_enh_dist_simul[0:50,])){
 segments(x0=row,y0=human_enh_dist_simul[row,]$CAGE_conflow,x1=row,y1=human_enh_dist_simul[row,]$CAGE_confup, col="red", lwd=0.4)}

points(human_enh_dist$ENCODE[0:50], type="l", col="dodgerblue3", lwd=LWD)
for (row in 1:nrow(human_enh_dist[0:50,])){
  segments(x0=row,y0=human_enh_dist[row,]$ENCODE_conflow,x1=row,y1=human_enh_dist[row,]$ENCODE_confup, col="dodgerblue3", lwd=LWD)}

points(human_enh_dist_simul$ENCODE[0:50], type="l",  lty=2, col="dodgerblue3", lwd=0.6)
for (row in 1:nrow(human_enh_dist_simul[0:50,])){
 segments(x0=row,y0=human_enh_dist_simul[row,]$ENCODE_conflow,x1=row,y1=human_enh_dist_simul[row,]$ENCODE_confup, col="dodgerblue3", lwd=0.2)}

points(human_enh_dist$RoadMap[0:50], type="l", col="forestgreen", lwd=LWD)
for (row in 1:nrow(human_enh_dist[0:50,])){
  segments(x0=row,y0=human_enh_dist[row,]$RoadMap_conflow,x1=row,y1=human_enh_dist[row,]$RoadMap_confup, col="forestgreen", lwd=LWD)}

points(human_enh_dist_simul$RoadMap[0:50], type="l", lty=2,  col="forestgreen", lwd=0.6)
for (row in 1:nrow(human_enh_dist_simul[0:50,])){
 segments(x0=row,y0=human_enh_dist_simul[row,]$RoadMap_conflow,x1=row,y1=human_enh_dist_simul[row,]$RoadMap_confup, col="forestgreen", lwd=0.2)}

points(human_enh_dist$GRO_seq[0:50], type="l", col="orange", lwd=LWD)
for (row in 1:nrow(human_enh_dist[0:50,])){
  segments(x0=row,y0=human_enh_dist[row,]$GRO_seq_conflow,x1=row,y1=human_enh_dist[row,]$GRO_seq_confup, col='orange', lwd=LWD)}

points(human_enh_dist_simul$GRO_seq[0:50], type="l", lty=2, col="orange", lwd=0.6)
for (row in 1:nrow(human_enh_dist_simul[0:50,])){
 segments(x0=row,y0=human_enh_dist_simul[row,]$GRO_seq_conflow,x1=row,y1=human_enh_dist_simul[row,]$GRO_seq_confup, col="orange", lwd=0.2)}

class_leg <- c("25Kb", "500Kb", "1Mb", "1.5Mb", "2Mb", "2.5Mb", "3Mb", "3.5Mb", "4Mb")
axis(1, at=seq(1,81,10), labels=F)
text(seq(1,81,10), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE, cex=CEX)
legend("topright", legend=c("RoadMap", "ENCODE", "GRO_seq","FANTOM5"), fill=c("forestgreen", "dodgerblue3", "orange","red"), bty='n', cex=CEX-0.1)

dev.off()


########### Diff
plot(human_enh_dist$CAGE[0:50]-human_enh_dist_simul$CAGE[0:50], type="l", col="red", lwd=LWD, ylab="Enhancers presence proportion",
     xlab="Linear distance to promoters regions (pb)", xaxt = "n", cex.lab=CEX, ylim=c(0,0.25), cex.axis=CEX)
for (row in 1:nrow(human_enh_dist[0:50,])){
  segments(x0=row,y0=human_enh_dist[row,]$CAGE_conflow-human_enh_dist_simul[row,]$CAGE_conflow,x1=row,y1=human_enh_dist[row,]$CAGE_confup-human_enh_dist_simul[row,]$CAGE_confup, col='red', lwd=LWD)}

points(human_enh_dist$ENCODE[0:50]-human_enh_dist_simul$ENCODE[0:50], type="l", col="dodgerblue3", lwd=LWD)
for (row in 1:nrow(human_enh_dist[0:50,])){
  segments(x0=row,y0=human_enh_dist[row,]$ENCODE_conflow-human_enh_dist_simul[row,]$ENCODE_conflow,x1=row,y1=human_enh_dist[row,]$ENCODE_confup-human_enh_dist_simul[row,]$ENCODE_confup, col="dodgerblue3", lwd=LWD)}

points(human_enh_dist$RoadMap[0:50]-human_enh_dist_simul$RoadMap[0:50], type="l", col="forestgreen", lwd=LWD)
for (row in 1:nrow(human_enh_dist[0:50,])){
  segments(x0=row,y0=human_enh_dist[row,]$RoadMap_conflow-human_enh_dist_simul[row,]$RoadMap_conflow,x1=row,y1=human_enh_dist[row,]$RoadMap_confup-human_enh_dist_simul[row,]$RoadMap_confup, col="forestgreen", lwd=LWD)}

points(human_enh_dist$GRO_seq[0:50]-human_enh_dist_simul$GRO_seq[0:50], type="l", col="orange", lwd=LWD)
for (row in 1:nrow(human_enh_dist[0:50,])){
  segments(x0=row,y0=human_enh_dist[row,]$GRO_seq_conflow-human_enh_dist_simul[row,]$GRO_seq_conflow,x1=row,y1=human_enh_dist[row,]$GRO_seq_confup-human_enh_dist_simul[row,]$GRO_seq_confup, col='orange', lwd=LWD)}

class_leg <- c("25Kb", "500Kb", "1Mb", "1.5Mb", "2Mb", "2.5Mb", "3Mb", "3.5Mb", "4Mb")
axis(1, at=seq(1,81,10), labels=F)
text(seq(1,81,10), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE, cex=CEX)
legend("topright", legend=c("RoadMap", "ENCODE", "GRO_seq","FANTOM5"), fill=c("forestgreen", "dodgerblue3", "orange","red"), bty='n', cex=CEX-0.1)
