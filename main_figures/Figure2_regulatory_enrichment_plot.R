############################################################ PLOT FIGURE 2 ###############################################
path <- "/home/laverre/Data/Regulatory_landscape/result/Main_figures/"
load(paste(path, "Fig2_human.Rdata", sep=""))

ref_sp = "human"
enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")


pdf(paste(path, "Figure2_human.pdf", sep=""), width=8.5, height=8)
par(mai = c(0.8, 0.8, 0.5, 0.2)) # bottom, left, top, right
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))

############################################   A - Global enhancer proportion ############################################ 
barcenter <- barplot(enh_prop$data, border=rep(c("darkgreen", "firebrick3", "white"),4), col="white", lwd=1.5, cex.names=0.8,
                     ylim=c(0,0.15), ylab="Enhancer length proportion (mean)", axisnames = F, main="Fig2.A.", las=2)


text(c(1.4,4.7,8.5,12.1), par("usr")[3]-0.005, labels = enhancers, pos = 1, xpd = TRUE, cex=0.8)
legend("topleft", legend = c("Observed", "Simulated"), border=c("darkgreen", "firebrick3"), fill="white", bty='n')
par(lwd=1)

segments(barcenter, enh_prop$conf_up, barcenter, enh_prop$conf_low, lwd = 3)
arrows(barcenter, enh_prop$conf_up, barcenter, enh_prop$conf_low, lwd = 1.5, angle = 90, code = 3, length = 0.05)

for (x in seq(1,length(barcenter)-1, by=3)){
  segments(barcenter[x], enh_prop$data[x]+0.01, barcenter[x+1], enh_prop$data[x]+0.01) 
  text("***", x=(barcenter[x]+barcenter[x+1])/2, y=enh_prop$data[x]+0.015)
}

############################################   B - Global enhancer proportion ############################################ 
color <- c("red", "navy", "forestgreen", "orange")

plot(prop_nb_sample[["obs"]]$CAGE, col="white", ylim=c(0,0.2), las=2,
     ylab="Enhancer length proportion (mean)", xlab="Number of sample", xaxt = "n", main="Fig2.B.")

col_nb = 1
for (enh in enhancers){
  points(prop_nb_sample[["obs"]][[enh]], type="l", col=color[col_nb])
  for (row in 1:nrow(prop_nb_sample[["obs"]])){
    segments(x0=row,y0=prop_nb_sample[["obs"]][row,paste0(enh, "_conflow")],
             x1=row,y1=prop_nb_sample[["obs"]][row,paste0(enh, "_confup")], col=color[col_nb], lwd=0.5)}
  
  points(prop_nb_sample[["simul"]][[enh]], type="l", lty=3, col=color[col_nb])
  for (row in 1:nrow(prop_nb_sample[["simul"]])){
    segments(x0=row,y0=prop_nb_sample[["simul"]][row,paste0(enh, "_conflow")],
             x1=row,y1=prop_nb_sample[["simul"]][row,paste0(enh, "_confup")], col=color[col_nb], lwd=0.3)}
  
  col_nb = col_nb + 1
}

if (ref_sp == "mouse"){class_leg <- c("1", "3", "5", "7", "9", "11", "13"); max_nb_cell=13
}else{class_leg <- c("1", "5", "10", "15", "20", "25"); max_nb_cell=25}

axis(1, at=seq(1,max_nb_cell+1,5), labels=F)
text(seq(1,max_nb_cell+1,5), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE)
legend("topleft", legend=enhancers, col=color, bty='n', lty=1, cex=0.8)

############################################   C - Global enhancer proportion ############################################ 
plot(prop_dist[["obs"]]$CAGE[0:50], type="l", col="white", ylab="Enhancer length proportion (mean)", main="Fig2.C.", las=2,
     xlab="Linear distance to promoters regions (pb)", xaxt = "n", ylim=c(0,0.2))

nb_col = 1
for (enh in enhancers){
  
  points(prop_dist[["obs"]][[paste0(enh)]][0:50], type="l", col=color[nb_col])
  for (row in 1:nrow(prop_dist[["obs"]][0:50,])){
    segments(x0=row,y0=prop_dist[["obs"]][row,paste0(enh, "_conflow")],
             x1=row,y1=prop_dist[["obs"]][row,paste0(enh, "_confup")], col=color[nb_col], lwd=0.5)}
  
  points(prop_dist[["simul"]][[paste0(enh)]][0:50], type="l", lty=2, col=color[nb_col], lwd=0.6)
  for (row in 1:nrow(prop_dist[["simul"]][0:50,])){
    segments(x0=row,y0=prop_dist[["simul"]][row,paste0(enh, "_conflow")],
             x1=row,y1=prop_dist[["simul"]][row,paste0(enh, "_confup")], col=color[nb_col], lwd=0.5)}
  
  nb_col = nb_col + 1
  
}

class_leg <- c("25Kb", "500Kb", "1Mb", "1.5Mb", "2Mb", "2.5Mb", "3Mb", "3.5Mb", "4Mb")
axis(1, at=seq(1,81,10), labels=F)
text(seq(1,81,10), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE)
legend("topright", legend=enhancers, col=color, bty='n', lty=1, cex=0.8)

dev.off()