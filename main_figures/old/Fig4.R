setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")

load("Fig4_mouse_mean.Rdata")

CEX = 1
CEX_lines = 1
png("Fig4_mouse_mean.png", width = 800, height = 800)
par(mfrow=c(4,1), mai = c(0.5, 0.7, 0.5, 0.2))

layout(matrix(c(1,1,2,2,3,3,4,5), 4, 2, byrow = TRUE))

####### A - Large scale conservation #######
plot(obs_conserv$score, type="b", col="red", cex=CEX_lines, main=paste(sp_origin, " to ", sp_target, " conservation", sep=""), ylim=c(0,0.6),
     xlab="", ylab="Ungapped Non-exonic Score", cex.lab=CEX, cex.axis=CEX, cex.main=CEX, xaxt='n')
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=3)}

points(obs_conserv_enh$score[2:8], type="b", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_conserv_enh)){
  segments(x0=row,y0=obs_conserv_enh[row,]$int_start,x1=row,y1=obs_conserv_enh[row,]$int_end, col='forestgreen', lwd=0.3)}

points(simul_conserv$score[2:8], type="b", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_conserv)){
  segments(x0=row,y0=simul_conserv[row,]$int_start,x1=row,y1=simul_conserv[row,]$int_end, col='blue', lwd=0.3)}


legend("topright", fill=c("forestgreen","red","blue"), legend = c("Observed with enh", "Observed", "Simulated"), bty='n', cex=CEX)
axis(1, at=seq(1,9,1), labels=F)
text(seq(1,9,1), par("usr")[3]-0.04, row.names(obs_conserv), pos = 1, xpd = TRUE, cex=CEX)



boxplot(origin_conserv_seq_simul$human, origin_conserv_seq$human, origin_conserv_seq[which(origin_conserv_seq$CAGE_count > 1),]$human,
        notch=T, col=c("dodgerblue2", "firebrick1", "forestgreen"), ylim=c(0,1.1),
        names=c("Simulated", "Observed", "With enhancers"), ylab="Alignment score")

segments(1, 1.04, 2, 1.04) 
text("***", x=1.5, y=1.07, cex=1.5)
segments(2, 1.08, 3, 1.08) 
text("***",x= 2.5, y=1.1, cex=1.5)

####### B - Conservation ~ Distance to promoters ALL SP #######

par(mfrow=c(3,3), mai = c(0.3, 0.7, 0.3, 0.2))
layout( matrix(c(1,2,3,4,5,6,7,8,9),nrow = 3,ncol = 3,byrow = TRUE))

for (sp_target in rev(species)){
  if (sp_target == "opossum"){
    xmin=0
    xmax=0.08
  }else if (sp_target == "human"){
    xmin=0.2
    xmax=0.45
  }else if (sp_target == "rat"){
    xmin=0.5
    xmax=0.8
  }else if (sp_target == "rabbit"){
    xmin=0.15
    xmax=0.35
  }else if (sp_target == "dog"){
    xmin=0.2
    xmax=0.4
  }else if (sp_target == "elephant"){
    xmin=0.15
    xmax=0.35
  }else if (sp_target == "cow"){
    xmin=0.15
    xmax=0.4
  }else if (sp_target == "macaque"){
    xmin=0.2
    xmax=0.45
  }else if (sp_target == "chicken"){
    xmin=0
    xmax=0.04}
  
  plot(conserv_dist_all_sp[[sp_target]][["obs"]][1:50,"inter"], type="l", col="red", cex=CEX_lines, main=paste(sp_origin, " to ", sp_target, sep=""),
                                     xlab="", ylab="Ungapped Non-exonic Score", xaxt = "n", ylim=c(xmin,xmax), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
  for (row in 1:nrow(conserv_dist_all_sp[[sp_target]][["obs"]])){
    segments(x0=row,y0=conserv_dist_all_sp[[sp_target]][["obs"]][row,"int_start"],x1=row,y1=conserv_dist_all_sp[[sp_target]][["obs"]][row,"int_end"], col='red', lwd=0.3)}
  
  lines(conserv_dist_all_sp[[sp_target]][["simul"]][1:50,"inter"], type="l", col="blue", cex=CEX_lines)
  for (row in 1:nrow(conserv_dist_all_sp[[sp_target]][["simul"]])){
    segments(x0=row,y0=conserv_dist_all_sp[[sp_target]][["simul"]][row,"int_start"],x1=row,y1=conserv_dist_all_sp[[sp_target]][["simul"]][row,"int_end"], col='blue', lwd=0.3)}
  
  lines(conserv_dist_all_sp[[sp_target]][["obs_enh"]][1:50,"inter"], type="l", col="forestgreen", cex=CEX_lines)
  for (row in 1:nrow(conserv_dist_all_sp[[sp_target]][["obs_enh"]])){
    segments(x0=row,y0=conserv_dist_all_sp[[sp_target]][["obs_enh"]][row,"int_start"],x1=row,y1=conserv_dist_all_sp[[sp_target]][["obs_enh"]][row,"int_end"], col='forestgreen', lwd=0.3)}
  
  axis(1, at=seq(1,51,10), labels=F)
  text(seq(1,51,10),par("usr")[3]-0.02, class_leg, xpd = TRUE, cex=CEX)
  }

####### B - Conservation ~ Distance to promoters (mouse/mouse) #######

plot(obs_dist$inter[1:50], type="l", col="red", cex=CEX_lines, main=paste(sp_origin, " to ", sp_target, " conservation", sep=""),
           xlab="", ylab="Ungapped Non-exonic Score", xaxt = "n", ylim=c(0.15,0.4), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}

lines(simul_dist$inter[1:50], type="l", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}


lines(obs_dist_enh$inter[1:50], type="l", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist_enh[row,]$int_start,x1=row,y1=obs_dist_enh[row,]$int_end, col='forestgreen', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.02, class_leg, xpd = TRUE, cex=CEX)


#######  Conservation ~ Nb cell ~ Distance
plot(obs_dist_1$inter[1:50], type="l", col="red", cex=CEX_lines, main=paste("Conservation by presence in samples", sep=""),
     xlab="", ylab="Ungapped Non-exonic Score", xaxt = "n", ylim=c(0.15,0.45), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:nrow(obs_dist_1[1:50,])){
  segments(x0=row,y0=obs_dist_1[row,]$int_start,x1=row,y1=obs_dist_1[row,]$int_end, col='red', lwd=0.3)}

lines(obs_dist_2$inter[1:50], type="l", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_dist_2[1:50,])){
  segments(x0=row,y0=obs_dist_2[row,]$int_start,x1=row,y1=obs_dist_2[row,]$int_end, col='forestgreen', lwd=0.3)}

lines(obs_dist_3$inter[1:50], type="l", col="orange", cex=CEX_lines)
for (row in 1:nrow(obs_dist_3[1:50,])){
  segments(x0=row,y0=obs_dist_3[row,]$int_start,x1=row,y1=obs_dist_3[row,]$int_end, col='orange', lwd=0.3)}

lines(obs_dist_4$inter[1:50], type="l", col="blue", cex=CEX_lines)
for (row in 1:nrow(obs_dist_4[1:50,])){
  segments(x0=row,y0=obs_dist_4[row,]$int_start,x1=row,y1=obs_dist_4[row,]$int_end, col='blue', lwd=0.3)}

legend("topleft", fill=c("red","forestgreen","orange", "blue"), legend = c("1-3 samples", "3-10 samples", "10-20 samples", ">20 samples"), bty='n', cex=CEX)

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.02, class_leg, xpd = TRUE, cex=CEX)



#######  C - Overlap with PhastCons element #######
plot(obs_phastcons$inter[1:50], type="l", col="red", cex=CEX_lines,main=paste("PhastCons non-exonic composition", sep=""),
     xlab="", ylab="Phastcons proportion", xaxt = "n", ylim=c(0,0.04), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:nrow(obs_phastcons[1:50,])){
  segments(x0=row,y0=obs_phastcons[row,]$int_start,x1=row,y1=obs_phastcons[row,]$int_end, col='red', lwd=0.3)}

lines(simul_phastcons$inter[1:50], type="l", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_phastcons[1:50,])){
  segments(x0=row,y0=simul_phastcons[row,]$int_start,x1=row,y1=simul_phastcons[row,]$int_end, col='blue', lwd=0.3)}

lines(obs_phastcons_enh$inter[1:50], type="l", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_phastcons_enh[1:50,])){
  segments(x0=row,y0=obs_phastcons_enh[row,]$int_start,x1=row,y1=obs_phastcons_enh[row,]$int_end, col='forestgreen', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.004, class_leg, xpd = TRUE, cex=CEX)



#######  D - Overlap with Repeat element #######
CEX_lines = 1
CEX = 1
plot(obs_repeat$inter[1:50], type="l", col="red", cex=CEX_lines,main=paste("Repeat composition", sep=""),
     xlab="", ylab="Repeat proportion", xaxt = "n", ylim=c(0.35,0.5), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:nrow(obs_repeat[1:50,])){
  segments(x0=row,y0=obs_repeat[row,]$int_start,x1=row,y1=obs_repeat[row,]$int_end, col='red', lwd=0.3)}

lines(simul_repeat$inter[1:50], type="l", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_repeat[1:50,])){
  segments(x0=row,y0=simul_repeat[row,]$int_start,x1=row,y1=simul_repeat[row,]$int_end, col='blue', lwd=0.3)}

lines(obs_repeat_enh$inter[1:50], type="l", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_repeat_enh[1:50,])){
  segments(x0=row,y0=obs_repeat_enh[row,]$int_start,x1=row,y1=obs_repeat_enh[row,]$int_end, col='forestgreen', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.02, class_leg, xpd = TRUE, cex=CEX)



#######  E - Conservation between mouse and mouse ~ Number of cell  #######
plot(obs_nb_cell$stats[3,], type="b", col="red", main=paste(sp_origin, " contacted sequences conservation",sep=""), cex=CEX_lines, ylim=c(0.2,0.4),
    xlab="", ylab="Ungapped Non-exonic Score", xaxt='n', cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:ncol(obs_nb_cell$stats)){segments(x0=row,y0=obs_nb_cell$conf[1,row],x1=row,y1=obs_nb_cell$conf[2,row], col='red')}

points(obs_nb_cell_enh$stats[3,], col="forestgreen", type="b", cex=CEX_lines)
for (row in 1:ncol(obs_nb_cell_enh$stats)){segments(x0=row,y0=obs_nb_cell_enh$conf[1,row],x1=row,y1=obs_nb_cell_enh$conf[2,row], col='forestgreen')}

axis(1, at=seq(1,7,1), labels=F)
text(seq(1,7,1), par("usr")[3]-0.01, labels = obs_nb_cell$names, pos = 1, xpd = TRUE, cex=CEX)



dev.off()


if (sp_target == "opossum"){
  xmin=0.04
  xmax=0.13
}else if (sp_target == "mouse"){
  xmin=0.2
  xmax=0.45
}else if (sp_target == "rat"){
  xmin=0.25
  xmax=0.45
}else if (sp_target == "rabbit"){
  xmin=0.3
  xmax=0.55
}else if (sp_target == "dog"){
  xmin=0.4
  xmax=0.6
}else if (sp_target == "elephant"){
  xmin=0.3
  xmax=0.5
}else if (sp_target == "cow"){
  xmin=0.35
  xmax=0.55
}else if (sp_target == "macaque"){
  xmin=0.8
  xmax=0.95
}else if (sp_target == "chicken"){
  xmin=0
  xmax=0.04}

