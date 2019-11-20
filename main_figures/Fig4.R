setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")

load("Fig4_mouse_mean.Rdata")


### Filtres ###
# Length
origin_conserv_seq <- origin_conserv_seq[which(origin_conserv_seq$length > 250 & origin_conserv_seq$length < 20000),]
origin_conserv_seq_simul$length <- origin_conserv_seq_simul$end - origin_conserv_seq_simul$start
origin_conserv_seq_simul <- origin_conserv_seq_simul[which(origin_conserv_seq_simul$length > 250 & origin_conserv_seq_simul$length < 20000),]

# Duplication
origin_conserv_seq <- origin_conserv_seq[which(origin_conserv_seq$duplication == 0),]
origin_conserv_seq_simul <- origin_conserv_seq_simul[which(origin_conserv_seq_simul$duplication == 0),]
# Baited
origin_conserv_seq <- origin_conserv_seq[which(origin_conserv_seq$baited == "unbaited"),]
origin_conserv_seq_simul <- origin_conserv_seq_simul[which(origin_conserv_seq_simul$baited == "unbaited"),]

CEX = 1.8
CEX_lines = 1.3
png("Fig4_mouse_mean.png", width = 800, height = 800)
par(mfrow=c(4,1), mai = c(0.5, 0.7, 0.5, 0.2))

layout(matrix(c(1,1,2,2,3,3,4,5), 4, 2, byrow = TRUE))

####### A - Large scale conservation #######
plot(obs_conserv$score, type="b", col="red", cex=CEX_lines, main=paste(sp_origin, " to ", sp_target, " conservation", sep=""), ylim=c(0,1),
     xlab="", ylab="Ungapped Non-exonic Score", cex.lab=CEX, cex.axis=CEX, cex.main=CEX, xaxt='n')
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}

points(obs_conserv_enh$score, type="b", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_conserv_enh)){
  segments(x0=row,y0=obs_conserv_enh[row,]$int_start,x1=row,y1=obs_conserv_enh[row,]$int_end, col='forestgreen', lwd=0.3)}

points(simul_conserv$score, type="b", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_conserv)){
  segments(x0=row,y0=simul_conserv[row,]$int_start,x1=row,y1=simul_conserv[row,]$int_end, col='blue', lwd=0.3)}


legend("topright", fill=c("forestgreen","red","blue"), legend = c("Observed with enh", "Observed", "Simulated"), bty='n', cex=CEX)
axis(1, at=seq(1,9,1), labels=F)
text(seq(1,9,1), par("usr")[3]-0.04, labels = row.names(obs_conserv), pos = 1, xpd = TRUE, cex=CEX)



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
     xlab="", ylab="Phastcons proportion", xaxt = "n", ylim=c(0.0,0.07), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:nrow(obs_phastcons[1:50,])){
  segments(x0=row,y0=obs_phastcons[row,]$int_start,x1=row,y1=obs_phastcons[row,]$int_end, col='red', lwd=0.3)}

lines(simul_phastcons$inter[1:50], type="l", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_phastcons[1:50,])){
  segments(x0=row,y0=simul_phastcons[row,]$int_start,x1=row,y1=simul_phastcons[row,]$int_end, col='blue', lwd=0.3)}

lines(obs_phastcons_enh$inter[1:50], type="l", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_phastcons_enh[1:50,])){
  segments(x0=row,y0=obs_phastcons_enh[row,]$int_start,x1=row,y1=obs_phastcons_enh[row,]$int_end, col='forestgreen', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.008, class_leg, xpd = TRUE, cex=CEX)



#######  D - Overlap with Repeat element #######
plot(obs_repeat$inter[1:50], type="l", col="red", cex=CEX_lines,main=paste("Repeat composition", sep=""),
     xlab="", ylab="Repeat proportion", xaxt = "n", ylim=c(0.2,0.52), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:nrow(obs_repeat[1:50,])){
  segments(x0=row,y0=obs_repeat[row,]$int_start,x1=row,y1=obs_repeat[row,]$int_end, col='red', lwd=0.3)}

lines(simul_repeat$inter[1:50], type="l", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_repeat[1:50,])){
  segments(x0=row,y0=simul_repeat[row,]$int_start,x1=row,y1=simul_repeat[row,]$int_end, col='blue', lwd=0.3)}

lines(obs_repeat_enh$inter[1:50], type="l", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_repeat_enh[1:50,])){
  segments(x0=row,y0=obs_repeat_enh[row,]$int_start,x1=row,y1=obs_repeat_enh[row,]$int_end, col='forestgreen', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.04, class_leg, xpd = TRUE, cex=CEX)



#######  E - Conservation between mouse and mouse ~ Number of cell  #######
# plot(obs_nb_cell$stats[3,], type="b", col="red", main=paste(sp_origin, " contacted sequences conservation",sep=""), cex=CEX_lines, ylim=c(0.2,0.4),
#      xlab="", ylab="Ungapped Non-exonic Score", xaxt='n', cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
# for (row in 1:ncol(obs_nb_cell$stats)){segments(x0=row,y0=obs_nb_cell$conf[1,row],x1=row,y1=obs_nb_cell$conf[2,row], col='red')}
# 
# points(obs_nb_cell_enh$stats[3,], col="forestgreen", type="b", cex=CEX_lines)
# for (row in 1:ncol(obs_nb_cell_enh$stats)){segments(x0=row,y0=obs_nb_cell_enh$conf[1,row],x1=row,y1=obs_nb_cell_enh$conf[2,row], col='forestgreen')}
# 
# axis(1, at=seq(1,7,1), labels=F)
# text(seq(1,7,1), par("usr")[3]-0.01, labels = obs_nb_cell$names, pos = 1, xpd = TRUE, cex=CEX)



dev.off()

