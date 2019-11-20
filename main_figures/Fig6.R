setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")

sp_origin = "human"
species <- c("macaque",  "dog", "cow", "elephant", "rabbit", "rat", "mouse", "opossum", "chicken")

load(paste("Fig6_", sp_origin,".Rdata", sep=""))

CEX = 2
CEX_lines = 3
png("Fig6_human2.png", width = 800, height = 800)
par(mfrow=c(3,1), mai = c(0.5, 0.7, 0.5, 0.2))

layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE))


# Conserv syntenie cis
null = rep(c(0,0,0,NA),9) 
plot(null, col = NA, ylim=c(50,100), xaxt='n', main="Human conserved syntenie (cis)",
     ylab='Conserved synteny (%)', xlab="",  cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
points(data, col=rep(c("blue", "red", "forestgreen", "black"),9), pch=3, cex=CEX_lines)
axis(1, at=seq(2,36,4), labels=F)
text(seq(2,36,4), par("usr")[3]-2, labels = species, pos = 1, xpd = TRUE, cex=CEX)
legend("bottomleft", fill=c("blue","red", "forestgreen"), legend = c("Simulated", "Oberved", "With enhancers"), bty='n', cex=CEX)



# Conserv syntenie < 2Mb
null = rep(c(0,0,0,NA),9) 

plot(null, col = NA, ylim=c(50,100), xaxt='n', ylab='Conserved synteny (%)', main="Human conserved syntenie (cis and <2MB)",
     xlab="", cex.lab=1.2, cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
points(data_2M, col=rep(c("blue", "red", "forestgreen", "black"),9), pch=3, cex=CEX_lines)
axis(1, at=seq(2,36,4), labels=F)
text(seq(2,36,4), par("usr")[3]-2, labels = species, pos = 1, xpd = TRUE,cex=CEX)

### C - Syntenie conservation (cis) ~~ genomic distance ###
plot(obs_dist$inter[1:50], type="l", col="red", cex=CEX_lines, main=paste(sp_origin, " to ", sp_target, " conserved syntenie (cis)", sep=""),
     xlab="", ylab="Ungapped Non-exonic Score", xaxt = "n", ylim=c(70,100), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}

lines(simul_dist$inter[1:50], type="l", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}


lines(obs_dist_enh$inter[1:50], type="l", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist_enh[row,]$int_start,x1=row,y1=obs_dist_enh[row,]$int_end, col='forestgreen', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-2, class_leg, xpd = TRUE, cex=CEX)

### D - Syntenie conservation (cis+2MB) ~~ genomic distance ###
plot(obs_dist_2M$inter, type="l", col="red", cex=CEX_lines, main=paste(sp_origin, " to ", sp_target, " conserved syntenie (cis and <2Mb)", sep=""),
     xlab="", ylab="", xaxt = "n", ylim=c(0,100), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
for (row in 1:nrow(obs_dist_2M)){
  segments(x0=row,y0=obs_dist_2M[row,]$int_start,x1=row,y1=obs_dist_2M[row,]$int_end, col='red', lwd=0.3)}

lines(simul_dist_2M$inter, type="l", col="blue", cex=CEX_lines)
for (row in 1:nrow(simul_dist_2M)){
  segments(x0=row,y0=simul_dist_2M[row,]$int_start,x1=row,y1=simul_dist_2M[row,]$int_end, col='blue', lwd=0.3)}

lines(obs_dist_enh_2M$inter, type="l", col="forestgreen", cex=CEX_lines)
for (row in 1:nrow(obs_dist)){
  segments(x0=row,y0=obs_dist_enh_2M[row,]$int_start,x1=row,y1=obs_dist_enh_2M[row,]$int_end, col='forestgreen', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-8, class_leg, xpd = TRUE, cex=CEX)

dev.off()

# 
# segments(1, 97, 2, 97) 
# text("***", x= 1.5, y=98, cex=1.3)
# segments(2, 99, 3, 99) 
# text("*", x= 2.5, y=100, cex=1.3)
# 
# segments(5, 99, 6, 99) 
# text("***",x=5.5, y=100, cex=1.3)
# segments(6, 101, 7, 101) 
# text("**", x=6.5, y=102, cex=1.3)
# 
# segments(9, 98, 10, 98) 
# text("***",x= 9.5, y=99, cex=1.3)
# segments(10, 100, 11, 100) 
# text("*", x= 10.5, y=101, cex=1.3)
# 
# segments(13, 96, 14, 96) 
# text("***",x= 13.5, y=97, cex=1.3)
# segments(14, 98, 15, 98) 
# text("NS", x=14.5, y=99, cex=1)
# 
# segments(17, 84, 18, 84)
# text("***",x= 17.5, y=85, cex=1.3)
# segments(18, 86, 19, 86) 
# text("NS", x= 18.5, y=87, cex=1)
# 
# segments(21, 73, 22, 73)
# text("***",x= 21.5, y=74, cex=1.4)
# segments(22, 75, 23, 75) 
# text("***", x= 22.5, y=76, cex=1.4)
