library(ggplot2)
library(gridExtra)
setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

sp_origin = 'human'
sp_target = 'mouse'

obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv.txt", sep=""), header=T)
simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_simul_with_notconserv.txt", sep=""), header=T)

obs$class <-cut(obs$origin_dist, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
simul$class <- cut(simul$origin_dist, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
class_leg <- c("25-75Kb", "1.97-2.02Mb", "3.97-4.02Mb","5.97-6.02Mb","7.97-8.02Mb","9.97-10.02Mb")

# Syntenie ~ distance sp origin
obs_conserv <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),])/ nrow(obs[which(obs$class == x & !is.na(obs$PIR_lift)),]))*100))
obs_conserv$int_start <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x & !is.na(obs$PIR_lift)),])+1, p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x & !is.na(obs$PIR_lift)),])+1, p=0.5)$conf.int[2])*100)

(nrow(obs[which(!is.na(obs$target_dist)),])/nrow(obs[which(!is.na(obs$PIR_lift)),]))*100
(nrow(simul[which(!is.na(simul$target_dist)),])/nrow(simul[which(!is.na(simul$PIR_lift)),]))*100

simul_conserv <- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),])/ nrow(simul[which(simul$class == x & !is.na(simul$PIR_lift)),]))*100))
simul_conserv$int_start <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),]), n=nrow(simul[which(simul$class == x & !is.na(simul$PIR_lift)),]), p=0.5)$conf.int[1])*100)
simul_conserv$int_end <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),]), n=nrow(simul[which(simul$class == x & !is.na(simul$PIR_lift)),]), p=0.5)$conf.int[2])*100)

plot(obs_conserv$inter, type="b", col="red", cex=0.4, main=paste("Syntenie", sp_origin,"to", sp_target),
     ylab="Proportion interaction intra (%)", xlab="", xaxt = "n")
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col=rgb(1,0,0,0.5), lwd=0.3)}
points(simul_conserv$inter, type="b", col="blue", cex=0.4)
for (row in 1:nrow(simul_conserv)){
  segments(x0=row,y0=simul_conserv[row,]$int_start,x1=row,y1=simul_conserv[row,]$int_end, col=rgb(0,0,1,0.5), lwd=0.3)}

axis(1, at=seq(1,201,40), labels=F)
text(seq(1,201,40), par("usr")[3]-0.6, labels = class_leg, pos = 1, xpd = TRUE, cex=0.8)
legend("bottomleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

# Proportion interaction > 10Mb sp target ~ distance sp origin
obs_conserv <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & obs$target_dist > 10000000),])/ nrow(obs[which(obs$class == x & !is.na(obs$PIR_lift)),]))*100))
obs_conserv$int_start <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & obs$target_dist > 10000000),]), n=nrow(obs[which(obs$class == x & !is.na(obs$PIR_lift)),])+1, p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & obs$target_dist > 10000000),]), n=nrow(obs[which(obs$class == x & !is.na(obs$PIR_lift)),])+1, p=0.5)$conf.int[2])*100)

(nrow(obs[which(obs$target_dist > 2000000),])/nrow(obs[which(!is.na(obs$PIR_lift)),]))*100
(nrow(simul[which(simul$target_dist > 2000000),])/nrow(simul[which(!is.na(simul$PIR_lift)),]))*100

simul_conserv <- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & simul$target_dist > 10000000),])/ nrow(simul[which(simul$class == x & !is.na(simul$PIR_lift)),]))*100))
simul_conserv$int_start <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & simul$target_dist > 10000000),]), n=nrow(simul[which(simul$class == x & !is.na(simul$PIR_lift)),]), p=0.5)$conf.int[1])*100)
simul_conserv$int_end <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & simul$target_dist > 10000000),]), n=nrow(simul[which(simul$class == x & !is.na(simul$PIR_lift)),]), p=0.5)$conf.int[2])*100)

plot(obs_conserv$inter, type="b", col="red", cex=0.4,
     main=paste(sp_origin,"interactions conserved at more than 10Mb in", sp_target),
     ylab="Proportion interaction >10Mb (%)", xlab="", xaxt = "n")
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col=rgb(1,0,0,0.5), lwd=0.3)}
points(simul_conserv$inter, type="b", col="blue", cex=0.4)
for (row in 1:nrow(simul_conserv)){
  segments(x0=row,y0=simul_conserv[row,]$int_start,x1=row,y1=simul_conserv[row,]$int_end, col=rgb(0,0,1,0.5), lwd=0.3)}

axis(1, at=seq(1,201,40), labels=F)
text(seq(1,201,40), par("usr")[3]-0.6, labels = class_leg, pos = 1, xpd = TRUE, cex=0.8)
legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

############## sp_origin distance ~ sp_target distance ##############
obs <- obs[which(!is.na(obs$target_dist)),]
simul <- simul[which(!is.na(simul$target_dist)),]

obs_int_conserv <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction.txt2", sep=""), header=T)
obs_int_conserv <- obs_int_conserv[which(!is.na(obs_int_conserv$target_dist)),]
obs$int_conserv <- obs$origin_interaction %in% obs_int_conserv$origin_interaction

simul_int_conserv <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction_simul.txt2", sep=""), header=T)
simul_int_conserv <- simul_int_conserv[which(!is.na(simul_int_conserv$target_dist)),]
simul$int_conserv <- simul$origin_interaction %in% simul_int_conserv$origin_interaction

simul$color <- 0
simul[which(simul$int_conserv == TRUE),]$color <- "red"
simul[which(simul$int_conserv == FALSE),]$color <- "grey50"
obs$color <- 0
obs[which(obs$int_conserv == TRUE),]$color <- "red"
obs[which(obs$int_conserv == FALSE),]$color <- "grey50"

par(mfrow=c(1,1))
## Observed plot
plot(log2(obs$origin_dist+1)~log2(obs$target_dist+1), col=obs$color, pch=19,
     cex=0.01, xlab=paste("log2(", sp_target,"dist) (pb)"), ylab=paste("log2(", sp_origin,"dist) (pb)"),  
     main=paste("Distance comparison of", sp_origin, "observed interactions with", sp_target), xlim=c(12,25))
abline(a=0,b=1)
abline(lm(log2(obs$origin_dist+1)~log2(obs$target_dist+1)), col='blue')
abline(lm(log2(obs[which(obs$int_conserv == TRUE),]$origin_dist+1)~log2(obs[which(obs$int_conserv == TRUE),]$target_dist+1)), col='green')
summary(lm(log2(obs$origin_dist+1)~log2(obs$target_dist+1)))

legend(x ="topleft",legend = c("Conserved interactions", "Other" ), col = c('red','grey50'), pch = 19, bty='n')
legend(x = 22.5, y=22,legend = "y = 0.81 x + 3.47", bty = 'n')

## Simulated plot
plot(log2(simul$origin_dist+1)~log2(simul$target_dist+1), col=simul$color, pch=19,
     cex=0.01, xlab=paste("log2(", sp_target,"dist) (pb)"), ylab=paste("log2(", sp_origin,"dist) (pb)"),  
     main=paste("Distance comparison of", sp_origin, "simulated interactions with", sp_target), xlim=c(12,25))
abline(a=0,b=1)
abline(lm(log2(simul$origin_dist+1)~log2(simul$target_dist+1)), col='blue')
abline(lm(log2(simul[which(simul$int_conserv == TRUE),]$origin_dist+1)~log2(simul[which(simul$int_conserv == TRUE),]$target_dist+1)), col='green')
summary(lm(log2(simul$origin_dist+1)~log2(simul$target_dist+1)))

legend(x ="topleft",legend = c("Conserved interactions", "Other" ), col = c('red','grey50'), pch = 19, bty='n')
legend(x = 22.5, y=22,legend = "y = 0.81 x + 3.47", bty = 'n')

############## Ratio log2 ~ PIR score ################
obs$ratio <- log2((obs$target_dist+1)/(obs$origin_dist+1))
simul$ratio <- log2((simul$target_dist+1)/(simul$origin_dist+1))
obs$conserv_class <-cut(obs$PIR_score, breaks=seq(from=0, to=1, by=0.1), include.lowest = T)
simul$conserv_class <-cut(simul$PIR_score, breaks=seq(from=0, to=1, by=0.1), include.lowest = T)

# boxplot(obs$ratio~obs$conserv_class, xlim = c(0.5, length(levels(obs$conserv_class))+0.5), outline=F,
#         boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), ylim=c(-1.05,0.5),
#         main=paste(sp_origin, 'to', sp_target), xlab='Score conservation',
#         ylab=paste("log2(dist", sp_target, "/ dist", sp_origin, ")")) #invisible boxes
# 
# boxplot(obs[which(obs$int_conserv == TRUE),]$ratio~obs[which(obs$int_conserv == TRUE),]$conserv_class,
#         xaxt = "n", yaxt='n', add = TRUE, boxfill="red", boxwex=0.2, outline=F, notch=T,
#         at = 1:length(levels(obs$conserv_class)) - 0.25) #shift these left by -0.15
# 
# boxplot(obs[which(obs$int_conserv == FALSE),]$ratio~obs[which(obs$int_conserv == FALSE),]$conserv_class,
#         xaxt = "n", yaxt='n', add = TRUE, boxfill="grey", boxwex=0.2, outline=F, notch=T,
#         at = 1:length(levels(obs$conserv_class))) #shift these left by -0.15
# 
# boxplot(simul$ratio~simul$conserv_class,
#         xaxt = "n", yaxt='n', add = TRUE, boxfill="blue", boxwex=0.2, outline=F, notch=T,
#         at = 1:length(levels(obs$conserv_class)) + 0.25) #shift these left by -0.15
# 
# legend("bottomright", legend = c("Observed and conserved", "Observed not conserved", "Simulated" ), 
#        col = c('red','grey', 'blue'), pch = 15, bty='n')

obs_box_cons <- boxplot(obs[which(obs$int_conserv == TRUE),]$ratio~obs[which(obs$int_conserv == TRUE),]$conserv_class, plot=F)
obs_box <- boxplot(obs$ratio~obs$conserv_class, plot=F)
simul_box_cons <- boxplot(simul[which(simul$int_conserv == TRUE),]$ratio~simul[which(simul$int_conserv == TRUE),]$conserv_class, plot=F)
simul_box <- boxplot(simul$ratio~simul$conserv_class, plot=F)

#ylim=c(0.1,0.35)
plot(obs_box_cons$stats[3,], type="b", col="red", main=paste(sp_origin, 'to', sp_target), ylim=c(-0.3,-0.1),
     xlab="Score conservation", ylab=paste("log2(dist", sp_target, "/ dist", sp_origin, ")"), xaxt='n')
for (row in 1:ncol(obs_box_cons$stats)){segments(x0=row,y0=obs_box_cons$conf[1,row],x1=row,y1=obs_box_cons$conf[2,row], lwd=0.3, col='red')}

points(obs_box$stats[3,], type="b")
for (row in 1:ncol(obs_box$stats)){segments(x0=row,y0=obs_box$conf[1,row],x1=row,y1=obs_box$conf[2,row], lwd=0.3)}

points(simul_box$stats[3,], col="blue", type="b")
for (row in 1:ncol(simul_box$stats)){segments(x0=row,y0=simul_box$conf[1,row],x1=row,y1=simul_box$conf[2,row], lwd=0.3, col='blue')}

points(simul_box_cons$stats[3,], col="green", type="b")
for (row in 1:ncol(simul_box_cons$stats)){segments(x0=row,y0=simul_box_cons$conf[1,row],x1=row,y1=simul_box_cons$conf[2,row], lwd=0.3, col='green')}

legend("topleft", fill=c("blue","black","red", "green"), legend = c("Simulated", "Observed", "Observed and conserved", "Simulated and conserved"), bty='n')
axis(1, at=seq(2,10,2), labels=F)
text(seq(2,10,2), par("usr")[3]-0.005, labels = c("0.1 - 0.2","0.3 - 0.4","0.5 - 0.6","0.7 - 0.8",'0.9 - 1'), pos = 1, xpd = TRUE, cex=0.8)


############## Ratio log2 ~ distance ############## 
intergenic <- read.table(paste("intergenic_distance_",sp_origin,"2",sp_target,".txt",sep=""), header=T)
colnames(intergenic) <- c("gene", "origin_dist", "target_dist")
intergenic$ratio <- log2(intergenic$target_dist/intergenic$origin_dist)
intergenic2 <- intergenic[which(intergenic$origin_dist < 3000000 ),]
obs2 <- obs[which(!is.na(obs$target_dist) & obs$origin_dist < 3000000 ),]
simul2 <- simul[which(!is.na(simul$target_dist) & simul$origin_dist < 3000000),]
obs2$dist_class <-cut(obs2$origin_dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
simul2$dist_class <-cut(simul2$origin_dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
intergenic2$dist_class <-cut(intergenic2$origin_dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
class_leg <- c("25-75Kb", "975Kb-1.02Mb", "1.97-2.02Mb","2.97-3.02Mb", "3.97-4.02Mb", "4.97-5.02Mb")


par(mfrow=c(1,1))
obs2_box_cons <- boxplot(obs2[which(obs2$int_conserv == TRUE),]$ratio~obs2[which(obs2$int_conserv == TRUE),]$dist_class, plot=F)
obs2_box <- boxplot(obs2$ratio~obs2$dist_class, plot=F)
simul2_box <- boxplot(simul2$ratio~simul2$dist_class, plot=F)
intergenic2_box <- boxplot(intergenic2$ratio~intergenic2$dist_class, plot=F)

# ylim=c(-0.3,0)
# ylim=c(-0.1,0.35)
plot(obs2_box_cons$stats[3,], type="b", col="red", main=paste(sp_origin, 'to', sp_target), ylim=c(-0.3,0), cex=0.7,
     xlab="", ylab=paste("log2(dist", sp_target, "/ dist", sp_origin, ")"), xaxt='n')
for (row in 1:ncol(obs2_box_cons$stats)){segments(x0=row,y0=obs2_box_cons$conf[1,row],x1=row,y1=obs2_box_cons$conf[2,row], lwd=0.3, col='red')}

points(obs2_box$stats[3,], type="b", cex=0.7)
for (row in 1:ncol(obs2_box$stats)){segments(x0=row,y0=obs2_box$conf[1,row],x1=row,y1=obs2_box$conf[2,row], lwd=0.3)}

points(simul2_box$stats[3,], col="blue", type="b", cex=0.7)
for (row in 1:ncol(simul2_box$stats)){segments(x0=row,y0=simul2_box$conf[1,row],x1=row,y1=simul2_box$conf[2,row], lwd=0.3, col='blue')}

points(intergenic2_box$stats[3,], col="green", type="b", cex=0.7)
for (row in 1:ncol(intergenic2_box$stats)){segments(x0=row,y0=intergenic2_box$conf[1,row],x1=row,y1=intergenic2_box$conf[2,row], lwd=0.3, col='blue')}


legend("topleft", fill=c("blue","black","red", "green"), legend = c("Simulated", "Observed", "Observed and conserved", "Intergenic"), bty='n')
axis(1, at=seq(1,101,20), labels=F)
text(seq(1,101,20), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE, cex=0.8)


# plot(obs2$ratio~obs2$origin_dist, ylim=c(-4, 4), cex=0.1, pch=19, main=paste(sp_origin, "to", sp_target), xlab=paste(sp_origin, "distance"), ylab=paste('log2(dist', sp_target, '/dist', sp_origin,')'))
# abline(h=0)
# abline(lm(obs2$ratio~obs2$origin_dist), col="red")
# summary(lm(obs2$ratio~obs2$origin_dist))
# med = signif(median(obs2$ratio), digits=3)
# legend("topright", legend=c(paste0("Med = ", med),"y = -0.29 + 1.8e-8 x"), bty="n")
# 
# plot(simul2$ratio~simul2$origin_dist, ylim=c(-4, 4), cex=0.1, pch=19, main=paste(sp_origin, "to", sp_target, 'simulation'), xlab=paste(sp_origin, "distance"), ylab='')
# abline(h=0)
# abline(lm(simul2$ratio~simul2$origin_dist), col="red")
# summary(lm(simul2$ratio~simul2$origin_dist))
# med = signif(median(simul2$ratio), digits=3)
# legend("topright", legend=c(paste0("Med = ", med), "y = -0.26 - 3.3e-8 x"), bty="n")
