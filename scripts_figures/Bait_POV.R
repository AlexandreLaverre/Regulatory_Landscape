setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/")

sp_origin = "human"
sp_target = "mouse"

# BAIT
obs_bait <- read.table(paste("Bait_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,"_merged.txt", sep=""), header=T)
obs_nomerged_bait <- read.table(paste("Bait_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,"_nomerged.txt", sep=""), header=T)
simul_bait <- read.table(paste("Bait_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,"_simul_merged.txt", sep=""), header=T)

#### Filtres ### ----------------------------------
# Length
obs_bait$length <- obs_bait$end - obs_bait$start
obs_bait <- obs_bait[which(obs_bait$length > 250 & obs_bait$length < 20000),]
obs_nomerged_bait$length <- obs_nomerged_bait$end - obs_nomerged_bait$start
obs_nomerged_bait <- obs_nomerged_bait[which(obs_nomerged_bait$length > 250 & obs_nomerged_bait$length < 20000),]
simul_bait$length <- simul_bait$end - simul_bait$start
simul_bait <- simul_bait[which(simul_bait$length > 250 & simul_bait$length < 20000),]
# Duplication
obs_bait <- obs_bait[which(obs_bait$duplication == 0),]
obs_nomerged_bait <- obs_nomerged_bait[which(obs_nomerged_bait$duplication == 0),]
simul_bait <- simul_bait[which(simul_bait$duplication == 0),]


## Classe nb contact
obs_bait_nbbaitcontacts=cut(obs_bait$nb_Baitcontact, breaks=c(1,5,10, 15,20,25,30, max(obs_bait$nb_Baitcontact)), include.lowest=T)
obs_bait_nbPIRcontacts=cut(obs_bait$nb_PIRcontact, breaks=c(1,2, 5, 10, max(obs_bait$nb_PIRcontact)), include.lowest=T)
obs_bait$Allexon_ungapped <- obs_bait$all_ungapped/obs_bait$all_length

obs_nomerged_bait_nbbaitcontacts=cut(obs_nomerged_bait$nb_Baitcontact, breaks=c(1,5,10, 15,20,25,30, max(obs_nomerged_bait$nb_Baitcontact)), include.lowest=T)
obs_nomerged_bait_nbPIRcontacts=cut(obs_nomerged_bait$nb_PIRcontact, breaks=c(1,2, 5, 10, max(obs_nomerged_bait$nb_PIRcontact)), include.lowest=T)
obs_nomerged_bait$Allexon_ungapped <- obs_nomerged_bait$all_ungapped/obs_nomerged_bait$all_length

simul_bait_nbbaitcontacts=cut(simul_bait$nb_Baitcontact, breaks=c(1,5,10, 15,20,25,30, max(simul_bait$nb_Baitcontact)), include.lowest=T)
simul_bait_nbPIRcontacts=cut(simul_bait$nb_PIRcontact, breaks=c(1,2, 5, 10, max(simul_bait$nb_PIRcontact)), include.lowest=T)
simul_bait$Allexon_ungapped <- simul_bait$all_ungapped/simul_bait$all_length

# Conserv ~ nb contact
boxplot(obs_bait$Allexon_ungapped~obs_bait_nbPIRcontacts, outline=F, notch = T, xlab="Other contact number", ylab="Score conserv Allexon ungapped", main="Bait")
boxplot(simul_bait$Allexon_ungapped~simul_bait_nbPIRcontacts, outline=F, notch = T, xlab="Other contact number", ylab="Score conserv Allexon ungapped", main="Bait")

boxplot(obs_bait$Allexon_ungapped~obs_bait_nbbaitcontacts, outline=F, notch = T, xlab="Bait contact number", ylab="Score conserv Allexon ungapped", main="Bait obs")
boxplot(simul_bait$Allexon_ungapped~simul_bait_nbbaitcontacts, outline=F, notch = T, xlab="Bait contact number", ylab="Score conserv Allexon ungapped", main="Bait simul")


################  Dist ~ Nb contact decile ################ 
obs_bait_nbbaitcontacts=cut(obs_bait$nb_Baitcontact, breaks=quantile(obs_bait$nb_Baitcontact, prob=seq(0,1, length=10)), include.lowest=T)
obs_nomerged_bait_nbbaitcontacts=cut(obs_nomerged_bait$nb_Baitcontact, breaks=quantile(obs_nomerged_bait$nb_Baitcontact, prob=seq(0,1, length=10)), include.lowest=T)
simul_bait_nbbaitcontacts=cut(simul_bait$nb_Baitcontact, breaks=quantile(simul_bait$nb_Baitcontact, prob=seq(0,1, length=10)), include.lowest=T)

test_obs <- boxplot(obs_bait$midist_obs~obs_bait_nbbaitcontacts, outline=F, notch = T, xlab="Bait contact number", ylab="Midist", main="Bait", plot=F)
test_nomerged_obs <- boxplot(obs_nomerged_bait$midist_obs~obs_nomerged_bait_nbbaitcontacts, outline=F, notch = T, xlab="Bait contact number", ylab="Midist", main="Bait", plot=F)
test_simul <- boxplot(simul_bait$midist_obs~simul_bait_nbbaitcontacts, outline=F, notch = T, xlab="Bait contact number", ylab="Midist", main="Bait", plot=F)

names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")
plot(test_obs$stats[3,], type='b', cex=0.7, col='red', xaxt='n', cex.lab=1.2, ylim=c(250000, 850000),
     xlab='Contact number decile', ylab="Median distance to contacted regions", main="Bait complexity")
points(test_simul$stats[3,], type='b', cex=0.7, col='navy')
points(test_nomerged_obs$stats[3,], type='b', cex=0.7, col='orange')
for (row in 1:length(test_obs$conf[1,])){
  segments(x0=row,y0=test_obs$conf[1,row],x1=row,y1=test_obs$conf[2,row], col='red')}
for (row in 1:length(test_nomerged_obs$conf[1,])){
  segments(x0=row,y0=test_nomerged_obs$conf[1,row],x1=row,y1=test_nomerged_obs$conf[2,row], col='orange')}
for (row in 1:length(test_simul$conf[1,])){
  segments(x0=row,y0=test_simul$conf[1,row],x1=row,y1=test_simul$conf[2,row], col='navy', lwd=0.3)}
axis(1, at=seq(1,9,1), labels=F)
text(seq(1,9,1),par("usr")[3]-30000, names, xpd = TRUE, cex=1.2)
legend("topleft", fill=c("red","orange", "navy"), legend = c("Observed", "No merged", "Simulated"), bty='n')

# Bait-bait vs Bait-other
obs_bait$nb_baited_contact <- obs_bait$nb_Baitcontact - obs_bait$nb_unbaited_contact
obs_bait_nbbaitedcontacts=cut(obs_bait$nb_baited_contact, breaks=quantile(obs_bait$nb_baited_contact, prob=seq(0,1, length=6)), include.lowest=T)
obs_bait_nbunbaitcontacts=cut(obs_bait$nb_unbaited_contact, breaks=quantile(obs_bait$nb_unbaited_contact, prob=seq(0,1, length=6)), include.lowest=T)

test_obs <- boxplot(obs_bait$midist_obs~obs_bait_nbbaitedcontacts, outline=F, notch = T, xlab="Bait contact number", ylab="Midist", main="Bait", plot=F)
test_obs_unbaited <- boxplot(obs_bait$midist_obs~obs_bait_nbunbaitcontacts, outline=F, notch = T, xlab="Bait contact number", ylab="Midist", main="Bait", plot=F)
names <- c("1", "2", "3", "4", "5")
plot(test_obs$stats[3,], type='b', cex=0.7, col='red', xaxt='n', cex.lab=1.2, ylim=c(250000, 600000),
     xlab='Contact number quintile', ylab="Median distance to contacted regions", main="Bait complexity")
points(test_obs_unbaited$stats[3,], type='b', cex=0.7, col='forestgreen')
for (row in 1:length(test_obs$conf[1,])){
  segments(x0=row,y0=test_obs$conf[1,row],x1=row,y1=test_obs$conf[2,row], col='red')}
for (row in 1:length(test_obs_unbaited$conf[1,])){
  segments(x0=row,y0=test_obs_unbaited$conf[1,row],x1=row,y1=test_obs_unbaited$conf[2,row], col='forestgreen')}
axis(1, at=seq(1,6,1), labels=F)
text(seq(1,6,1),par("usr")[3]-20000, names, xpd = TRUE, cex=1.2)
legend("topleft", fill=c("red","forestgreen"), legend = c("Bait-other", "Bait-bait"), bty='n')



################ Nb contact ~ dist ################ 
obs_bait$class <-cut(obs_bait$midist_obs, breaks=seq(from=25000, to=2500000, by=50000), include.lowest = T)
obs_nomerged_bait$class <-cut(obs_nomerged_bait$midist_obs, breaks=seq(from=25000, to=2500000, by=50000), include.lowest = T)
simul_bait$class <-cut(simul_bait$midist_obs, breaks=seq(from=25000, to=2500000, by=50000), include.lowest = T)
class_leg <- c("50Kb", "500Kb", "1Mb", "1.5Mb", "2Mb","2.5Mb")

obs_dist <- data.frame(inter = sapply(levels(obs_bait$class), function(x) mean(obs_bait[which(obs_bait$class == x),]$nb_Baitcontact)))
obs_dist$int_start <- sapply(levels(obs_bait$class), function(x) t.test(obs_bait[which(obs_bait$class == x),]$nb_Baitcontact)[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(obs_bait$class), function(x) t.test(obs_bait[which(obs_bait$class == x),]$nb_Baitcontact)[["conf.int"]][2])

obs_nomerged_dist <- data.frame(inter = sapply(levels(obs_nomerged$class), function(x) mean(obs_nomerged[which(obs_nomerged$class == x),]$nb_Baitcontact)))
obs_nomerged_dist$int_start <- sapply(levels(obs_nomerged$class), function(x) t.test(obs_nomerged[which(obs_nomerged$class == x),]$nb_Baitcontact)[["conf.int"]][1])
obs_nomerged_dist$int_end <- sapply(levels(obs_nomerged$class), function(x) t.test(obs_nomerged[which(obs_nomerged$class == x),]$nb_Baitcontact)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul_bait$class), function(x) mean(simul_bait[which(simul_bait$class == x),]$nb_Baitcontact)))
simul_dist$int_start <- sapply(levels(simul_bait$class), function(x) t.test(simul_bait[which(simul_bait$class == x),]$nb_Baitcontact)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul_bait$class), function(x) t.test(simul_bait[which(simul_bait$class == x),]$nb_Baitcontact)[["conf.int"]][2])

plot(obs_nomerged_dist$inter[1:50], type="l", col="orange", cex=0.7, cex.lab=1.3,
     main="Bait complexity", xlab="", ylab="Mean nb contact", xaxt = "n")
for (row in 1:nrow(obs_nomerged_dist[1:50,])){
  segments(x0=row,y0=obs_nomerged_dist[row,]$int_start,x1=row,y1=obs_nomerged_dist[row,]$int_end, col='orange', lwd=0.3)}

points(obs_dist$inter[1:50], type="l", col="red", cex=0.7, cex.lab=1.3)
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="l", col="navy", cex=0.7, cex.lab=1.3)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='navy', lwd=0.3)}
axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-15, class_leg, xpd = TRUE, cex=1.2)
legend("topleft", fill=c("red","orange", "navy"), legend = c("Observed", "No merged", "Simul"), bty='n')


## Bait-bait vs Bait-other
obs_bait$nb_baited_contact <- obs_bait$nb_Baitcontact - obs_bait$nb_unbaited_contact

obs_dist <- data.frame(inter = sapply(levels(obs_bait$class), function(x) mean(obs_bait[which(obs_bait$class == x),]$nb_baited_contact)))
obs_dist$int_start <- sapply(levels(obs_bait$class), function(x) t.test(obs_bait[which(obs_bait$class == x),]$nb_baited_contact)[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(obs_bait$class), function(x) t.test(obs_bait[which(obs_bait$class == x),]$nb_baited_contact)[["conf.int"]][2])

obs_dist_other <- data.frame(inter = sapply(levels(obs_bait$class), function(x) mean(obs_bait[which(obs_bait$class == x),]$nb_unbaited_contact)))
obs_dist_other$int_start <- sapply(levels(obs_bait$class), function(x) t.test(obs_bait[which(obs_bait$class == x),]$nb_unbaited_contact)[["conf.int"]][1])
obs_dist_other$int_end <- sapply(levels(obs_bait$class), function(x) t.test(obs_bait[which(obs_bait$class == x),]$nb_unbaited_contact)[["conf.int"]][2])

plot(obs_dist$inter[1:50], type="l", col="forestgreen", cex=0.7, cex.lab=1.3, ylim=c(0,50),
     main="Bait complexity", xlab="", ylab="Mean nb contact", xaxt = "n")
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='forestgreen', lwd=0.3)}

points(obs_dist_other$inter[1:50], type="l", col="red", cex=0.7, cex.lab=1.3, ylim=c(0,70))
for (row in 1:nrow(obs_dist_other[1:50,])){
  segments(x0=row,y0=obs_dist_other[row,]$int_start,x1=row,y1=obs_dist_other[row,]$int_end, col='red', lwd=0.3)}
axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-5, class_leg, xpd = TRUE, cex=1.2)
legend("topleft", fill=c("red","forestgreen"), legend = c("Bait-other", "Bait-bait"), bty='n')


# Contact enhancer
obs_bait$CAGE_portion <- obs_bait$CAGE_contact / obs_bait$nb_Baitcontact
obs_nomerged_bait$CAGE_portion <- obs_nomerged_bait$CAGE_contact / obs_nomerged_bait$nb_Baitcontact
simul_bait$CAGE_portion <- simul_bait$CAGE_contact / simul_bait$nb_Baitcontact

obs_dist <- data.frame(inter = sapply(levels(obs_bait$class), function(x) mean(obs_bait[which(obs_bait$class == x ),]$CAGE_portion)))
obs_dist$int_start <- sapply(levels(obs_bait$class), function(x) t.test(obs_bait[which(obs_bait$class == x),]$CAGE_portion)[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(obs_bait$class), function(x) t.test(obs_bait[which(obs_bait$class == x),]$CAGE_portion)[["conf.int"]][2])

obs_nomerged_dist <- data.frame(inter = sapply(levels(obs_nomerged_bait$class), function(x) mean(obs_nomerged_bait[which(obs_nomerged_bait$class == x ),]$CAGE_portion)))
obs_nomerged_dist$int_start <- sapply(levels(obs_nomerged_bait$class), function(x) t.test(obs_nomerged_bait[which(obs_nomerged_bait$class == x),]$CAGE_portion)[["conf.int"]][1])
obs_nomerged_dist$int_end <- sapply(levels(obs_nomerged_bait$class), function(x) t.test(obs_nomerged_bait[which(obs_nomerged_bait$class == x),]$CAGE_portion)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul_bait$class), function(x) mean(simul_bait[which(simul_bait$class == x ),]$CAGE_portion)))
simul_dist$int_start <- sapply(levels(simul_bait$class), function(x) t.test(simul_bait[which(simul_bait$class == x),]$CAGE_portion)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul_bait$class), function(x) t.test(simul_bait[which(simul_bait$class == x),]$CAGE_portion)[["conf.int"]][2])


plot(obs_dist$inter[1:50], type="l", col="red", cex=0.7, cex.lab=1.3, ylim=c(0,0.3),
     main="Bait", xlab="", ylab="CAGE proportion", xaxt = "n")
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}

points(obs_nomerged_dist$inter[1:50], type="l", col="orange", cex=0.7,)
for (row in 1:nrow(obs_nomerged_dist[1:50,])){
  segments(x0=row,y0=obs_nomerged_dist[row,]$int_start,x1=row,y1=obs_nomerged_dist[row,]$int_end, col='orange', lwd=0.3)}

points(simul_dist$inter[1:50], type="l", col="navy", cex=0.7,)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='navy', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.01, class_leg, xpd = TRUE, cex=1.2)
legend("topright", fill=c("red","orange", "navy"), legend = c("Observed", "No merged", "Simul"), bty='n')


