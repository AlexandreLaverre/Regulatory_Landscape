
####################################################################################################################
############################################ Résultats Alignements TBA  ############################################
####################################################################################################################
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/human2mouse/")
data_nonull <- read.table("AlignmentStatistics_TBA_human2mouse_withoutnull.txt", header=T)
baited_nonull <- read.table("Align_baited_fragments_withoutnull.txt", header=F)
contacted_nonull <- read.table("Align_contacted_fragments_withoutnull.txt", header=F)
contacted_simul_nonull <- read.table("Align_contacted_simul_withoutnull.txt", header=F)

par(mfrow=c(2,1))
hist(data_nonull$LengthIdentical/data_nonull$TotalLength.human, main="Whitout no_align", xlab="Alignement score (%)")
abline(v = mean(data_nonull$LengthIdentical/data_nonull$TotalLength.human), col="red", lwd=3, lty=2)
mean = signif(mean(data_nonull$LengthIdentical/data_nonull$TotalLength.human), digits=3)
length = length(data_nonull$LengthIdentical/data_nonull$TotalLength.human)
legend("topright", legend=paste0("Mean = ",mean, "\nN = ", length), bty="n")

hist(contacted_nonull$V6/contacted_nonull$V3, main="Contacted fragments", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = median(contacted_nonull$V6/contacted_nonull$V3), col="red", lwd=3, lty=2)
med = signif(mean(contacted_nonull$V6/contacted_nonull$V3), digits=3)
length = length(contacted_nonull$V6/contacted_nonull$V3)
legend("topright", legend=c(paste0("Med = ",med, "\nN = ", length)), bty="n")

hist(contacted_simul_nonull$V6/contacted_simul_nonull$V3, main="Contacted simul fragments", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = median(contacted_simul_nonull$V6/contacted_simul_nonull$V3), col="red", lwd=3, lty=2)
med = signif(median(contacted_simul_nonull$V6/contacted_simul_nonull$V3), digits=3)
length = length(contacted_simul_nonull$V6/contacted_simul_nonull$V3)
legend("topright", legend=c(paste0("Med = ",med, "\nN = ", length)), bty="n")

hist(baited_nonull$V6/baited_nonull$V3, main="Baited fragments", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = median(baited_nonull$V6/baited_nonull$V3), col="red", lwd=3, lty=2)
med = signif(mean(baited_nonull$V6/baited_nonull$V3), digits=3)
length = length(baited_nonull$V6/baited_nonull$V3)
legend("topright", legend=c(paste0("Med = ",med, "\nN = ", length)), bty="n")

# Compare boxplot distribution
par(mfrow=c(1,1))
par(mgp = c(3, 1.5,0))
length_all = paste0("\n All frag \n N = ", length(data_nonull$LengthIdentical/data_nonull$TotalLength.human))
length_cont_simul = paste0("\n simul frag \n N = ", length(contacted_simul_nonull$V6/contacted_simul_nonull$V3))
length_cont = paste0("\n contacted frag \n N = ", length(contacted_nonull$V6/contacted_nonull$V3))
length_bait = paste0("\n baited frag \n N = ", length(baited_nonull$V6/baited_nonull$V3))
boxplot(data_nonull$LengthIdentical/data_nonull$TotalLength.human, contacted_simul_nonull$V6/contacted_simul_nonull$V3,
        contacted_nonull$V6/contacted_nonull$V3, baited_nonull$V6/baited_nonull$V3, 
        notch=T, outline=F, names=c(length_all, length_cont_simul, length_cont, length_bait), ylab="Alignment score",
        main="Human2Mouse")

# Distance ~ conservation
inter_conserv <- read.table("interaction_conserv.txt", header=T)
inter_simul <- read.table("interaction_conserv_simul.txt", header=T)

hist(inter_conserv$PIR_score, main="Mean interactions score ", xlab="align score (%)")
abline(v = median(inter_conserv[!is.na(inter_conserv$PIR_score),]$PIR_score), col="red", lwd=3, lty=2)
med = signif(median(inter_conserv[!is.na(inter_conserv$PIR_score),]$PIR_score), digits=3)
length = length(inter_conserv$PIR_score)
legend("topright", legend=c(paste0("Med = ",med, "\nN = ", length)), bty="n")

hist(inter_simul$PIR_score, main="Mean simul score", xlab="align score (%)")
abline(v = median(inter_simul[!is.na(inter_simul$PIR_score),]$PIR_score), col="red", lwd=3, lty=2)
med = signif(median(inter_simul[!is.na(inter_simul$PIR_score),]$PIR_score), digits=3)
length = length(inter_simul$PIR_score)
legend("topright", legend=c(paste0("Med = ",med, "\nN = ", length)), bty="n")

hist(inter_simul$dist_obs, breaks=50000, xlim=c(0,1000000))
hist(inter_simul$dist_conserv, breaks=50000, xlim=c(0,1000000))

hist(log2(inter_conserv$dist_conserv/inter_conserv$dist_obs), xlim=c(-2,2), breaks=1000)
hist(log2(inter_simul$dist_conserv/inter_simul$dist_obs), xlim=c(-2,2), breaks=1000)

par(mfrow=c(1,1))
boxplot(log2(inter_conserv$dist_conserv/inter_conserv$dist_obs), log2(inter_simul$dist_conserv/inter_simul$dist_obs),
        notch=T, outline=F, names=c("Real", "Simul"), ylab="fold change (log2)", main="Mouse2Human")

inter_conserv$class <-cut(inter_conserv$dist_obs, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
inter_simul$class <-cut(inter_simul$dist_obs, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)

conserv_dist <- data.frame(inter = sapply(levels(inter_conserv$class), function(x) (nrow(inter_conserv[which(inter_conserv$class == x & !is.na(inter_conserv$dist_conserv)),])/ nrow(inter_conserv[which(inter_conserv$class == x),]))*100))
conserv_dist$int_start <- sapply(levels(inter_conserv$class), function(x)  (prop.test(x = nrow(inter_conserv[which(inter_conserv$class == x & !is.na(inter_conserv$dist_conserv)) ,]), n=nrow(inter_conserv[which(inter_conserv$class == x),]), p=0.5)$conf.int[1])*100)
conserv_dist$int_end <- sapply(levels(inter_conserv$class), function(x)  (prop.test(x = nrow(inter_conserv[which(inter_conserv$class == x & !is.na(inter_conserv$dist_conserv)) ,]), n=nrow(inter_conserv[which(inter_conserv$class == x),]), p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(inter_simul$class), function(x) (nrow(inter_simul[which(inter_simul$class == x & !is.na(inter_simul$dist_conserv)),])/ nrow(inter_simul[which(inter_simul$class == x),]))*100))
simul_dist$int_start <- sapply(levels(inter_simul$class), function(x)  (prop.test(x = nrow(inter_simul[which(inter_simul$class == x & !is.na(inter_simul$dist_conserv)) ,]), n=nrow(inter_simul[which(inter_simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(inter_simul$class), function(x)  (prop.test(x = nrow(inter_simul[which(inter_simul$class == x & !is.na(inter_simul$dist_conserv)) ,]), n=nrow(inter_simul[which(inter_simul$class == x),]), p=0.5)$conf.int[2])*100)

plot(conserv_dist$inter[1:50], type="b", col="red",  ylim=c(32,63), main="Human2Mouse", xlab="",ylab="Séquences conservées (%)", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

###### MidDistance ~ Conservation #####
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/mouse2human/")
midist_conserv <- read.table("PIR_midist_conserv.txt", header=T)
midist_simul <- read.table("PIR_midist_conserv_simul.txt", header=T)

midist_conserv$class <-cut(midist_conserv$midist_obs, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
midist_simul$class <-cut(midist_simul$midist_obs, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
class <- c("25-75Kb", "225-275Kb", "475-525Kb", "725-775Kb", "975Kb-1.02Mb", "1.02-1.07Mb", "1.22-1.27Mb", "1.42-1.47Mb", "1.62-1.67Mb", "1.82-1.87Mb", "2.02-2.07Mb", "2.22-2.27Mb", "2.42-2.47Mb", "2.62-2.67Mb", "2.82-2.87Mb", "3.02-3.07Mb")

midist_conserv$class_score <-cut(midist_conserv$PIR_score, breaks=c(0, seq(from=0.001, to=1, by=0.33)), include.lowest = T, include.highest=T)
midist_simul$class_score <-cut(midist_simul$PIR_score, breaks=c(0, seq(from=0.001, to=1, by=0.33)), include.lowest = T, include.highest=T)

midist_conserv$part_exon <- (midist_conserv$exon_pb/midist_conserv$length)*100
midist_conserv$part_exon_class <-cut(midist_conserv$part_exon, breaks=seq(from=0, to=100, by=10), include.lowest = T, include.highest=T)
midist_simul$part_exon <- (midist_simul$exon_pb/midist_simul$length)*100
midist_simul$part_exon_class <-cut(midist_simul$part_exon, breaks=seq(from=0, to=100, by=10), include.lowest = T, include.highest=T)

par(mfrow=c(1,2))
# Nb contact ~ Align score
plot(midist_conserv$PIR_score~midist_conserv$nb_contact, cex=0.5, main = "Observation", ylab="Score align", xlab="Nombre de contact")
plot(midist_simul$PIR_score~midist_simul$nb_contact, cex=0.5, main = "Simulation", ylab="Score align", xlab="Nombre de contact")

# Distance ~ Nb contact
plot(midist_conserv$nb_contact~midist_conserv$midist_obs, cex=0.5, main = "Observation", xlab="Distance médiane", ylab="Nombre de contact")
plot(midist_simul$nb_contact~midist_simul$midist_obs, cex=0.5, main = "Simulation", xlab="Distance médiane", ylab="Nombre de contact")

# Distance ~ Align score
plot(midist_conserv$PIR_score~midist_conserv$midist_obs, cex=0.5, main = "Observation", xlab="Distance médiane", ylab="Align score")
plot(midist_simul$PIR_score~midist_simul$midist_obs, cex=0.5, main = "Simulation", xlab="Distance médiane", ylab="Align score")

# Distance ~ Align score moyen
par(mfrow=c(1,1))
real <- data.frame(mean = sapply(levels(midist_conserv$class), function(x) mean(midist_conserv[which(midist_conserv$class == x & !is.na(midist_conserv$PIR_score)),]$PIR_score)))
simul <- data.frame(mean = sapply(levels(midist_simul$class), function(x) mean(midist_simul[which(midist_simul$class == x & !is.na(midist_simul$PIR_score)),]$PIR_score)))
plot(real$mean[1:80], type='l', cex=0.5, main="Human2Mouse", xlab="", ylab="Score align moyen", ylim=c(0.14,0.30), xaxt = "n")
points(simul$mean[1:80], col='red', type='l', cex=0.5)
axis(1, at=seq(1,80,5), labels=F)
text(seq(1,80,5), par("usr")[3]-0.007, labels = class, srt = 45, pos = 1, xpd = TRUE, cex=0.8)

par(mfrow=c(1,1))
conserv_dist <- data.frame(inter = sapply(levels(midist_conserv$class), function(x) (nrow(midist_conserv[which(midist_conserv$class == x & midist_conserv$PIR_score >0),])/ nrow(midist_conserv[which(midist_conserv$class == x),]))*100))
conserv_dist$int_start <- sapply(levels(midist_conserv$class), function(x)  (prop.test(x = nrow(midist_conserv[which(midist_conserv$class == x & midist_conserv$PIR_score >0),]), n=nrow(midist_conserv[which(midist_conserv$class == x),])+1, p=0.5)$conf.int[1])*100)
conserv_dist$int_end <- sapply(levels(midist_conserv$class), function(x)  (prop.test(x = nrow(midist_conserv[which(midist_conserv$class == x & midist_conserv$PIR_score >0),]), n=nrow(midist_conserv[which(midist_conserv$class == x),])+1, p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(midist_simul$class), function(x) (nrow(midist_simul[which(midist_simul$class == x & midist_simul$PIR_score >0),])/ nrow(midist_simul[which(midist_simul$class == x),]))*100))
simul_dist$int_start <- sapply(levels(midist_simul$class), function(x)  (prop.test(x = nrow(midist_simul[which(midist_simul$class == x & midist_simul$PIR_score >0),]), n=nrow(midist_simul[which(midist_simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(midist_simul$class), function(x)  (prop.test(x = nrow(midist_simul[which(midist_simul$class == x & midist_simul$PIR_score >0),]), n=nrow(midist_simul[which(midist_simul$class == x),]), p=0.5)$conf.int[2])*100)

class <- c("25-75Kb", "225-275Kb", "475-525Kb", "725-775Kb", "975Kb-1.02Mb", "1.02-1.07Mb", "1.22-1.27Mb", "1.42-1.47Mb", "1.62-1.67Mb")

plot(conserv_dist$inter[1:40], type="b", col="red", ylim=c(40,70), main="Human2Mouse", xlab="",ylab="Séquences conservées (%)", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:40,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:40], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:40,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,41,5), labels=F)
text(seq(1,41,5), par("usr")[3] - 1.6, labels = class, srt = 45, pos = 1, xpd = TRUE, cex=0.8)

# Conservation ~ Distance ~ Score alignement
library(forcats)
library(ggplot2)

midist_conserv$cat <- as.factor(paste(midist_conserv$class, midist_conserv$class_score))
midist_conserv$nb_cat <- table(midist_conserv$cat)[midist_conserv$cat]
midist_conserv$nb_class <- table(midist_conserv$class)[midist_conserv$class]
midist_conserv$prop_conserv <- ((midist_conserv$nb_cat/midist_conserv$nb_class)*100)

midist_simul$cat <- as.factor(paste(midist_simul$class, midist_simul$class_score))
midist_simul$nb_cat <- table(midist_simul$cat)[midist_simul$cat]
midist_simul$nb_class <- table(midist_simul$class)[midist_simul$class]
midist_simul$prop_conserv <- ((midist_simul$nb_cat/midist_simul$nb_class)*100)

conserv <- midist_conserv[!duplicated(midist_conserv$cat),c("class","class_score","prop_conserv")]
conserv$data <-  "Observed"
conserv_simul <- midist_simul[!duplicated(midist_simul$cat),c("class","class_score","prop_conserv")]
conserv_simul$data <-"Simulated"

result <- rbind(conserv, conserv_simul)
result <- subset(result, result$class %in% levels(midist_conserv$class)[1:50])
result <- subset(result, !is.na(result$class_score))

ggplot(result, aes(x=class, y=prop_conserv, group = paste(class_score, data))) + ylab("Seq. conservée (%)") + xlab("") +
  geom_line(aes(linetype=data, color=class_score)) + theme_minimal() + ylim(0,60) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, 
                                   color=rep(c("black", "transparent","transparent","transparent", "transparent", "transparent"), length(result$class)/6)))

####################################################################################################################
##################################### Comparaison simulations vs observations ######################################
####################################################################################################################
library(ggplot2)
setwd("/home/laverre/Documents/Regulatory_Landscape/data/mouse/Simulations/")
dist_real <- read.table("dist_real.txt")
dist_simul2 <- read.table("dist_simulations_20-10Mb_bin5kb_fragoverbin.txt")
dist_simul <- read.table("dist_simulations_mouse_10Mb_bin5kb_fragoverbin.txt")
cumul_real <- ecdf(dist_real$V1)
cumul_siml <- ecdf(dist_simul2$V1)


par(mfrow=c(1,2))
boxplot(dist_real$V1, dist_simul2$V1, main="Frag over bin",
        notch=T, outline=F, border=c("black", "red"), names=c("Real", "Simul"))
plot(cumul_real, main="Cumulative Distribution Function", xlab="Distance (pb)")
lines(cumul_siml, col="red")
legend("bottomright", legend=c( "Wasserstein \n = 118 671\n\n\n\n"), bty="n")


# Distribution
real <- data.frame(Distance=dist_real$V1)
simul <- data.frame(Distance=dist_simul$V1)
simul2 <- data.frame(Distance=dist_simul2$V1)
real$Data <- 'Real'
simul$Data <- 'Simul'
simul2$Data <- 'Simul_new'
dataLengths <- rbind(real, simul, simul2)

ggplot(dataLengths, aes(Distance, fill = Data)) + 
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 1000) + xlim(c(-1000000, +1000000)) +
  theme_minimal() + ggtitle("Comparaison distribution Mouse") 


####################################################################################################################
##################################### Corrélation activité peaks/enhancers #########################################
####################################################################################################################
setwd("/home/laverre/Documents/Regulatory_Landscape/result/expression_correlation")
cor_mouse <- read.table("expression_correlation_mouse.txt", header=T)
cor_human <- read.table("expression_correlation_human.txt", header=T)
cor_mouse_valid <- read.table("expression_correlation_mouse_valid.txt", header=T)

par(mfrow=c(2,1))
hist(cor_human$spearman_cor, xlim=c(-0.4, 1), breaks=50, main="Human - peaks-enh activity correlation", xlab="Spearman")
abline(v = median(cor_human$spearman_cor), col="red", lwd=3, lty=2)
hist(cor_human$pearson_cor, xlim=c(-0.4, 1), breaks=100, main="", xlab="Pearson")
abline(v = median(cor_human$pearson_cor), col="red", lwd=3, lty=2)


hist(cor_mouse$spearman_cor, xlim=c(-0.4, 1), breaks=50, main="mouse - peaks-enh activity correlation", xlab="Spearman")
abline(v = median(cor_mouse$spearman_cor), col="red", lwd=3, lty=2)
hist(cor_mouse_valid$spearman_cor, xlim=c(-0.4, 1), breaks=50, main="mouse - peaks-enh activity correlation", xlab="Spearman")
abline(v = median(cor_mouse_valid$spearman_cor), col="red", lwd=3, lty=2)

hist(cor_mouse$pearson_cor, xlim=c(-0.4, 1), breaks=100, main="", xlab="Pearson")
abline(v = median(cor_mouse$pearson_cor), col="red", lwd=3, lty=2)
hist(cor_mouse_valid$pearson_cor, xlim=c(-0.4, 1), breaks=100, main="", xlab="Pearson")
abline(v = median(cor_mouse_valid$pearson_cor), col="red", lwd=3, lty=2)



vect_enh <- read.table("vect_enh")
vect_peak <- read.table("vect_peak")
vect_enh <- t(vect_enh[,-1])
vect_peak <- t(vect_peak[,-c(1:15)])
cor <- as.data.frame(cbind(vect_enh, vect_peak))
colnames(cor) <- c("vect_enh", "vect_peak")

par(mfrow=c(1,1))
plot(log(vect_enh+1), log(vect_peak+1))
abline(lm(log(vect_enh+1) ~ log(vect_peak+1)))

plot(log(cor[cor$vect_enh != 0,]$vect_enh+1),log(cor[cor$vect_enh != 0,]$vect_peak+1))
abline(lm(log(cor[cor$vect_enh != 0,]$vect_enh+1) ~ log(cor[cor$vect_enh != 0,]$vect_peak+1)))


plot(log(cor[cor$vect_peak != 0 & cor$vect_enh != 0,]$vect_enh+1),log(cor[cor$vect_peak != 0 & cor$vect_enh != 0,]$vect_peak+1))
abline(lm(log(cor[cor$vect_peak != 0 & cor$vect_enh != 0,]$vect_enh+1) ~ log(cor[cor$vect_peak != 0 & cor$vect_enh != 0,]$vect_peak+1)))

cor(log(cor[cor$vect_enh != 0,]$vect_enh+1), log(cor[cor$vect_enh != 0,]$vect_peak+1), method = "pearson")
