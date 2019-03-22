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
midist_conserv <- read.table("PIR_midist_conserv_0.4.txt", header=T)
midist_simul <- read.table("PIR_midist_conserv_simul_0.4.txt", header=T)

midist_conserv$class <-cut(midist_conserv$midist_obs, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
midist_simul$class <-cut(midist_simul$midist_obs, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
class <- c("25-75Kb", "225-275Kb", "475-525Kb", "725-775Kb", "975Kb-1.02Mb", "1.02-1.07Mb", "1.22-1.27Mb", "1.42-1.47Mb", "1.62-1.67Mb", "1.82-1.87Mb", "2.02-2.07Mb", "2.22-2.27Mb", "2.42-2.47Mb", "2.62-2.67Mb", "2.82-2.87Mb", "3.02-3.07Mb")

midist_conserv$class_score <-cut(midist_conserv$PIR_score, breaks=seq(from=0.01, to=1, by=0.33), include.lowest = T)
midist_simul$class_score <-cut(midist_simul$PIR_score, breaks=seq(from=0.01, to=1, by=0.33), include.lowest = T)

midist_conserv$class_score <- as.character(midist_conserv$class_score)
midist_conserv$class_score[is.na(midist_conserv$class_score)] <- 0
midist_conserv$class_score <- as.factor(midist_conserv$class_score)
midist_simul$class_score <- as.character(midist_simul$class_score)
midist_simul$class_score[is.na(midist_simul$class_score)] <- 0
midist_simul$class_score <- as.factor(midist_simul$class_score)

midist_conserv$part_all_exon <- (midist_conserv$all_exon_pb/midist_conserv$length)*100
midist_conserv$part_all_exon_class <-cut(midist_conserv$part_all_exon, breaks=seq(from=0, to=100, by=20), include.lowest = T, include.highest=T)
midist_simul$part_all_exon <- (midist_simul$all_exon_pb/midist_simul$length)*100
midist_simul$part_all_exon_class <-cut(midist_simul$part_all_exon, breaks=seq(from=0, to=100, by=20), include.lowest = T, include.highest=T)

midist_conserv$part_coding_exon <- (midist_conserv$coding_exon_pb/midist_conserv$length)*100
midist_conserv$part_coding_exon_class <-cut(midist_conserv$part_coding_exon, breaks=seq(from=0, to=100, by=20), include.lowest = T, include.highest=T)
midist_simul$part_coding_exon <- (midist_simul$coding_exon_pb/midist_simul$length)*100
midist_simul$part_coding_exon_class <-cut(midist_simul$part_coding_exon, breaks=seq(from=0, to=100, by=20), include.lowest = T, include.highest=T)

midist_conserv$part_nocoding_exon <- (midist_conserv$nocoding_exon_pb/midist_conserv$length)*100
midist_conserv$part_nocoding_exon_class <-cut(midist_conserv$part_nocoding_exon, breaks=seq(from=0, to=100, by=20), include.lowest = T, include.highest=T)
midist_simul$part_nocoding_exon <- (midist_simul$nocoding_exon_pb/midist_simul$length)*100
midist_simul$part_nocoding_exon_class <-cut(midist_simul$part_nocoding_exon, breaks=seq(from=0, to=100, by=20), include.lowest = T, include.highest=T)

midist_conserv$part_repeat <- (midist_conserv$nocoding_exon_pb/midist_conserv$length)*100
midist_conserv$part_repeat_class <-cut(midist_conserv$part_nocoding_exon, breaks=seq(from=0, to=100, by=20), include.lowest = T, include.highest=T)
midist_simul$part_repeat <- (midist_simul$nocoding_exon_pb/midist_simul$length)*100
midist_simul$part_repeat_class <-cut(midist_simul$part_nocoding_exon, breaks=seq(from=0, to=100, by=20), include.lowest = T, include.highest=T)

par(mfrow=c(1,2))
# % exon ~ align_score
plot(midist_conserv$PIR_score~midist_conserv$part_all_exon, cex=0.2, main = "Observation", ylab="Score align", xlab="% pb all_exon")
abline(lm(midist_conserv$PIR_score~midist_conserv$part_all_exon), col="red")
plot(midist_simul$PIR_score~midist_simul$part_all_exon, cex=0.2, main = "Simulation", ylab="Score align", xlab="%pb all_exon")
abline(lm(midist_simul$PIR_score~midist_simul$part_all_exon), col="red")

par(mfrow=c(1,1))
# % exon_class ~ align_score
plot(midist_conserv$PIR_score~midist_conserv$part_all_exon, cex=0.2, main = "All exons", ylab="Score align", xlab="% pb all_exon")
abline(lm(midist_conserv$PIR_score~midist_conserv$part_all_exon), col="red")

plot(midist_conserv$PIR_score~midist_conserv$part_coding_exon, cex=0.2, main = "Coding exons", ylab="Score align", xlab="%pb exon")
abline(lm(midist_conserv$PIR_score~midist_conserv$part_coding_exon), col="red")

plot(midist_conserv$PIR_score~midist_conserv$part_nocoding_exon, cex=0.2, main = "No coding exons", ylab="Score align", xlab="%pb exon")
abline(lm(midist_conserv$PIR_score~midist_conserv$part_nocoding_exon), col="red")

plot(midist_conserv$PIR_score~midist_conserv$part_repeat, cex=0.2, main = "Repeat elements", ylab="Score align", xlab="%pb exon")
abline(lm(midist_conserv$PIR_score~midist_conserv$part_repeat), col="red")

# % exon ~ distance 
plot(midist_conserv$part_all_exon~midist_conserv$midist_obs, cex=0.2,  col=midist_conserv$chr, main = "Observation", xlab="Distance médiane", ylab="% pb exons")

plot(midist_simul$part_all_exon~midist_simul$midist_obs, cex=0.2,  col=midist_simul$chr, main = "Simulation", xlab="Distance médiane", ylab="% pb exons")

plot(midist_conserv$part_coding_exon~midist_conserv$midist_obs, cex=0.5, main = "Coding exons", xlab="Distance médiane", ylab="% pb exons")
plot(midist_conserv$part_nocoding_exon~midist_conserv$midist_obs, cex=0.5, main = "No coding exons", xlab="Distance médiane", ylab="% pb exons")


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

class <- c("25-75Kb", "225-275Kb", "475-525Kb", "725-775Kb", "975Kb-1.02Mb", "1.22-1.27Mb", "1.47-1.52Mb", "1.72-1.77Mb", "1.97-2.02Mb")
class <- c("25-75Kb", "475-525Kb", "975Kb-1.02Mb", "1.47-1.52Mb", "1.97-2.02Mb","2.47-2.52Mb")

plot(conserv_dist$inter[1:40], type="b", col="red", ylim=c(30,55), main="Mouse2Human 0.1", xlab="",ylab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:40,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:40], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:40,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,41,10), labels=F)
text(seq(1,41,10), par("usr")[3] - 2, labels = class, srt = 45, pos = 1, xpd = TRUE, cex=0.8)

legend("bottomleft", cex=1.2, fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

# Conservation ~ Distance ~ Score alignement
library(forcats)
library(ggplot2)

midist_conserv$cat <- as.factor(paste(midist_conserv$class, midist_conserv$part_exon_class))
midist_conserv$nb_cat <- table(midist_conserv$cat)[midist_conserv$cat]
midist_conserv$nb_class <- table(midist_conserv$class)[midist_conserv$class]
midist_conserv$prop_conserv <- ((midist_conserv$nb_cat/midist_conserv$nb_class)*100)

midist_simul$cat <- as.factor(paste(midist_simul$class, midist_simul$part_exon_class))
midist_simul$nb_cat <- table(midist_simul$cat)[midist_simul$cat]
midist_simul$nb_class <- table(midist_simul$class)[midist_simul$class]
midist_simul$prop_conserv <- ((midist_simul$nb_cat/midist_simul$nb_class)*100)

conserv <- midist_conserv[!duplicated(midist_conserv$cat),c("class","part_exon_class","prop_conserv")]
conserv$data <-  "Observed"
conserv_simul <- midist_simul[!duplicated(midist_simul$cat),c("class","part_exon_class","prop_conserv")]
conserv_simul$data <-"Simulated"

result <- rbind(conserv, conserv_simul)
result <- subset(result, result$class %in% levels(midist_conserv$class)[1:50])
result <- subset(result, !is.na(result$part_exon_class))

ggplot(result, aes(x=class, y=prop_conserv, group = paste(part_exon_class, data))) + ylab("Seq. conservée (%)") + xlab("") +
  geom_line(aes(linetype=data, color=part_exon_class)) + theme_minimal() + ylim(0,10) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, 
                                   color=rep(c("black", "transparent","transparent","transparent", "transparent", "transparent"), length(result$class)/6)))
