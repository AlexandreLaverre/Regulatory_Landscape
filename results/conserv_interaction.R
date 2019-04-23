library(ggplot2)
library(gridExtra)
setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

sp_origin = 'mouse'
sp_target = 'human'

obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction.txt", sep=""), header=T)
simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction_simul.txt", sep=""), header=T)
 
obs$conserv_class <-cut(obs$PIR_score, breaks=seq(from=0, to=1, by=0.1), include.lowest = T)
simul$conserv_class <-cut(simul$PIR_score, breaks=seq(from=0, to=1, by=0.1), include.lowest = T)

color_pallete_function <- colorRampPalette(colors = c("blue", "orange", "red"), space = "Lab")
conserv_color <- color_pallete_function(nlevels(obs$conserv_class))

############## Interaction conserv ~ distance ############## 
obs$class <-cut(obs$origin_dist, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
simul$class <- cut(simul$origin_dist, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)

class_leg <- c("25-75Kb", "1.97-2.02Mb", "3.97-4.02Mb","5.97-6.02Mb","7.97-8.02Mb","9.97-10.02Mb")

obs_conserv <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),])/ nrow(obs[which(obs$class == x),]))*100))
obs_conserv$int_start <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x),])+1, p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x),])+1, p=0.5)$conf.int[2])*100)

simul_conserv <- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),])/ nrow(simul[which(simul$class == x),]))*100))
simul_conserv$int_start <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_conserv$int_end <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[2])*100)

obs_conserv$nb <- sapply(levels(obs$class), function(x) nrow(obs[which(obs$class == x),]))

plot(obs_conserv$inter, type="b", col="red", cex=0.7, main=paste(sp_origin,"interactions conserved in", sp_target),
     ylab="Proportion interaction conserv (%)", xlab="", xaxt = "n")
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}
points(simul_conserv$inter, type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_conserv)){
  segments(x0=row,y0=simul_conserv[row,]$int_start,x1=row,y1=simul_conserv[row,]$int_end, col='blue', lwd=0.3)}

#polygon(x=c(seq(1,199,1),seq(199,1,-1)), y=c(obs_conserv$int_start, rev(obs_conserv$int_end)), col=rgb(1,0,0,0.1), border = NA)
#polygon(x=c(seq(1,199,1),seq(199,1,-1)), y=c(simul_conserv$int_start, rev(simul_conserv$int_end)), col=rgb(0,0,1,0.1), border = NA)

axis(1, at=seq(1,201,40), labels=F)
text(seq(1,201,40), par("usr")[3]-0.3, labels = class_leg, pos = 1, xpd = TRUE, cex=0.8)
legend("topright", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

############## Score conserv ~ distance ############## 
obs$class <-cut(obs$origin_dist, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
simul$class <-cut(simul$origin_dist, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
class_leg <- c("25-75Kb", "1.97-2.02Mb", "3.97-4.02Mb","5.97-6.02Mb","7.97-8.02Mb","9.97-10.02Mb")

obs_conserv <- data.frame(inter = sapply(levels(obs$class), function(x) mean(obs[which(obs$class == x),]$PIR_score)))
obs_conserv$int_start <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$PIR_score)[["conf.int"]][1])
obs_conserv$int_end <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$PIR_score)[["conf.int"]][2])

simul_conserv <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$PIR_score)))
simul_conserv$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$PIR_score)[["conf.int"]][1])
simul_conserv$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$PIR_score)[["conf.int"]][2])

plot(obs_conserv$inter, type="b", col="red", cex=0.5, ylim=c(0.1,0.45),
     main=paste(sp_origin, "contacted regions conserved in", sp_target), ylab="Score conserv", xlab="", xaxt = "n")
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}
points(simul_conserv$inter, type="b", col="blue", cex=0.5)
for (row in 1:nrow(simul_conserv)){
  segments(x0=row,y0=simul_conserv[row,]$int_start,x1=row,y1=simul_conserv[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,201,40), labels=F)
text(seq(1,201,40), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE, cex=0.8)
legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

### Relation strength ~ frequency
boxplot(obs$origin_strength~obs$origin_nb_tissu, outline=F, notch=T, 
        xlab="Sample number", ylab="Contact strength (median)", main=paste(sp_origin, "interactions"))

# Conserv ~ frequency
barplot(table(obs$origin_nb_tissu), xlab="Nombre d'échantillon", ylab="Fréquence",main=paste(sp_origin, "interactions"))
obs_conserv <- data.frame(inter = sapply(levels(as.factor(obs$origin_nb_tissu)), function(x) (nrow(obs[which(obs$origin_nb_tissu == x & !is.na(obs$target_dist)),])/ nrow(obs[which(obs$origin_nb_tissu == x),]))*100))
obs_conserv$int_start <- sapply(levels(as.factor(obs$origin_nb_tissu)), function(x)  (prop.test(x = nrow(obs[which(obs$origin_nb_tissu == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$origin_nb_tissu == x),]), p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(as.factor(obs$origin_nb_tissu)), function(x)  (prop.test(x = nrow(obs[which(obs$origin_nb_tissu == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$origin_nb_tissu == x),]), p=0.5)$conf.int[2])*100)

plot(obs_conserv$inter, type="b", col="red", main=paste(sp_origin,"interactions conserved in", sp_target),
     ylab="Proportion interaction conserv (%)", xlab="Nombre d'échantillon")
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}

# Conserv ~ strength
obs$origin_strength <- as.integer(obs$origin_strength)
obs[which(obs$origin_strength >= 20),]$origin_strength <- 20

obs$target_strength <- as.integer(obs$target_strength)
obs[which(obs$target_strength >= 20),]$target_strength <- 20

barplot(table(obs$origin_strength), xlab="Force de contact", ylab="Fréquence", main=paste(sp_origin, "interactions"))

obs_conserv <- data.frame(inter = sapply(levels(as.factor(obs$origin_strength)), function(x) (nrow(obs[which(obs$origin_strength == x & !is.na(obs$target_dist)),])/ nrow(obs[which(obs$origin_strength == x),]))*100))
obs_conserv$int_start <- sapply(levels(as.factor(obs$origin_strength)), function(x)  (prop.test(x = nrow(obs[which(obs$origin_strength == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$origin_strength == x),]), p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(as.factor(obs$origin_strength)), function(x)  (prop.test(x = nrow(obs[which(obs$origin_strength == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$origin_strength == x),]), p=0.5)$conf.int[2])*100)

plot(obs_conserv$inter, type="b", col="red", cex=0.5, main=paste(sp_origin,"interactions conserved in", sp_target),
     ylab="Proportion interaction conserv (%)", xlab="Force du contact (médiane)", xaxt='n', ylim=c(10,35))
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}

axis(1, at=seq(1,16,5), labels=F)
text(seq(1,16,5), par("usr")[3]-0.5, labels = c("5", "10", "15", "20+"), pos = 1, xpd = TRUE)

# Nb tissu target ~ origin
boxplot(obs[which(!is.na(obs$target_nb_tissu)),]$origin_nb_tissu~obs[which(!is.na(obs$target_nb_tissu)),]$target_nb_tissu, 
        outline=F, notch=T, xlab=paste("Nombre d'échantillon", sp_target), ylab=paste("Nombre d'échantillon", sp_origin),
        main=paste(sp_origin, "interactions conserved in", sp_target))

boxplot(obs[which(!is.na(obs$target_strength)),]$origin_strength~obs[which(!is.na(obs$target_strength)),]$target_strength, 
        outline=F, notch=T, xlab=paste("Force de contact", sp_target), ylab=paste("Force de contact", sp_origin),
        main=paste(sp_origin, "interactions conserved in", sp_target))


