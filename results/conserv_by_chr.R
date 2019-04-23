library(ggplot2)
library(gridExtra)
setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

sp_origin = 'human'
sp_target = 'mouse'

obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction.txt2", sep=""), header=T)
simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction_simul.txt2", sep=""), header=T)
intergenic <- read.table(paste("intergenic_distance_",sp_origin, "_",sp_target,".txt", sep=""), header=T)
colnames(intergenic) <- c("gene", "origin_dist", "target_dist")

obs$conserv_class <-cut(obs$PIR_score, breaks=seq(from=0, to=1, by=0.1), include.lowest = T)
simul$conserv_class <-cut(simul$PIR_score, breaks=seq(from=0, to=1, by=0.1), include.lowest = T)

# Interaction according to chromosome
par(mfrow=c(4,3))
obs$chr <- as.factor(sub(":.*", "", obs$origin_interaction))
obs$origin_bait <- as.factor(sub("-.*", "", obs$origin_interaction))
obs$origin_PIR <- as.factor(sub(".*-", "", obs$origin_interaction))
obs$origin_PIR_end <- as.integer(sub("chr.*:", "", obs$origin_PIR))
obs$origin_bait_end <- as.integer(sub("chr.*:", "", obs$origin_bait))

color_pallete_function <- colorRampPalette(colors = c("blue", "violet", "orange", "red"), space = "Lab")
conserv_color <- color_pallete_function(nlevels(obs$chr))

for (each_chr in levels(obs$chr)){
  plot(obs[which(obs$chr == each_chr),]$bait_score~obs[which(obs$chr == each_chr),]$origin_dist, cex=0.1, 
       ylab="Bait score", xlab="Origin_dist", main=each_chr)
  #plot(obs[which(obs$chr == each_chr),]$origin_bait_end~obs[which(obs$chr == each_chr),]$origin_PIR_end, cex=0.1, 
  #    ylab='Bait position', xlab="PIR position", main=each_chr)
  #abline(a=0,b=1)
}

plot(obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]$bait_score~obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]$origin_dist, cex=0.1, 
     ylab="Bait score", xlab="Origin_dist", main="chr10", xlim=c(0,10000000))
plot(obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]$origin_bait_end~obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]$origin_PIR_end, cex=0.1, 
     ylab='Bait position', xlab="PIR position", main="chr10")
abline(a=0,b=1)

data_chr10 <- obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]
data_chr10$origin_bait <- as.factor(sub("-.*", "", data_chr10$origin_interaction))

#### Test comparaison chromo obs vs simul ###"
par(mfrow=c(2,1))
simul$chr <- as.factor(sub(":.*", "", simul$origin_interaction))
simul$origin_bait <- as.factor(sub("-.*", "", simul$origin_interaction))
simul$origin_PIR <- as.factor(sub(".*-", "", simul$origin_interaction))
simul$origin_PIR_end <- as.integer(sub("chr.*:", "", simul$origin_PIR))
simul$origin_bait_end <- as.integer(sub("chr.*:", "", simul$origin_bait))
color_pallete_function <- colorRampPalette(colors = c("blue", "violet", "orange", "red"), space = "Lab")
conserv_color <- color_pallete_function(nlevels(simul$chr))

par(mfrow=c(1,2))
for (each_chr in levels(obs$chr)){
  plot(obs[which(obs$chr == each_chr),]$origin_bait_end~obs[which(obs$chr == each_chr),]$origin_PIR_end, cex=0.1, 
       ylab='Bait position', xlab="PIR position", main=c(each_chr, "observation"))
  legend("topleft", legend = c('Distance médiane:',median(obs[which(obs$chr == each_chr),]$origin_dist)), bty='n')
  abline(a=0,b=1)
  
  plot(simul[which(simul$chr == each_chr),]$origin_bait_end~simul[which(simul$chr == each_chr),]$origin_PIR_end, cex=0.1, 
       ylab='Bait position', xlab="PIR position", main=c(each_chr, "simulation"))
  legend("topleft", legend = c('Distance médiane:',median(simul[which(simul$chr == each_chr),]$origin_dist)), bty='n')
  abline(a=0,b=1)
}
