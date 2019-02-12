
####################################################################################################################
############################################ Résultats Alignements TBA  ############################################
####################################################################################################################


setwd("/home/laverre/Documents/Regulatory_Landscape/result/human2mouse")
data_nonull <- read.table("AlignmentStatistics_TBA_human2mouse_withoutnull.txt", header=T)
data <- read.table("AlignmentStatistics_TBA_human2mouse.txt", header=T)
baited <- read.table("Align_baited_fragments.txt", header=F)
#baited_nonull <- read.table("Align_baited_fragments_withoutnull.txt", header=F)
contacted <- read.table("Align_contacted_fragments.txt", header=F)
#contacted_nonull <- read.table("Align_contacted_fragments_withoutnull.txt", header=F)

par(mfrow=c(1,2))
hist(data$LengthIdentical/data$TotalLength.human, main="Human to Mouse all_fragments", xlab="Alignement score (%)")
abline(v = mean(data$LengthIdentical/data$TotalLength.human), col="red", lwd=3, lty=2)
mean = signif(mean(data$LengthIdentical/data$TotalLength.human), digits=3)
length = length(data$LengthIdentical/data$TotalLength.human)
legend("topright", legend=paste0("Mean = ",mean, "\nN = ", length), bty="n")

hist(data_nonull$LengthIdentical/data_nonull$TotalLength.human, main="Whitout no_align", xlab="Alignement score (%)")
abline(v = mean(data_nonull$LengthIdentical/data_nonull$TotalLength.human), col="red", lwd=3, lty=2)
mean = signif(mean(data_nonull$LengthIdentical/data_nonull$TotalLength.human), digits=3)
length = length(data_nonull$LengthIdentical/data_nonull$TotalLength.human)
legend("topright", legend=paste0("Mean = ",mean, "\nN = ", length), bty="n")


hist(contacted$V6/contacted$V3, main="Contacted fragments", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = median(contacted$V6/contacted$V3), col="red", lwd=3, lty=2)
med = signif(median(contacted$V6/contacted$V3), digits=3)
length = length(contacted$V6/contacted$V3)
legend("topright", legend=c(paste0("Med = ",med, "\nN = ", length)), bty="n")

hist(baited$V6/baited$V3, main="Baited fragments", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = median(baited$V6/baited$V3), col="red", lwd=3, lty=2)
med = signif(median(baited$V6/baited$V3), digits=3)
length = length(baited$V6/baited$V3)
legend("topright", legend=c(paste0("Med = ",med, "\nN = ", length)), bty="n")

# Compare boxplot distribution
par(mfrow=c(1,1))
par(mgp = c(3, 1.5,0))
length_all = paste0("\n All frag \n N = ", length(data$LengthIdentical/data$TotalLength.human))
length_cont = paste0("\n contacted frag \n N = ", length(contacted$V6/contacted$V3))
length_bait = paste0("\n baited frag \n N = ", length(baited$V6/baited$V3))
boxplot(data$LengthIdentical/data$TotalLength.human, contacted$V6/contacted$V3, baited$V6/baited$V3, 
        notch=T, outline=F, names=c(length_all, length_cont, length_bait), ylab="Alignment score",
        main="Human2Mouse")

# Test violin plot
all <- data.frame(Alignment_score=data$LengthIdentical/data$TotalLength.human)
contact <- data.frame(Alignment_score=contacted$V6/contacted$V3)
bait <- data.frame(Alignment_score=baited$V6/baited$V3)
all$Data <- 'All'
contact$Data <- 'Contacted'
bait$Data <- 'Baited'
dataAlign <- rbind(all, bait, contact)

ggplot(dataAlign, aes(x = Data, y = Alignment_score, fill= Data)) + 
  scale_x_discrete(limits=c("All", "Contacted", "Baited"))+
  geom_violin() + scale_color_brewer(palette="Dark2") + theme_linedraw() +
  geom_boxplot(width=0.05, fill="white", notch=T) +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  ggtitle("Human2Mouse")

####################################################################################################################
##################################### Comparaison simulations vs observations ######################################
####################################################################################################################

library(ggplot2)
setwd("/home/laverre/Documents/Regulatory_Landscape/data/mouse/Simulations/")
dist_real <- read.table("dist_real.txt")
dist_simul <- read.table("dist_simul_simulations_5Mb_bin5kb_uniq10_replace_prob.txt")
cumul_real <- ecdf(dist_real$V1)
cumul_siml <- ecdf(dist_simul$V1)


par(mfrow=c(1,2))
boxplot(dist_real$V1, dist_simul$V1, main="Uniq*10 + proba",
        notch=T, outline=F, border=c("black", "red"), names=c("Real", "Simul"))
plot(cumul_real, main="Cumulative Distribution Function", xlab="Distance (pb)")
lines(cumul_siml, col="red")
legend("bottomright", legend=c( "Wasserstein \n = 310 739\n\n\n\n"), bty="n")


# Distribution
real <- data.frame(Distance=dist_real$V1)
simul <- data.frame(Distance=dist_simul$V1)
real$Data <- 'Real'
simul$Data <- 'Simul'
dataLengths <- rbind(real, simul)

ggplot(dataLengths, aes(Distance, fill = Data)) + 
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 10000) + 
  scale_color_brewer(palette="Dark2") + theme_minimal() + ggtitle("No uniq")
