
setwd("/home/laverre/Documents/Regulatory_Landscape/result")
data_nonull <- read.table("AlignmentStatistics_TBA_mouse2human_withoutnull.txt", header=F)
data <- read.table("AlignmentStatistics_TBA_mouse2human.txt", header=T)
baited <- read.table("Align_baited_fragments.txt", header=F)
contacted <- read.table("Align_contacted_fragments.txt2", header=F)


hist(data$LengthIdentical/data$TotalLength.mouse, main="Alignement score", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = mean(data$LengthIdentical/data$TotalLength.mouse), col="red", lwd=3, lty=2)
legend("topright", legend=c("Mean = 0.268"), bty="n")


par(mfrow=c(2,2))
hist(data_nonull$V6/data_nonull$V3, , main="Whitout no_align", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = mean(data_nonull$V6/data_nonull$V3), col="red", lwd=3, lty=2)
legend("topright", legend=c("Mean = 0.356"), bty="n")

hist(baited$V6/baited$V3, , main="Baited fragments", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = mean(baited$V6/baited$V3), col="red", lwd=3, lty=2)
legend("topright", legend=c("Mean = 0.375"), bty="n")

hist(contacted$V6/contacted$V3, , main="Contacted fragments", xlab="Séq. Human aligned on mouse fragment (%)")
abline(v = mean(contacted$V6/contacted$V3), col="red", lwd=3, lty=2)
legend("topright", legend=c("Mean = 0.368"), bty="n")

par(mfrow=c(1,1))
boxplot(data$LengthIdentical/data$TotalLength.mouse, contacted$V6/contacted$V3, baited$V6/baited$V3, 
        notch=T, outline=F, names=c("All frag", "Contacted frag", "Baited frag"), ylab="Alignment score",
        main="Mouse fragments aligned with lifted human seq.")

library(ggplot2)
par(mfrow=c(1,2))
setwd("/home/laverre/Documents/Regulatory_Landscape/data/mouse")
dist_real <- read.table("dist_real.txt")
dist_simul <- read.table("dist_simul.txt")
cumul_real <- ecdf(dist_real$V1)
cumul_siml <- ecdf(dist_simul$V1)

par(mfrow=c(1,2))
boxplot(dist_real$V1, dist_simul$V1, main="Distribution",
        notch=T, outline=F, border=c("black", "red"), names=c("Real", "Simul"))
plot(cumul_real, main="Cumulative Distribution Function", xlab="Distance (pb)")
lines(cumul_siml, col="red")
legend("bottomright", legend=c( "Wasserstein \n = 312354.403\n\n\n\n"), bty="n")

hist(abs(dist_real$V1), abs(dist_simul$V1), main="Distribution")
abline(col="red", v=median())





hist(dist_real$V1, breaks=800, xlim=c(-1500000,1500000), main="Mouse real distance distribution",
     xlab="Distance (pb)")

real <- data.frame(Distance=dist_real$V1)
simul <- data.frame(Distance=dist_simul$V1)
real$Data <- 'Real'
simul$Data <- 'Simul'
dataLengths <- rbind(real, simul)

ggplot(dataLengths, aes(Distance, fill = Data)) + 
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 10000) + 
  scale_color_brewer(palette="Dark2") + theme_minimal()


geom_density(alpha = 0.2)
       
       



