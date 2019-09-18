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
