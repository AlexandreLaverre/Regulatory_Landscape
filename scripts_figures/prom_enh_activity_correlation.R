####################################################################################################################
##################################### Corrélation activité peaks/enhancers #########################################
####################################################################################################################
setwd("/home/laverre/Documents/Regulatory_Landscape/result/expression_correlation")
cor_mouse <- read.table("expression_correlation.txt_new", header=T)
cor_human <- read.table("expression_correlation.txt_new_simul", header=T)
#cor_mouse_valid <- read.table("expression_correlation_mouse_valid.txt", header=T)

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


## Max cor
max_human <- read.table("expression_correlation.txt_new_maxcor", header=F)
hist(max_human$V2, breaks=50)
max_human2 <- read.table("expression_correlation.txt_new_simul_maxcor", header=F)
hist(max_human2$V2, breaks=50)
