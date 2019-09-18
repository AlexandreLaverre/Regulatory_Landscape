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
