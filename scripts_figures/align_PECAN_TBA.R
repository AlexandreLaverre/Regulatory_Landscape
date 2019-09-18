setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/human2mouse/")
pecan <- read.table("AlignmentStatistics_pecan_Excluding_all_Exons.txt_sort", header=T)
tba <- read.table("AlignmentStatistics_TBA_Excluding_all_Exons.txt_sort", header=T)

pecan1 <- pecan
tba1 <- tba

par(mfrow=c(2,2))
######### Ungapped #######
# Total
pecan1$conserv <- pecan1$TotalUngappedLength/pecan1$human.length
hist(pecan1$conserv, main="PECAN", xlab='Total Ungapped / human length')
abline(v=median(pecan1$conserv), col='red')

tba1$conserv <- tba1$TotalUngappedLength/tba1$human.length
hist(tba1$conserv, main="TBA", xlab='Total Ungapped / human length')
abline(v=median(tba1$conserv), col='red')

# Filtered
pecan1$conserv <- pecan1$FilteredUngappedLength/pecan1$human.length
hist(pecan1$conserv, main="", xlab='Filtered Ungapped / human length')
abline(v=median(pecan1$conserv), col='red')

tba1$conserv <- tba1$FilteredUngappedLength/tba1$human.length
hist(tba1$conserv, main="", xlab='Filtered Ungapped / human length')
abline(v=median(tba1$conserv), col='red')

####### Identical ###### 
# Total
pecan1$conserv <- pecan1$TotalIdenticalLength/pecan1$human.length
hist(pecan1$conserv, main="PECAN", xlab='Total Identical / human length')
abline(v=median(pecan1$conserv), col='red')

tba1$conserv <- tba1$TotalIdenticalLength/tba1$human.length
hist(tba1$conserv, main="TBA", xlab='Total Identical / human length')
abline(v=median(tba1$conserv), col='red')

# Filtered
pecan1$conserv <- pecan1$FilteredIdenticalLength/pecan1$human.length
hist(pecan1$conserv, main="", xlab='Filtered Identical / human length')
abline(v=median(pecan1$conserv), col='red')

tba1$conserv <- tba1$FilteredIdenticalLength/tba1$human.length
hist(tba1$conserv, main="", xlab='Filtered Identical / human length')
abline(v=median(tba1$conserv), col='red')



###### Score comparison ###### 
compare <- data.frame(pecan_conserv = pecan$TotalUngappedLength/pecan$human.length)
compare$tba_conserv <- tba$TotalUngappedLength/tba$human.length
compare$frag <- tba$ID.human
compare$human.length <- tba$human.length
compare <- compare[which(compare$human.length > 250),]

par(mfrow=c(1,1))
plot(compare[which(compare$human.length > 250),]$pecan_conserv~compare[which(compare$human.length > 250),]$tba_conserv,xlab='TBA score', ylab='PECAN score', cex=0.05)
abline(0,1, col="red")
abline(lm(compare$pecan_conserv~compare$tba_conserv), col='red')


tba$score <- tba$TotalUngappedLength/tba$human.length
hist(tba[which(tba$score == 0),]$human.length, breaks=100, xlim=c(0,5000))



