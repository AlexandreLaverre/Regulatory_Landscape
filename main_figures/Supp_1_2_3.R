
cell <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/human_samples/human_MK.ibed", header=T)
list_bait <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/human/Digest_hg38_HindIII_None.txt.baitmap")
list_bait$name <- paste(list_bait$V1, list_bait$V2, list_bait$V3)

cell$bait_name <- paste(cell$bait_chr, cell$bait_start, cell$bait_end)
cell$other_name <- paste(cell$otherEnd_chr, cell$otherEnd_start, cell$otherEnd_end)
cell$other_baited <- cell$other_name %in% list_bait$name

cell$midbait <- (cell$bait_start + cell$bait_end)/2
cell$midother <- (cell$otherEnd_start + cell$otherEnd_end)/2
cell$distance <- abs(cell$midbait-cell$midother)

cell$bait_chr <- factor(cell$bait_chr, levels=levels(cell$otherEnd_chr))
cell$cis <- ifelse(cell$bait_chr == cell$otherEnd_chr, "TRUE", "FALSE")


par(mfrow=c(1,4))
# Distance distribution
cell[which(cell$distance > 3000000),]$distance <- 3000000
hist(cell[which(cell$cis == TRUE & cell$other_baited == FALSE),]$dist,breaks=40, xlim=c(0,3000000),xlab="Distance (pb)", main="Megakaryocytes (MK)")
hist(cell[which(cell$cis == TRUE & cell$other_baited == TRUE),]$dist, breaks=40, add=TRUE, col="red")

# Score distribution
cell[which(cell$score > 40),]$score <- 40
hist(cell[which(cell$cis == TRUE & cell$other_baited == FALSE),]$score,breaks=45, xlim=c(0,40), xlab="CHICAGO score", main="")
hist(cell[which(cell$cis == TRUE & cell$other_baited == TRUE),]$score, breaks=45, xlim=c(0,40), add=TRUE, col="red")

# Nb contact per bait distribution
hist(table(cell[which(cell$cis == TRUE & cell$other_baited == TRUE),]$bait_name), xlim=(c(0,80)), xlab="Contact per bait", main="", col='red')
hist(table(cell[which(cell$cis == TRUE & cell$other_baited == FALSE),]$bait_name), breaks=80, xlim=(c(0,100)), add=TRUE)

# Nb bait per contact distribution
hist(table(cell[which(cell$cis == TRUE & cell$other_baited == FALSE),]$other_name), breaks=50, xlim=c(0,15), xlab="Bait per contact", main="")
hist(table(cell[which(cell$cis == TRUE & cell$other_baited == TRUE),]$other_name), breaks=80,add=TRUE, col='red')
legend("topright", legend=c("Bait-other", "Bait-bait"), fill=c("white", "red"), bty="n")


# N reads distribution
cell[which(cell$N_reads > 300),]$N_reads <- 300
hist(cell[which(cell$cis == TRUE & cell$other_baited == FALSE),]$N_reads, breaks=100, xlim=c(0,300), xlab="reads number")
hist(cell[which(cell$cis == TRUE & cell$other_baited == TRUE),]$N_reads, breaks=100, xlim=c(0,300), add=TRUE, col="red")
