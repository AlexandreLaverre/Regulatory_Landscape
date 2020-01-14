setwd("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/")


contact_human <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions_chr_merged.txt_cell_names", header=T)
contact_mouse <- read.table("/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/all_interactions_chr_merged.txt_cell_names", header=T)
  
contact_human$bait_chr <- factor(contact_human$bait_chr, levels=levels(contact_human$chr))
contact_human$cis <- ifelse(contact_human$bait_chr == contact_human$chr, "TRUE", "FALSE")
contact_mouse$bait_chr <- factor(contact_mouse$bait_chr, levels=levels(contact_mouse$chr))
contact_mouse$cis <- ifelse(contact_mouse$bait_chr == contact_mouse$chr, "TRUE", "FALSE")

contact_human$bait_name <- paste(contact_human$bait_chr, contact_human$bait_start, contact_human$bait_end)
contact_human$other_name <- paste(contact_human$chr, contact_human$start, contact_human$end)
contact_mouse$bait_name <- paste(contact_mouse$bait_chr, contact_mouse$bait_start, contact_mouse$bait_end)
contact_mouse$other_name <- paste(contact_mouse$chr, contact_mouse$start, contact_mouse$end)


contact_human <- contact_human[which(contact_human$dist < 10000000 & contact_human$dist > 25000),]
contact_mouse <- contact_mouse[which(contact_mouse$dist < 10000000 & contact_mouse$dist > 25000),]

CEX = 1.2
LWD = 1.5
  
png("Fig1.png",width = 800, height = 600)
par(mfrow=c(2,2))

# Distribution of genomic distances -- Density plot
x <- density(contact_mouse[which(contact_mouse$cis == TRUE & contact_mouse$baited_frag == "baited"),]$dist)
plot(x ,xlab="Distance (pb)", main="",  col="orange", lwd=LWD, xlim=c(0,3000000), cex.axis=CEX, cex.lab=CEX)
lines(density(contact_mouse[which(contact_mouse$cis == TRUE & contact_mouse$baited_frag == "unbaited"),]$dist, bw=x$bw),
      col="red", lwd=LWD)

lines(density(contact_human[which(contact_human$cis == TRUE & contact_human$baited_frag == "baited"),]$dist, bw=x$bw),
      col="deepskyblue", lwd=LWD)

lines(density(contact_human[which(contact_human$cis == TRUE & contact_human$baited_frag == "unbaited"),]$dist, bw=x$bw),
      col="blue", lwd=LWD, xlim=c(0,5000000))

legend("topright", legend=c("Mouse bait-bait", "Mouse bait-other", "Human bait-bait", "Human bait-other"),
       fill=c("orange", "red", "deepskyblue", "blue"), bty='n', cex=CEX)

abline(v=median(contact_human$dist), col="blue")
abline(v=median(contact_mouse$dist), col="red")

# Distribution of Nb contact per bait
x <- density(table(contact_human[which(contact_human$cis == TRUE & contact_human$baited_frag == "unbaited"),]$bait_name))
plot(x, xlab="Contact per bait", main="", col="blue", xlim=c(0,50), ylim=c(0,0.15), cex.axis=CEX, cex.lab=CEX)
lines(density(table(contact_human[which(contact_human$cis == TRUE & contact_human$baited_frag == "baited"),]$bait_name), bw=x$bw),
              col='deepskyblue', lwd=LWD)
lines(density(table(contact_mouse[which(contact_mouse$cis == TRUE & contact_mouse$baited_frag == "baited"),]$bait_name), bw=x$bw),
              col="orange", lwd=LWD)
lines(density(table(contact_mouse[which(contact_mouse$cis == TRUE & contact_mouse$baited_frag == "unbaited"),]$bait_name), bw=x$bw),
              col="red", lwd=LWD)

legend("topright", legend=c("Mouse bait-bait", "Mouse bait-other", "Human bait-bait", "Human bait-other"),
       fill=c("orange", "red", "deepskyblue", "blue"), bty='n', cex=CEX)


# Distribution of Nb cell 
par(mfrow=c(1,1))
hist(contact_human[which(contact_human$cis == TRUE & contact_human$baited_frag == "unbaited"),]$nb_type, breaks=50,
     main='', xlab='Cell number', col="blue", cex.axis=CEX, cex.lab=CEX)

hist(contact_mouse[which(contact_mouse$cis == TRUE & contact_mouse$baited_frag == "unbaited"),]$nb_type, breaks=50,
     main='', xlab='Cell number', col="red", cex.axis=CEX, cex.lab=CEX)

dev.off()

