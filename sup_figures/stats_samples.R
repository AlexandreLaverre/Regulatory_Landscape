### Stats par echantillons
library(plotrix)

sp="human"

#setwd(paste("/home/laverre/Documents/Regulatory_Landscape/data/", sp, "/all_interactions/", sp, "_samples/", sep=""))
setwd(paste("/beegfs/data/alaverre/Regulatory_landscape/result/simulations/", sp, "_samples/", sep=""))

samples = list.files(pattern="_infos")

for (sample in samples){
  cell <- read.table(sample, header=T)
  
  cell$bait_chr <- factor(cell$bait_chr, levels=levels(cell$chr))
  cell$cis <- ifelse(cell$bait_chr == cell$chr, "TRUE", "FALSE")
  
  cell$bait_name <- paste(cell$bait_chr, cell$bait_start, cell$bait_end)
  cell$other_name <- paste(cell$chr, cell$start, cell$end)
  
  cell <- cell[which(cell$dist < 10000000 & cell$dist > 25000),]
  
  CEX = 1.2
  LWD = 1.5
  
  png(paste(unlist(strsplit(sample, "[.]"))[1], "_stats.png", sep=""),width = 800, height = 600)
  par(mfrow=c(2,2))
  
  # Distribution of genomic distances -- Density plot
  dens_bait <- density(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$dist)
  dens_unbait <- density(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$dist, bw=dens_bait$bw)
  
  if (max(dens_unbait$y) > max(dens_bait$y)){
    plot(dens_unbait ,xlab="Distance (pb)", main="",  col="red", lwd=LWD, xlim=c(0,1000000), cex.axis=CEX, cex.lab=CEX)
    lines(dens_bait, col="orange", lwd=LWD)
  }else{
    plot(dens_bait ,xlab="Distance (pb)", main="",  col="orange", lwd=LWD, xlim=c(0,1000000), cex.axis=CEX, cex.lab=CEX)
    lines(dens_unbait, col="red", lwd=LWD)
  }
  
  legend("topright", legend=c("Bait-bait", "Bait-other"),fill=c("orange", "red"), bty='n', cex=CEX)
  abline(v=median(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$dist), col="red")
  
  # Distribution of Nb contact per bait
  dens_unbait <- density(table(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$bait_name))
  dens_bait <- density(table(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$bait_name), bw=dens_unbait$bw)

  if (max(dens_unbait$y) > max(dens_bait$y)){
    plot(dens_unbait, xlab="Contact per bait", main="", col="red", xlim=c(0,15), cex.axis=CEX, cex.lab=CEX)
    lines(dens_bait, col='orange', lwd=LWD)
  }else {
    plot(dens_bait, xlab="Contact per bait", main="", col="orange", xlim=c(0,15), cex.axis=CEX, cex.lab=CEX)
    lines(dens_unbait, col='red', lwd=LWD)
  }
  
  legend("topright", legend=c("Bait-bait", "Bait-other"),fill=c("orange", "red"), bty='n', cex=CEX)
  
  contact_bait <- as.data.frame(table(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$bait_name))
  abline(v=median(contact_bait$Freq), col="red")
  
  # Distribution of Nb cell 
  l <- list(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$nb_type,cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$nb_type)
  
  if (sp=="mouse"){
    breaks=seq(0,7,1)
    args=seq(1,7,1)
    maxy=0.3
  }else{
    breaks=seq(0,16,1)
    args=seq(1,16,1)
    maxy=0.2
  }
  
  multhist(l, breaks=breaks, names.arg=args, ylim=c(0,maxy), space=c(0.2,1),
           main='', xlab='Cell number', ylab="Frequency", col=c("orange", "red"), cex.axis=CEX, cex.lab=CEX, cex.names=CEX, freq=F)
  
  plot.new()
  legend("center", "center", paste(unlist(strsplit(sample, "[.]"))[1]), cex=2, bty="n")
  dev.off()
}


