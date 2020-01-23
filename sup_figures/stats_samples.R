### Stats par echantillons
library(plotrix)

sp="human"

#setwd(paste("/home/laverre/Documents/Regulatory_Landscape/data/", sp, "/all_interactions/", sp, "_samples/", sep=""))
setwd(paste("/home/laverre/Data/Regulatory_landscape/result/simulations/", sp, "_samples/bait_all/", sep=""))

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
  par("oma"=c(0,0,0,0))
  
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
  abline(v=median(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$dist), col="orange")
  
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
  
  contact_unbait <- as.data.frame(table(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$bait_name))
  contact_bait <- as.data.frame(table(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$bait_name))
  abline(v=median(contact_unbait$Freq), col="red")
  abline(v=median(contact_bait$Freq), col="orange")
  
  # Distribution of Nb cell 
  nb_cell <- list(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$nb_type,cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$nb_type)
  
  if (sp=="mouse"){
    breaks=seq(0,7,1)
    args=seq(1,7,1)
    maxy=0.3 # 0.3
  }else{
    breaks=seq(0,16,1)
    args=seq(1,16,1)
    maxy=0.2 #0.2
  }
  
  multhist(nb_cell, breaks=breaks, names.arg=args, space=c(0.2,1),
           main='', xlab='Cell number', ylab="Frequency", col=c("orange", "red"), cex.axis=CEX, cex.lab=CEX, cex.names=CEX, freq=F)
  
  # Distribution of contacted fragment length
  cell$length <- cell$end-cell$start
  
  dens_bait <- density(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$length)
  dens_unbait <- density(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$length, bw=dens_bait$bw)
  
  if (max(dens_unbait$y) > max(dens_bait$y)){
    plot(dens_unbait ,xlab="Length (pb)", main="",  col="red", lwd=LWD, xlim=c(0,20000), cex.axis=CEX, cex.lab=CEX)
    lines(dens_bait, col="orange", lwd=LWD)
  }else{
    plot(dens_bait ,xlab="Length (pb)", main="",  col="orange", lwd=LWD, xlim=c(0,20000), cex.axis=CEX, cex.lab=CEX)
    lines(dens_unbait, col="red", lwd=LWD)
  }
  
  legend("topright", legend=c("Bait-bait", "Bait-other"),fill=c("orange", "red"), bty='n', cex=CEX)
  abline(v=median(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$length), col="red")
  abline(v=median(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$length), col="orange")
    
  par("oma"=c(0,0,3,0))
  #mtext(paste(unlist(strsplit(sample, "[.]"))[1]), outer=T, side=3, line=1, cex=CEX)
  #if simulation : 
  mtext(paste(unlist(strsplit(unlist(strsplit(sample, "[.]"))[1], "_"))[5]), outer=T, side=3, line=1, cex=CEX)

  dev.off()
}


