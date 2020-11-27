### Stats par echantillons
library(plotrix)
options(stringsAsFactors = FALSE)
sp="mouse"
path <- "/home/laverre/Manuscript/"

setwd(paste(path,"SupplementaryDataset1/", sp, "/interactions_samples/", sep=""))
samples = list.files(pattern=".ibed")

#############################################################################################################################################
#################################################### Distribution of genomic distances ######################################################
if(sp=="human"){pdf_name="Figures/Sup_Figure1.pdf"; window=c(7,4); H=8
}else{pdf_name="Figures/Sup_Figure4.pdf"; window=c(4,4); H=5}

pdf(paste(path, pdf_name, sep=""), width = 8, height = H)
par(mfrow=window)
par(mai = c(0.4, 0.2, 0.2, 0.1)) #bottom, left, top and right 
par("oma"=c(0,0,0,0))

x=0
class_leg <- c("0", "0.5", "1", "1.5", "2")
for (sample in samples){
  cell <- read.table(sample, header=T)
  cell$cis <- ifelse(cell$bait_chr == cell$chr, "TRUE", "FALSE")
  
  cell$bait_name <- paste(cell$bait_chr, cell$bait_start, cell$bait_end)
  cell$other_name <- paste(cell$chr, cell$start, cell$end)
  
  cell <- cell[which(cell$dist < 10000000 & cell$dist > 25000),]
  
  dens_bait <- density(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$dist)
  dens_unbait <- density(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$dist, bw=dens_bait$bw)
  
  cell_name = unlist(strsplit(sample, "[.]"))[1]
  
  if (max(dens_unbait$y) > max(dens_bait$y)){
    plot(dens_unbait, xlab="",  col="red", xaxt = "n", yaxt="n", lwd=LWD, xlim=c(0,2000000), cex.axis=CEX, cex.lab=CEX, main = cell_name)
    lines(dens_bait, col="orange", lwd=LWD)
  }else{
    plot(dens_bait ,xlab="", col="orange", xaxt = "n", yaxt="n", lwd=LWD, xlim=c(0,2000000), cex.axis=CEX, cex.lab=CEX, main= cell_name)
    lines(dens_unbait, col="red", lwd=LWD)
  }
  
  med_other = median(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$dist)
  med_bait = median(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$dist)
  abline(v=med_other, col="red")
  abline(v=med_bait, col="orange")
  legend("topright", legend=c(paste("med =", round(med_other/1000, digit=0), "Kb", sep=" "), paste("med =", round(med_bait/1000, digit=0), "Kb", sep=" ")),
         col=c("red", "orange"), lty=1, bty='n', cex=0.9)
  
  axis(1, at=seq(1,2000001,500000), labels=F)
  mtext(class_leg, side=1, line=0.7, at=seq(1,2000001,500000), cex=0.8)
  
  if (x%%4 == 0){
    mtext("Density", side=2, line=0.5, cex=0.7)
    mtext("Distance (Mb)", side=1, line=1.5, cex=0.7)}
  
  x=x+1
  
}

plot.new()
legend("center", legend=c("Bait-other contacts", "Bait-bait contacts"), col=c("red", "orange"), lty=1, bty='n', cex=CEX)
dev.off()

#############################################################################################################################################
#################################################### Distribution of Nb contact per bait ####################################################
if(sp=="human"){pdf_name="Figures/Sup_Figure2.pdf"; window=c(7,4); H=8
}else{pdf_name="Figures/Sup_Figure5.pdf"; window=c(4,4); H=5}

pdf(paste(path, pdf_name, sep=""), width = 8, height = H)
par(mfrow=window)
par(mai = c(0.4, 0.2, 0.2, 0.1)) #bottom, left, top and right 
par("oma"=c(0,0,0,0))

x=0
for (sample in samples){
  cell <- read.table(sample, header=T)
  cell$cis <- ifelse(cell$bait_chr == cell$chr, "TRUE", "FALSE")
  
  cell$bait_name <- paste(cell$bait_chr, cell$bait_start, cell$bait_end)
  cell$other_name <- paste(cell$chr, cell$start, cell$end)
  
  cell <- cell[which(cell$dist < 10000000 & cell$dist > 25000),]
  cell_name = unlist(strsplit(sample, "[.]"))[1]

  dens_unbait <- density(table(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$bait_name))
  dens_bait <- density(table(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$bait_name), bw=dens_unbait$bw)

  if (max(dens_unbait$y) > max(dens_bait$y)){
    plot(dens_unbait, xlab="Contact per bait", xaxt="n", yaxt="n", main=cell_name, col="red", xlim=c(0,15), cex.axis=CEX, cex.lab=CEX)
    lines(dens_bait, col='orange', lwd=LWD)
  }else {
    plot(dens_bait, xlab="Contact per bait", xaxt="n", yaxt="n", main=cell_name, col="orange", xlim=c(0,15), cex.axis=CEX, cex.lab=CEX)
    lines(dens_unbait, col='red', lwd=LWD)
  }
  
  contact_unbait <- as.data.frame(table(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$bait_name))
  contact_bait <- as.data.frame(table(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$bait_name))
  
  med_other = median(contact_unbait$Freq)
  med_bait = median(contact_bait$Freq)
  abline(v=med_other, col="red")
  abline(v=med_bait, col="orange")
  
  legend("topright", legend=c(paste("med =", med_other, sep=" "), paste("med =", med_bait, sep=" ")),
         col=c("red", "orange"), lty=1, bty='n', cex=0.9)
  
  axis(1, at=seq(0,15,5), labels=F)
  mtext(seq(0,15,5), side=1, line=0.7, at=seq(0,15,5), cex=0.8)
  
  if (x%%4 == 0){
    mtext("Density", side=2, line=0.5, cex=0.7)
    mtext("Contact per bait", side=1, line=1.5, cex=0.7)}
  x=x+1
  
}

plot.new()
legend("center", legend=c("Bait-other contacts", "Bait-bait contacts"), col=c("red", "orange"), lty=1, bty='n', cex=CEX)
dev.off()

#############################################################################################################################################
############################################### Distribution of Nb baits per contacted fragments ############################################
if(sp=="human"){pdf_name="Figures/Sup_Figure3.pdf"; window=c(7,4); H=8
}else{pdf_name="Figures/Sup_Figure6.pdf"; window=c(4,4); H=5}

pdf(paste(path, pdf_name, sep=""), width = 8, height = H)
par(mfrow=window)
par(mai = c(0.3, 0.5, 0.2, 0.1)) #bottom, left, top and right 
par("oma"=c(0,0,0,0))

x=0
for (sample in samples){
  cell <- read.table(sample, header=T)
  cell$cis <- ifelse(cell$bait_chr == cell$chr, "TRUE", "FALSE")
  
  cell$bait_name <- paste(cell$bait_chr, cell$bait_start, cell$bait_end)
  cell$other_name <- paste(cell$chr, cell$start, cell$end)
  
  cell <- cell[which(cell$dist < 10000000 & cell$dist > 25000),]
  cell_name = unlist(strsplit(sample, "[.]"))[1]
  
  dens_unbait <- density(table(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$other_name))
  dens_bait <- density(table(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$other_name), bw=dens_unbait$bw)
  
  nb_bait <- list(table(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$other_name),
                  table(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$other_name))
  
  breaks=seq(0, range(nb_bait)[2], 1)
  
  multhist(nb_bait, xlim=c(0,31), names.arg=rep("",range(nb_bait)[2]), space=c(0.2,1), breaks=breaks,
           main=cell_name, xlab='Cell number', ylab="", col=c("red", "orange"), freq=F)
  
  mtext(seq(1,10,2), side=1, line=0.3, at=seq(2,33,6.4), cex=0.7)
  
  med_other = median(nb_bait[[1]])
  med_bait = median(nb_bait[[2]])
  
  legend("topright", legend=c(paste("med =", med_other, sep=" "), paste("med =", med_bait, sep=" ")),
         col=c("red", "orange"), lty=1, bty='n', cex=0.9)
  
  if (x%%4 == 0){
    mtext("Density", side=2, line=2.5, cex=0.7)
    mtext("Bait per contacted fragment", side=1, line=1, cex=0.7)
    }
  x=x+1
}

plot.new()
legend("center", legend=c("Bait-other contacts", "Bait-bait contacts"), col=c("red", "orange"), lty=1, bty='n', cex=CEX)
dev.off()


#############################################################################################################################################
########################################################## Distribution of Nb cell  #########################################################
if(sp=="human"){pdf_name="Figures/Sup_FigureX.pdf"; window=c(7,4); H=8
}else{pdf_name="Figures/Sup_FigureY.pdf"; window=c(4,4); H=5}

path_local <- paste("/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset1_original_interactions/", sp, "/interactions_samples/", sep="")

pdf(paste(path, pdf_name, sep=""), width = 8, height = H)
par(mfrow=window)
par(mai = c(0.3, 0.5, 0.2, 0.1)) #bottom, left, top and right 
par("oma"=c(0,0,0,0))

x=0
for (sample in samples){
  cell <- read.table(paste(path_local, sample, sep="/"), header=T)
  
  cell$cis <- ifelse(cell$bait_chr == cell$chr, "TRUE", "FALSE")
  
  cell$bait_name <- paste(cell$bait_chr, cell$bait_start, cell$bait_end)
  cell$other_name <- paste(cell$chr, cell$start, cell$end)
  
  cell <- cell[which(cell$dist < 10000000 & cell$dist > 25000),]
  cell_name = unlist(strsplit(sample, "[.]"))[1]
  
  nb_cell <- list(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$nb_type,
                  cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$nb_type)
  
  breaks=seq(0, range(nb_cell)[2], 1)
  
  multhist(nb_cell, breaks=breaks, names.arg=rep("",range(nb_cell)[2]), space=c(0.2,1),
                main=cell_name, xlab='', ylab="", col=c("orange", "red"), freq=F)
  
  mtext(seq(1,range(nb_cell)[2],2), side=1, line=0.3, at=seq(2,range(nb_cell)[2]*3,6.4), cex=0.7)
  
  if (x%%4 == 0){
    mtext("Density", side=2, line=2.5, cex=0.7)
    mtext("Cell type number", side=1, line=1, cex=0.7)
  }
  x=x+1
}

plot.new()
legend("center", legend=c("Bait-other contacts", "Bait-bait contacts"), fill=c("red", "orange"), bty='n', cex=CEX)
dev.off()


# # Distribution of contacted fragment length
# cell$length <- cell$end-cell$start
# 
# dens_bait <- density(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$length)
# dens_unbait <- density(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$length, bw=dens_bait$bw)
# 
# if (max(dens_unbait$y) > max(dens_bait$y)){
#   plot(dens_unbait ,xlab="Length (pb)", main="",  col="red", lwd=LWD, xlim=c(0,20000), cex.axis=CEX, cex.lab=CEX)
#   lines(dens_bait, col="orange", lwd=LWD)
# }else{
#   plot(dens_bait ,xlab="Length (pb)", main="",  col="orange", lwd=LWD, xlim=c(0,20000), cex.axis=CEX, cex.lab=CEX)
#   lines(dens_unbait, col="red", lwd=LWD)
# }
# 
# legend("topright", legend=c("Bait-bait", "Bait-other"),fill=c("orange", "red"), bty='n', cex=CEX)
# abline(v=median(cell[which(cell$cis == TRUE & cell$baited_frag == "unbaited"),]$length), col="red")
# abline(v=median(cell[which(cell$cis == TRUE & cell$baited_frag == "baited"),]$length), col="orange")
#   
# par("oma"=c(0,0,3,0))
# #mtext(paste(unlist(strsplit(sample, "[.]"))[1]), outer=T, side=3, line=1, cex=CEX)
# #if simulation : 
# mtext(paste(unlist(strsplit(unlist(strsplit(sample, "[.]"))[1], "_"))[5]), outer=T, side=3, line=1, cex=CEX)


