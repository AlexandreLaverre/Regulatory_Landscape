setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/Sequence_conservation/")
library("RColorBrewer")

Align_score <- function(sp_origin, sp_target, data){
  dataframe <- read.table(paste(sp_origin, "/", enhancer, "/", sp_origin,"2",sp_target,data,"_merged.txt", sep=""), header=T)
  dataframe$Allexon_ungapped <- dataframe$exclude_ungapped/dataframe$all_exclude
  dataframe[which(dataframe$Allexon_ungapped == "NaN"),]$Allexon_ungapped <- 0
  
  return(dataframe$Allexon_ungapped)
}

#######  E - Enhancer sequence conservation  #######

sp_origin = "human"
type = "bait"

LWD = 1.8
CEX = 1.8
png(paste(sp_origin, "_conserv_", type, "_median_simul.png",sep=""), width = 800, height = 800)
par(mfrow=c(3,2), mai = c(0.5, 0.7, 0.5, 0.2)) # bottom, left, top, right

layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), heights = c(0.40,0.40,0.2))

enh = c("CAGE", "ENCODE", "RoadMap", "GRO_seq")

for (enhancer in enh){
  species <- c("chicken", "opossum", "rat", "mouse", "rabbit", "dog", "elephant", "cow", "macaque")
  enh_cons <- read.table(paste(sp_origin, "/", enhancer, "/", sp_origin,"2rabbit_simul_merged.txt", sep=""), header=T)
  enh_cons <- enh_cons[,c(1,8,9,10)]
  for (sp_target in species){enh_cons[[sp_target]] <- Align_score(sp_origin, sp_target, "")}
  
  y = paste("nb_", type, sep="")
  color <- brewer.pal(n = 9, name = 'Set1')
  
  ## Mean ###
  # cut_y = cut(enh_cons[[y]], breaks=c(0, 1 , 2, 10, 20, max(enh_cons[[y]])), include.lowest=T)
  # means <- data.frame(sapply(species, function(x) tapply(enh_cons[[x]],cut_y,mean)))
  # conf_inf <- data.frame(sapply(species, function(x) tapply(enh_cons[[x]],cut_y, function(y) t.test(y)[["conf.int"]][1])))
  # conf_sup <- data.frame(sapply(species, function(x) tapply(enh_cons[[x]],cut_y, function(y) t.test(y)[["conf.int"]][2])))
  # 
  # plot(means$chicken, ylim=c(0,1), xaxt='n', ylab="Score conserv Allexon ungapped", xlab=paste('Contacted ', type, sep=""), col=color[1], type='b',
  #      main=paste(sp_origin, enhancer,"conservation (mean)", sep=" "), lwd=LWD, cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
  # for (row in 1:nrow(means)){segments(x0=row,y0=conf_inf$chicken[row],x1=row,y1=conf_sup$chicken[row], col=color[1])}
  # 
  # c = 2
  # species <- c("opossum", "rat", "mouse", "rabbit", "dog", "elephant", "cow", "macaque")
  # for (sp_target in species){
  #   points(means[[sp_target]], col=color[c], type='b', lwd=LWD)
  #   for (row in 1:nrow(means)){segments(x0=row,y0=conf_inf[[sp_target]][row],x1=row,y1=conf_sup[[sp_target]][row], col=color[c])}
  #   c = c+1}
  # 
  # axis(1, at=seq(1,5,1), labels=F)
  # text(seq(1,5,1),par("usr")[3]-0.05, levels(cut_y), xpd = TRUE, cex=CEX)
   
  ## Median ###
  cut_y=cut(enh_cons[[y]], breaks=c(0, 1 , 2, 5, 10, max(enh_cons[[y]])), include.lowest=T)
  enh_bait <- boxplot(enh_cons$chicken~cut_y, plot=F)
  plot(enh_bait$stats[3,],ylim=c(0,1), xaxt='n', ylab="", xlab=paste('Contacted ', type, sep=""), col=color[1], type='b',
       main=paste(sp_origin, enhancer,"conservation (median)", sep=" "), lwd=LWD, cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
  for (row in 1:ncol(enh_bait$conf)){segments(x0=row,y0=enh_bait$conf[1,row],x1=row,y1=enh_bait$conf[2,row], col=color[1])}

  c = 2
  species <- c("opossum", "rat", "mouse", "rabbit", "dog", "elephant", "cow", "macaque")
  for (sp_target in species){
    cut_y = cut(enh_cons[[y]], breaks=c(0, 1, 2, 5, 10, max(enh_cons[[y]])), include.lowest=T)
    enh_bait <- boxplot(enh_cons[[sp_target]]~cut_y, plot=F)
    points(enh_bait$stats[3,], col=color[c], type='b', lwd=LWD)
    for (row in 1:ncol(enh_bait$conf)){segments(x0=row,y0=enh_bait$conf[1,row],x1=row,y1=enh_bait$conf[2,row], col=color[c])}
    c = c+1
    }
  axis(1, at=seq(1,5,1), labels=F)
  text(seq(1,5,1),par("usr")[3]-0.07, levels(cut_y), xpd = TRUE,  cex=CEX)
}

species <- c("chicken", "opossum", "rat", "mouse", "rabbit", "dog", "elephant", "cow", "macaque")
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0, legend = rev(species), col=rev(color), lwd=3, cex=1.5, ncol=3)

dev.off()
