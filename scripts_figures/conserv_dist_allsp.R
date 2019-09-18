##### DATA #####
sp_origin = "human"
species <- c("mouse", "dog", "cow", "elephant", "opossum", "chicken")

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.4,0.1))

for (sp_target in species){
  setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")
  obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv.txt2", sep=""), header=T)
  simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv_simul.txt2", sep=""), header=T)
  obs$class <-cut(obs$origin_dist, breaks=seq(from=25000, to=5000000, by=50000), include.lowest = T)
  simul$class <- cut(simul$origin_dist, breaks=seq(from=25000, to=5000000, by=50000), include.lowest = T)
  class_leg <- c("50Kb", "2Mb", "4Mb","6Mb","8Mb","10Mb")
  
  #### Filtres
  # Length
  obs <- obs[which(obs$bait_length > 250 & obs$bait_length < 20000 & obs$PIR_length >250 & obs$PIR_length < 20000),]
  simul <- simul[which(simul$bait_length > 250 & simul$bait_length < 20000 & simul$PIR_length >250 & simul$PIR_length < 20000),]
  # Duplication
  obs <- obs[which(obs$bait_dupli == 0 & obs$PIR_dupli == 0),]
  simul <- simul[which(simul$bait_dupli == 0 & simul$PIR_dupli == 0),]
  
  obs$ratio <- log2((obs$target_dist+1)/(obs$origin_dist+1))
  simul$ratio <- log2((simul$target_dist+1)/(simul$origin_dist+1))
  
  obs <- obs[which(!is.na(obs$target_dist)),]
  simul <- simul[which(!is.na(simul$target_dist)),]
  
  if (sp_target == "mouse"){
    obs_int_conserv <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction_pecan.txt", sep=""), header=T)
    obs_int_conserv <- obs_int_conserv[which(!is.na(obs_int_conserv$target_dist)),]
    obs$int_conserv <- obs$origin_interaction %in% obs_int_conserv$origin_interaction
    
    simul_int_conserv <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction_pecan_simul.txt", sep=""), header=T)
    simul_int_conserv <- simul_int_conserv[which(!is.na(simul_int_conserv$target_dist)),]
    simul$int_conserv <- simul$origin_interaction %in% simul_int_conserv$origin_interaction
    }
  
  intergenic <- read.table(paste("interfrag_distance_",sp_origin,"2",sp_target,"_PECAN.txt",sep=""), header=T)
  intergenic$ratio <- log2(intergenic$target_dist/intergenic$origin_dist)
  intergenic2 <- intergenic[which(intergenic$origin_dist < 3000000 ),]
  obs2 <- obs[which(!is.na(obs$target_dist) & obs$origin_dist < 3000000 ),]
  simul2 <- simul[which(!is.na(simul$target_dist) & simul$origin_dist < 3000000),]
  obs2$dist_class <-cut(obs2$origin_dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  simul2$dist_class <-cut(simul2$origin_dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  intergenic2$dist_class <-cut(intergenic2$origin_dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  class_leg <- c("50Kb", "1Mb", "2Mb","3Mb", "4Mb", "5Mb")
  
  # Overlap > 1 potential enh
  setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/")
  obs2_enh <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,".txt",sep=""), header=T)
  simul2_enh <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,"_simul.txt",sep=""), header=T)
  
  obs2$PIR <- as.factor(sub(".*-", "", obs2$origin_interaction))
  simul2$PIR <- as.factor(sub(".*-", "", simul2$origin_interaction))
  
  obs2$CAGE <- obs2$PIR %in%  obs2_enh[which(obs2_enh$CAGE_count > 0),]$PIR
  obs2$ENCODE <- obs2$PIR %in%  obs2_enh[which(obs2_enh$ENCODE_count > 0),]$PIR
  obs2$GRO_seq <- obs2$PIR %in%  obs2_enh[which(obs2_enh$GRO_seq_count > 0),]$PIR
  obs2$RoadMap <- obs2$PIR %in%  obs2_enh[which(obs2_enh$RoadMap_count > 0),]$PIR
  simul2$CAGE <- simul2$PIR %in%  simul2_enh[which(simul2_enh$CAGE_count > 0),]$PIR
  simul2$ENCODE <- simul2$PIR %in%  simul2_enh[which(simul2_enh$ENCODE_count > 0),]$PIR
  simul2$GRO_seq <- simul2$PIR %in%  simul2_enh[which(simul2_enh$GRO_seq_count > 0),]$PIR
  simul2$RoadMap <- simul2$PIR %in%  simul2_enh[which(simul2_enh$RoadMap_count > 0),]$PIR
  
  obs2_box_enh <-  boxplot(obs2[which(obs2$CAGE == TRUE),]$ratio~obs2[which(obs2$CAGE == TRUE),]$dist_class, plot=F)
  obs2_box_cons <- boxplot(obs2[which(obs2$int_conserv == TRUE),]$ratio~obs2[which(obs2$int_conserv == TRUE),]$dist_class, plot=F)
  obs2_box <- boxplot(obs2$ratio~obs2$dist_class, plot=F)
  simul2_box <- boxplot(simul2$ratio~simul2$dist_class, plot=F)
  intergenic2_box <- boxplot(intergenic2$ratio~intergenic2$dist_class, plot=F)
  
  # ylim=c(-0.3,0)
  # ylim=c(-0.1,0.35)
  if (sp_target == "mouse"){
    par(mar = c(2,4.2,3,1))
    plot(obs2_box_cons$stats[3,], type="b", col="orange", main=paste(sp_origin, 'to', sp_target), cex=0.7, 
         ylim=c(min(obs2_box$conf[1,],simul2_box$conf[1,], obs2_box_cons$stats[3,]), max(obs2_box$conf[2,],simul2_box$conf[2,], obs2_box_cons$stats[3,])),
         xlab="", ylab=paste("log2(dist", sp_target, "/ dist", sp_origin, ")"), xaxt='n', cex.lab=1.2)
    }else{
    par(mar = c(2,4.2,3,1))
    plot(obs2_box_cons$stats[3,], type="b", col="orange", main=paste(sp_origin, 'to', sp_target), cex=0.7, 
         ylim=c(min(obs2_box$conf[1,],simul2_box$conf[1,]), max(obs2_box$conf[2,],simul2_box$conf[2,])),
         xlab="", ylab=paste("log2(dist", sp_target, "/ dist", sp_origin, ")"), xaxt='n', cex.lab=1.2)
    }
  for (row in 1:ncol(obs2_box_cons$stats)){segments(x0=row,y0=obs2_box_cons$conf[1,row],x1=row,y1=obs2_box_cons$conf[2,row], lwd=0.3, col='orange')}
  
  points(obs2_box$stats[3,], col="red", type="b", cex=0.7)
  for (row in 1:ncol(obs2_box$stats)){segments(x0=row,y0=obs2_box$conf[1,row],x1=row,y1=obs2_box$conf[2,row], lwd=0.3, col='red')}
  
  points(simul2_box$stats[3,], col="blue", type="b", cex=0.7)
  for (row in 1:ncol(simul2_box$stats)){segments(x0=row,y0=simul2_box$conf[1,row],x1=row,y1=simul2_box$conf[2,row], lwd=0.3, col='blue')}
  
  points(intergenic2_box$stats[3,], col="black", type="b", cex=0.7)
  for (row in 1:ncol(intergenic2_box$stats)){segments(x0=row,y0=intergenic2_box$conf[1,row],x1=row,y1=intergenic2_box$conf[2,row], lwd=0.3, col='black')}
  
  points(obs2_box_enh$stats[3,], col="forestgreen", type="b", cex=0.7)
  for (row in 1:ncol(obs2_box_enh$stats)){segments(x0=row,y0=obs2_box_enh$conf[1,row],x1=row,y1=obs2_box_enh$conf[2,row], lwd=0.3, col='forestgreen')}
  
  axis(1, at=seq(1,101,20), labels=F)
  if (sp_target == "elephant") {
    text(seq(1,101,20), par("usr")[3]-0.005, labels = class_leg, pos = 1, xpd = TRUE)
    }else{
    text(seq(1,101,20), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE)
    }
}

par(mar = c(0,0,3,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Inter-fragments", "Simulated", "Observed", "Observed with enh", "Conserved interaction"), 
       col=c("black", "blue", "red", "forestgreen", "orange"), lwd=5, horiz = TRUE, cex = 1.2)

