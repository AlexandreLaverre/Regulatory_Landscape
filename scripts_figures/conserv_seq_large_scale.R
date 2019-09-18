################################### conservation de séquence par classe de distance ###################################
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/")

sp_origin = "human"
sp_target <- c("mouse", "dog", "cow", "elephant", "opossum", "chicken")
m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.4,0.1))

sp_target = "mouse"
for (sp_target in species){
  obs <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,".txt", sep=""), header=T)
  simul <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,"_simul.txt", sep=""), header=T)

  #### Filtres ###
  # Length
  obs <- obs[which(obs$length > 250 & obs$length < 20000),]
  simul <- simul[which(simul$length > 250 & simul$length < 20000),]
  # Duplication
  obs <- obs[which(obs$Duplication == 0),]
  simul <- simul[which(simul$Duplication == 0),]

  obs <- obs[which(obs$midist_obs < 3000000),]
  simul <- simul[which(simul$midist_obs < 3000000),]
  obs$class <-cut(obs$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  simul$class <-cut(simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  class_leg <- c("50Kb", "500Kb", "1Mb", "1.5Mb", "2Mb","2.5Mb")
  
  ##### Score obs par classe de distance #### 
  # Extract all exons
  obs_dist <- data.frame(inter = sapply(levels(obs$class), function(x) mean(obs[which(obs$class == x),]$Allexon_ungapped)))
  obs_dist$int_start <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$Allexon_ungapped)[["conf.int"]][1])
  obs_dist$int_end <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$Allexon_ungapped)[["conf.int"]][2])
  
  obs_dist_genes <- data.frame(inter = sapply(levels(obs$class), function(x) mean(obs[which(obs$class == x),]$Allgene_ungapped)))
  obs_dist_genes$int_start <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$Allgene_ungapped)[["conf.int"]][1])
  obs_dist_genes$int_end <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$Allgene_ungapped)[["conf.int"]][2])
  
  enh_dist <- data.frame(inter = sapply(levels(obs$class), function(x) mean(obs[which(obs$class == x & obs$CAGE_count > 0),]$Allexon_ungapped)))
  enh_dist$int_start <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x & obs$CAGE_count > 0),]$Allexon_ungapped)[["conf.int"]][1])
  enh_dist$int_end <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x & obs$CAGE_count > 0),]$Allexon_ungapped)[["conf.int"]][2])

  simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$Allexon_ungapped)))
  simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$Allexon_ungapped)[["conf.int"]][1])
  simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$Allexon_ungapped)[["conf.int"]][2])
  
  if ((sp_target == "mouse") || (sp_target == "elephant")) {
    par(mar = c(2,4.2,3,0.7))
    plot(obs_dist$inter[1:50], type="l", col="red", cex=0.7, 
         ylim=c(min(simul_dist$inter, enh_dist$inter, obs_dist$inter),max(simul_dist$inter, enh_dist$inter, obs_dist$inter)),
         main=paste(sp_origin,"to",sp_target, sep=" "), xlab="", ylab="Mean alignment score", xaxt = "n", cex.lab=1.3)
  }else{
    par(mar = c(2,3,3,1.9))
    plot(obs_dist$inter[1:50], type="l", col="red", cex=0.7, 
         ylim=c(min(simul_dist$inter, enh_dist$inter, obs_dist$inter),max(simul_dist$inter, enh_dist$inter, obs_dist$inter)),
         main=paste(sp_origin,"to",sp_target, sep=" "),  xlab="", ylab="", xaxt = "n")
    }
  
    for (row in 1:nrow(obs_dist[1:50,])){
    segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}
  
  #points(obs_dist_genes$inter[1:50], type="l", col="orange", cex=0.7)
  #for (row in 1:nrow(enh_dist[1:50,])){
  #  segments(x0=row,y0=obs_dist_genes[row,]$int_start,x1=row,y1=obs_dist_genes[row,]$int_end, col='orange', lwd=0.3)}
  
  points(enh_dist$inter[1:50], type="l", col="forestgreen", cex=0.7)
  for (row in 1:nrow(enh_dist[1:50,])){
    segments(x0=row,y0=enh_dist[row,]$int_start,x1=row,y1=enh_dist[row,]$int_end, col='forestgreen', lwd=0.3)}
  
  points(simul_dist$inter[1:50], type="l", col="blue", cex=0.7)
  for (row in 1:nrow(simul_dist[1:50,])){
    segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
  
  axis(1, at=seq(1,51,10), labels=F)
  if (sp_target == "opossum"){
    text(seq(1,51,10),par("usr")[3]-min(obs_dist_genes$inter)/1.5, class_leg, xpd = TRUE, cex=1.2)}
  else if (sp_target == "chicken") {
    text(seq(1,51,10),par("usr")[3]-min(obs_dist_genes$inter)*3, class_leg, xpd = TRUE, cex=1.2)}
  else if (sp_target == "dog") {
    text(seq(1,51,10),par("usr")[3]-min(obs_dist_genes$inter)/10, class_leg, xpd = TRUE, cex=1.2)}
  else{text(seq(1,51,10),par("usr")[3]-min(obs_dist_genes$inter)/5, class_leg, xpd = TRUE, cex=1.2)}
}

par(mar = c(0,0,3,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Observed without genes", "Simulated non-exonic", "Observed non-exonic","Observed non-exonic with enh"), 
       col=c("orange","blue", "red", "forestgreen"), lwd=5, horiz = TRUE, cex = 1.2)

text(seq(1,51,10),par("usr")[3]-0.01, class_leg, xpd = TRUE)
legend(x = "topleft",
       legend = c("Simulated", "Observed","With enhancers"), 
       fill=c("blue", "red", "forestgreen"), cex = 1.2, bty='n')

## Boxplot effet des enhancers
m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.4,0.1))

for (sp_target in species){
  obs <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,".txt", sep=""), header=T)
  simul <- read.table(paste("PIR_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,"_simul.txt", sep=""), header=T)
  
  #### Filtres ###
  # Length
  obs <- obs[which(obs$length > 250 & obs$length < 20000),]
  simul <- simul[which(simul$length > 250 & simul$length < 20000),]
  # Duplication
  obs <- obs[which(obs$Duplication == 0),]
  simul <- simul[which(simul$Duplication == 0),]
  
  obs <- obs[which(obs$midist_obs < 3000000),]
  simul <- simul[which(simul$midist_obs < 3000000),]
  obs$class <-cut(obs$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  simul$class <-cut(simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  class_leg <- c("50Kb", "500Kb", "1Mb", "1.5Mb", "2Mb","2.5Mb")
  
  if ((sp_target == "mouse") || (sp_target == "elephant")) {
    par(mar = c(2,4.2,3,0.7))
    boxplot(simul$Allexon_ungapped, obs$Allexon_ungapped,
            obs[which(obs$GRO_seq_count > 0),]$Allexon_ungapped, obs[which(obs$RoadMap_count > 0),]$Allexon_ungapped, 
            obs[which(obs$ENCODE_count > 0),]$Allexon_ungapped, obs[which(obs$CAGE_count > 0),]$Allexon_ungapped, ylim=c(0,1.1),
            outline=F, notch=T, col=c(rgb(0,0,1, 0.7), rgb(1,0,0, 0.7), 'green1', 'green2', 'green3', 'green4'), width=proportion, boxwex=c(1,1,1,1,1,1), at=c(1,2,2.65,3.1,3.55,3.85),
            main=paste(sp_origin,"to",sp_target, sep=" "), xlab="", ylab="Mean align score", xaxt = "n", cex.lab=1.2)
    }else{
      par(mar = c(2,3,3,1.9))
      boxplot(simul[which(simul$CAGE_count == 0),]$Allexon_ungapped, obs[which(obs$CAGE_count == 0 & obs$RoadMap_count == 0 & obs$ENCODE_count == 0),]$Allexon_ungapped,
              obs[which(obs$GRO_seq_count > 0),]$Allexon_ungapped, obs[which(obs$RoadMap_count > 0),]$Allexon_ungapped, 
              obs[which(obs$ENCODE_count > 0),]$Allexon_ungapped, obs[which(obs$CAGE_count > 0),]$Allexon_ungapped,
              outline=F, notch=T, col=c('blue', 'firebrick2', 'green1', 'green2', 'green3', 'green4'), 
              main=paste(sp_origin,"to",sp_target, sep=" "), xlab="", ylab="", xaxt = "n")
      }
}

par(mar = c(1,0,1,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Simulated non-exonic", "Observed non-exonic without enh","Observed non-exonic with GRO_seq",
                  "Observed non-exonic with RoadMap", "Observed non-exonic with ENCODE", "Observed non-exonic with CAGE"), 
       col=c('blue', 'firebrick2', 'green1', 'green2', 'green3', 'green4'), ncol=2, lwd=5, cex = 1.2)




proportion <- c(length(simul$Allexon_ungapped), length(obs$Allexon_ungapped),
                length(obs[which(obs$GRO_seq_count > 0),]$Allexon_ungapped), length(obs[which(obs$RoadMap_count > 0),]$Allexon_ungapped), 
                length(obs[which(obs$ENCODE_count > 0),]$Allexon_ungapped), length(obs[which(obs$CAGE_count > 0),]$Allexon_ungapped))


boxplot(simul$Allexon_ungapped, obs$Allexon_ungapped, obs[which(obs$CAGE_count > 0),]$Allexon_ungapped, ylim=c(0,1.1),
        outline=F, notch=T, col=c("dodgerblue3", "firebrick1","forestgreen"), box.wex=0.1,
        main=paste(sp_origin,"to",sp_target, sep=" "), xlab="", ylab="Mean alignment score", xaxt = "n", cex.lab=1.3)

axis(1, at=c(1,2,3), labels = c("Simulated", "Observed", "With enhancers"), cex.axis=1.2)

segments(1, 1.08, 2, 1.08) 
text("***", x= 1.5, y=1.1, cex=1.7)
segments(2, 1.04, 3, 1.04) 
text("***",x= 2.5, y=1.06, cex=1.7)


