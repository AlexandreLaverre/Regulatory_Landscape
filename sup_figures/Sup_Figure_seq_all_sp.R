library(data.table)
options(stringsAsFactors = FALSE)

ref_sp = "human"
target_sp = "mouse"

path <- "/home/laverre/Data/Regulatory_landscape/result/"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, "restriction_fragments", sep="/")

####################################### Restriction fragments datas  ####################################### 
path_data <- "/home/laverre/Manuscript/"
obs <- read.table(paste(path_data, "SupplementaryDataset5/", ref_sp, "/statistics_contacted_sequence_original.txt", sep=""), header=T)
simul <- read.table(paste(path_data, "SupplementaryDataset5/", ref_sp,"/statistics_contacted_sequence_simulated.txt", sep=""), header=T)

obs <- obs[which(obs$baited == "unbaited"),]
simul <- simul[which(simul$baited == "unbaited"),]
obs <- obs[which(obs$BLAT_match < 2),]
simul <- simul[which(simul$BLAT_match < 2),]

obs$ID <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
simul$ID <-  do.call(paste,c(simul[c("chr","start","end")],sep=":"))

################# Alignment score vs all species #################
align <- read.table(paste(path_evol,"sequence_conservation/restriction_fragments/Alignments_stats_all_species_total_ungapped.txt", sep="/"), header=T)

species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")

align_obs <- align[which(align$enh %in% obs$ID), c("enh", species)]
align_obs_enh <- align[which(align$enh %in% obs[which(obs$ENCODE_bp > 0),]$ID), c("enh", species)]
align_simul <- align[which(align$enh %in% simul$ID), c("enh", species) ]

################# Sequence conservation human 2 all species vs distance to promoters #################
CEX=1.2
CEX_lines=1

obs$dist_class <- cut(obs$median_dist, breaks=seq(from=25000, to=2025000, by=50000), include.lowest = T)
simul$dist_class <-cut(simul$median_dist, breaks=seq(from=25000, to=2025000, by=50000), include.lowest = T)
obs <- obs[order(obs$ID), ]
simul <- simul[order(simul$ID), ]

path <- "/home/laverre/Manuscript/Figures/"

if(ref_sp == "human"){pdf_name="Sup_Figure_align_all_human.pdf"}else{pdf_name="Sup_Figure_align_all_mouse.pdf"}

pdf(paste(path, pdf_name, sep=""))
par(mai = c(0.4, 0.6, 0.4, 0.1)) #bottom, left, top and right 
par(mfrow=c(3,3))

for (sp in species){
  obs$align_target <- align_obs[order(align_obs$enh), sp]
  simul$align_target <- align_simul[order(align_simul$enh), sp]
  
  obs_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$align_target)))
  obs_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][1])
  obs_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][2])
  
  obs_enh_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)))
  obs_enh_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][1])
  obs_enh_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][2])
  
  simul_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$align_target)))
  simul_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][1])
  simul_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][2])
  
  class_leg <- c("0", "0.5", "1", "1.5", "2")
  
  xmin=min(simul_dist$int_start)-0.02
  xmax=max(obs_enh_dist$int_end)+0.02

  if (sp == "chicken" | sp == "opossum"){xmin=abs(xmin)+0.01; xmax=xmax-0.02}
  
  print(sp)
  plot(obs_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines, main=sp,
       xlab="", ylab="Alignment score", xaxt = "n", ylim=c(xmin,xmax), cex.lab=CEX, cex.axis=CEX, las=2)
  
  for (row in 1:nrow(obs_dist)){
    segments(x0=row,y0=obs_dist[row,"int_start"],x1=row,y1=obs_dist[row,"int_end"], col='forestgreen', lwd=0.3)}
  
  lines(simul_dist[,"inter"], type="l", col="firebrick1", cex=CEX_lines)
  for (row in 1:nrow(simul_dist)){
    segments(x0=row,y0=simul_dist[row,"int_start"],x1=row,y1=simul_dist[row,"int_end"], col='firebrick1', lwd=0.3)}
  
  lines(obs_enh_dist[,"inter"], type="l", col="dodgerblue", cex=CEX_lines)
  for (row in 1:nrow(obs_enh_dist)){
    segments(x0=row,y0=obs_enh_dist[row,"int_start"],x1=row,y1=obs_enh_dist[row,"int_end"], col='dodgerblue', lwd=0.3)}
  
  axis(1, at=seq(1,nrow(obs_dist)+1,10), labels=F)
  mtext(class_leg, side=1, line=0.7, at=seq(1,nrow(obs_dist)+1,10), cex=0.8)
  mtext("Distance to promoters (Mb)", side=1, line=2, cex=0.8)
  
}

dev.off()