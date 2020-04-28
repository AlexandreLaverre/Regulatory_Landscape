##### Difference between alignment score calcul #####

ref_sp = "human"
target_sp = "mouse" 
enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")
col <- c("red", "navy", "forestgreen", "orange")

path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, sep="/")

enh <- "CAGE"

################################################# Conserv ~ all species   #################################################
pdf(paste(path_evol, "/", ref_sp, "2all_species_enhancers_difference_between_scores.pdf", sep=""), width=6, height=6)

par(mfrow=c(2,2))
par(mai = c(0.5, 0.8, 0.8, 0.2))
for (enh in enhancers){
  obs_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/original_stats.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)

  # Total ungapped
  align <- read.table(paste(path_evol,"enhancers_conservation", enh, "Alignments_stats_all_species_total_ungapped.txt", sep="/"), header=T)
  align_obs <- align[which(align$enh %in% obs_stats$enh),]
  align_simul <- align[which(align$enh %in% simul_stats$enh),]
  species <- c("macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
  align <- align[, c("enh", species)]
  
  plot(sapply(species, function(x) mean(align_obs[,x])-mean(align_simul[,x]) ), ylim=c(-0.02, 0.03),
       main=enh, ylab="Difference mean(obs)-mean(simul)", xlab='', xaxt="n", cex.lab=0.7)

  
  # Total identity
  align <- read.table(paste(path_evol,"enhancers_conservation", enh, "Alignments_stats_all_species_total_identity.txt", sep="/"), header=T)
  align_obs <- align[which(align$enh %in% obs_stats$enh),]
  align_simul <- align[which(align$enh %in% simul_stats$enh),]
  species <- c("macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
  align <- align[, c("enh", species)]
  
  points(sapply(species, function(x) mean(align_obs[,x])-mean(align_simul[,x])), col="orange")
  
  
  # Filtered ungapped
  align <- read.table(paste(path_evol,"enhancers_conservation", enh, "Alignments_stats_all_species.txt", sep="/"), header=T)
  align_obs <- align[which(align$enh %in% obs_stats$enh),]
  align_simul <- align[which(align$enh %in% simul_stats$enh),]
  species <- c("macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
  align <- align[, c("enh", species)]
  
  points(sapply(species, function(x) mean(align_obs[,x])-mean(align_simul[,x])), col="green")
  
  
  # Filtered identity
  align <- read.table(paste(path_evol,"enhancers_conservation", enh, "Alignments_stats_all_species_identity.txt", sep="/"), header=T)
  align_obs <- align[which(align$enh %in% obs_stats$enh),]
  align_simul <- align[which(align$enh %in% simul_stats$enh),]
  species <- c("macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
  align <- align[, c("enh", species)]
  
  points(sapply(species, function(x) mean(align_obs[,x])-mean(align_simul[,x])), col="red")

  abline(h=0, col="black", lty=2)
  if (enh =="CAGE"){legend("bottomleft", fill=c("black","orange", "green","red"), 
                             legend = c("Total ungapped", "Total identity", "Filtered ungapped", "Filtered identity"), bty='n', cex=0.7)}
  axis(1, at=seq(1,9,1),labels=F)
  text(seq(1,9,1), par("usr")[3]-0.001, labels = species, pos = 1, xpd = TRUE, srt=20, cex=0.8)

}

dev.off()
