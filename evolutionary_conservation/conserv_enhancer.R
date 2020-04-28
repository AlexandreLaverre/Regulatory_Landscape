ref_sp = "human"
target_sp = "mouse" 
enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")
col <- c("red", "navy", "forestgreen", "orange")

path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, sep="/")

enh <- "CAGE"

################################################# Conserv ~ all species   #################################################
library(ape)
library(vioplot)

#pdf(paste(path_evol, "/", ref_sp, "2all_species_enhancers_conservation_identity.pdf", sep=""), width=6, height=6)
pdf(paste(path_evol, "/", ref_sp, "2all_species_enhancers_identity_diff_median_total_ungapped.pdf", sep=""), width=6, height=6)
x = 1
par(mfrow=c(2,2))
par(mai = c(0.5, 0.8, 0.8, 0.2))
for (enh in enhancers){
  enh_annot <- read.table(paste(path_annot,"/", enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
  align <- read.table(paste(path_evol,"enhancers_conservation", enh, "Alignments_stats_all_species_total_ungapped.txt", sep="/"), header=T)
  obs_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/original_stats.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)
  
  tree <- read.tree("/home/laverre/Documents/Regulatory_Landscape/data/ensembl_tree")
  tree <- keep.tip(tree, c("Mus_musculus", "Homo_sapiens", "Rattus_norvegicus", "Macaca_mulatta", "Oryctolagus_cuniculus", "Canis_lupus_familiaris", "Bos_taurus", "Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus"))

  # Filters
  align_obs <- align[which(align$enh %in% obs_stats$enh),]
  align_simul <- align[which(align$enh %in% simul_stats$enh),]
  species <- c("macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
  align <- align[, c("enh", species)]
  
  plot(sapply(species, function(x) median(align_obs[,x])-median(align_simul[,x]) ), ylim=c(-0.02, 0.03),
       main=enh, ylab="Difference total ungapped \n median(obs)-median(simul)", xlab='', xaxt="n", cex.lab=0.7)
  
  enh_annot <- enh_annot[which(enh_annot$nb_match == 1),]
  align_obs <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% obs_stats$enh),]
  align_simul <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% simul_stats$enh),]
  align <- align[, c("enh", species)]
  points(sapply(species, function(x) median(align_obs[,x])-median(align_simul[,x])), col="red")
  
  
  enh_annot <- enh_annot[which(enh_annot$nb_match == 1 & (enh_annot$nb_N/enh_annot$length) < 0.1),] # & (enh_annot$nb_N/enh_annot$length) < 0.2),]
  align_obs <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% obs_stats$enh),]
  align_simul <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% simul_stats$enh),]
  align <- align[, c("enh", species)]
  points(sapply(species, function(x) median(align_obs[,x])-median(align_simul[,x])), col="orange")
  
  abline(h=0, col="black", lty=2)
  if (enh =="ENCODE"){legend("topleft", fill=c("black","red","orange"), legend = c("no_filter", "BLAT_match=1", "repet_part < 0.1"), bty='n', cex=0.7)}
  axis(1, at=seq(1,9,1),labels=F)
  text(seq(1,9,1), par("usr")[3]-0.001, labels = species, pos = 1, xpd = TRUE, srt=20, cex=0.8)
  
  
  plot(tree, cex=0.8, y.lim=c(0.4,10.3), x.lim=c(0,1.07), show.tip.label=F)
  tiplabels(c("human", species), bg=F, frame = "n", adj = -0.15, cex=1)

  vioplot(align_obs[,species], at=c(0,3,6,9,12,15,18,21,24), wex=1.1,
          border='darkgreen', col="palegreen", colMed=c("white", rep("forestgreen",9)),
          las=2, cex.names = 0.5, names=species, ylab="Alignment score",
          main=paste(ref_sp, enh, "conservation", sep=" "))
  
  vioplot(align_simul[,species], at=c(1,4,7,10,13,16,19,22,25), add=T, wex=1.1,
          border='firebrick3', col="darksalmon", colMed=c("white", rep("firebrick3",9)),
          axes=F, yaxt='n', horizontal = F, las=1, cex.names = 0.8,
          main=paste(ref_sp, enh, "conservation", sep=" "))
  
  }

dev.off()



################################################# PhastCons  #################################################
by_sample="F"

if (by_sample == "TRUE"){#pdf(paste(path_evol, "/", ref_sp, "PhastCons_Placental_enhancers_conservation_by_sample.pdf", sep=""), width=9, height=5)
}else{
  #pdf(paste(path_evol, "/", ref_sp, "PhastCons_Placental_enhancers_conservation.pdf", sep=""), width=6, height=6)
  par(mai = c(0.7, 0.7, 0.7, 0.2))
  layout(matrix(c(1,2,3,4), nrow = 2, byrow = TRUE))}

for (enh in enhancers){
  enh_annot <- read.table(paste(path_annot,"/", enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
  align <- read.table(paste(path_evol,"/enhancers_conservation/PhastCons/30way/PhastCons_vertebrates_", enh, "_MaskedExons_Ensembl94.txt", sep=""), header=T)
  #mouse = align <- read.table(paste(path_evol,"/enhancers_conservation/PhastCons/PhastCons_placental_",enh,"_MaskedExons_Ensembl94.txt", sep=""), header=T)
  colnames(align) <- c("enh", "chr", "start", "end", "score", "coveredlength", "analyzedlength")
  obs_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/original_stats.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)
  
  # Filters
  align[which(is.na(align$score)),]$score <- 0
  align$score <- (align$score*align$coveredlength)/align$analyzedlength
  #enh_annot <- enh_annot[which(enh_annot$Nb_match == 1 & (enh_annot$nb_N/enh_annot$length) < 0.2),]
  
  if (by_sample == "TRUE"){
    # Plot empty boxplot
    samples = colnames(obs_stats)[7:length(obs_stats)]
    align$sample <- sample(samples, size=nrow(align), replace=T)
    
    boxplot(align$score~align$sample,
            boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), names=samples,
            main=enh, ylab="Alignment score", xlab="", ylim=c(0,0.3))
    x = 1
    
    for (sample in samples){
      obs_stats_sample <- obs_stats[which(obs_stats[[sample]] == 1),]
      simul_stats_sample <- simul_stats[which(simul_stats[[sample]] == 1),]
      
      align_obs <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% obs_stats_sample$enh),]
      align_obs$data <- "Original"
      align_simul <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% simul_stats_sample$enh),]
      align_simul$data <- "Simulated"
      align_sample <- rbind(align_obs, align_simul)
      
      boxstats <- boxplot(align_sample$score~align_sample$data, plot=F)
      if (boxstats$stats[3,1] > boxstats$stats[3,2]){if (boxstats$conf[1,1] > boxstats$conf[2,2]){test="obs>simul Significant"}else{test="obs>simul NS"}
      }else if (boxstats$conf[2,1] < boxstats$conf[1,2]){test="obs<simul Significant"}else{test="obs<simul NS"}
      
      if (enh == "ENCODE" || enh == "GRO_seq"){YLAB=""}else{YLAB="Score whithout exons"}
      
      boxplot(align_sample$score~align_sample$data, outline=F, notch=T, boxwex=0.2, add=TRUE,
              at = c(x-0.17,x+0.17), border=c("darkgreen", "firebrick3"), xaxt="n", yaxt="n")
      
      message(paste(enh, sample, ": OBS median=", boxstats$stats[3,1], "; SIMUL median=", boxstats$stats[3,2], test, sep=" "))
      x = x+1
    }
    
  }else{
    align_obs <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% obs_stats$enh),]
    align_obs$data <- "Original"
    align_simul <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% simul_stats$enh),]
    align_simul$data <- "Simulated"
    align <- rbind(align_obs, align_simul)
    
    boxstats <- boxplot(align$score~align$data, plot=F)
    if (boxstats$stats[3,1] > boxstats$stats[3,2]){if (boxstats$conf[1,1] > boxstats$conf[2,2]){test="obs>simul Significant"}else{test="obs>simul NS"}
    }else if (boxstats$conf[2,1] < boxstats$conf[1,2]){test="obs<simul Significant"}else{test="obs<simul NS"}

    boxplot(align$score~align$data,outline=F, notch=T, border=c("darkgreen", "firebrick3"), ylim=c(0,max(boxstats$stats[5,])),
            xlab="", ylab=YLAB, main=paste("PhasCons ", enh, "\n", test, sep=""))
    
    message(paste(enh, ": OBS median=", boxstats$stats[3,1], "; SIMUL median=", boxstats$stats[3,2], test, sep=" "))
  }
}

#dev.off()

################################################# Conserv human2mouse   #################################################
#pdf(paste(path_evol, "/", ref_sp, "2mouse_enhancers_conservation_identity.pdf", sep=""), width=6, height=6)
x = 1
par(mfrow=c(2,2))
par(mai = c(0.5, 0.8, 0.8, 0.2))
for (enh in enhancers){
  enh_annot <- read.table(paste(path_annot,"/", enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
  align <- read.table(paste(path_evol,"enhancers_conservation", enh, "Alignments_stats_all_species_identity.txt", sep="/"), header=T)
  obs_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/original_stats.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)
  
  # Filters
  #enh_annot <- enh_annot[which(enh_annot$Nb_match == 1),] # & (enh_annot$nb_N/enh_annot$length) < 0.2),]
  align <- align[, c("enh", "mouse")]

  align_obs <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% obs_stats$enh),]
  align_obs$data <- "Original"
  align_simul <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% simul_stats$enh),]
  align_simul$data <- "Simulated"

  align <- rbind(align_obs, align_simul)
  
  x <- boxplot(align$mouse~align$data, plot=F)
  if (x$stats[3,1] > x$stats[3,2]){if (x$conf[1,1] > x$conf[2,2]){test="obs>simul Significant"}else{test="obs>simul NS"}
  }else if (x$conf[2,1] < x$conf[1,2]){test="obs<simul Significant"}else{test="obs<simul NS"}
  
  if (enh == "ENCODE" || enh == "GRO_seq"){YLAB=""}else{YLAB="Identity Score \n whithout exons"}
  
  boxplot(align$mouse~align$data,outline=F, notch=T, border=c("darkgreen", "firebrick3"), ylim=c(0,max(x$stats[5,])),
          xlab="", ylab=YLAB, main=paste("human2mouse ", enh, "\n", test, sep=""))
}
#dev.off()

  
