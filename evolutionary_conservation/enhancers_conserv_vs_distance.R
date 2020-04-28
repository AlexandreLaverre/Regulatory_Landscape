options(stringsAsFactors = FALSE)
ref_sp = "mouse"
target_sp = "human" 
enhancers <- c("CAGE", "ENCODE") #, "RoadMap", "GRO_seq")
col <- c("red", "navy", "forestgreen", "orange")

path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, sep="/")
path_local <- "/home/laverre//Documents/Regulatory_Landscape/data/"

############ Function for align score #########
dist_vs_align_score <- function(stats_data, enh, enh_annot, align){
  #align <- read.table(paste(path_evol,"enhancers_conservation", enh, "Alignments_stats_all_species.txt", sep="/"), header=T)
  align <- read.table(paste(path_evol,"/enhancers_conservation/PhastCons/PhastCons_vertebrates_", enh, "_MaskedExons_Ensembl94.txt", sep=""), header=T)
  colnames(align) <- c("enh", "chr", "start", "end", "mouse", "coveredlength", "analyzedlength")
  align[which(is.na(align$mouse)),]$mouse <- 0
  
  align <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% stats_data$enh),]
  align <- align[order(align$enh),]
  
  conserv <- data.frame(result = sapply(levels(stats_data$class_dist), function(x) 
    mean(align[which(align$enh %in% stats_data[which(stats_data$class_dist == x),]$enh),]$mouse)))
  
  conserv$int_start <- sapply(levels(stats_data$class_dist), function(x)
    tryCatch(t.test(align[which(align$enh %in% stats_data[which(stats_data$class_dist == x),]$enh),]$mouse)[["conf.int"]][1],error=function(e) 0))
  conserv$int_end <- sapply(levels(stats_data$class_dist), function(x)
    tryCatch(t.test(align[which(align$enh %in% stats_data[which(stats_data$class_dist == x),]$enh),]$mouse)[["conf.int"]][2],error=function(e) 0))
  
  
  #conserv$int_start <- sapply(levels(stats_data$class_dist), function(x)
    #boxplot.stats(align[which(align$enh %in% stats_data[which(stats_data$class_dist == x),]$enh),]$mouse)[["conf"]][1])
  #conserv$int_end <- sapply(levels(stats_data$class_dist), function(x)
    #boxplot.stats(align[which(align$enh %in% stats_data[which(stats_data$class_dist == x),]$enh),]$mouse)[["conf"]][2])
  
  return(conserv)
}

############ Function for other variable results #########
dist_vs_other_var <- function(stats_data, enh, enh_annot, var){
  if (var == "GC_rate"){
    conserv <- data.frame(result = sapply(levels(stats_data$class_dist), function(x) 
      mean(enh_annot[which(enh_annot[,1] %in% stats_data[which(stats_data$class_dist == x),]$enh),]$GC_rate)))
    
    # conserv$int_start <- sapply(levels(stats_data$class_dist), function(x) 
    #   t.test(enh_annot[which(enh_annot[,1] %in% stats_data[which(stats_data$class_dist == x),]$enh),]$GC_rate)[["conf.int"]][1])
    # conserv$int_end <- sapply(levels(stats_data$class_dist), function(x) 
    #   t.test(enh_annot[which(enh_annot[,1] %in% stats_data[which(stats_data$class_dist == x),]$enh),]$GC_rate)[["conf.int"]][2])

  }else if (var == "enh_class"){
    enh_class <- read.table(paste(path_local,"/", ref_sp, "/", enh, "_classification.txt", sep=""), header=T, sep="\t")
    enh_class$window <- (enh_class$end + 500000) - (enh_class$start - 500000)
    enh_class$prop_exonic_bp <- enh_class$X500kb_exonic_bp/enh_class$window
    
    conserv <- data.frame(result = sapply(levels(stats_data$class_dist), function(x) 
      mean(enh_class[which(enh_class$enh %in% stats_data[which(stats_data$class_dist == x),]$enh),]$prop_exonic_bp)))
    
  }else{
    conserv <- data.frame(result = sapply(levels(stats_data$class_dist), function(x) 
      mean(stats_data[which(stats_data$class_dist == x),][[var]])))
    
    #conserv$int_start <- sapply(levels(stats_data$class_dist), function(x)
    #   boxplot.stats(stats_data[which(stats_data$class_dist == x),][[var]])[["conf"]][1])
    #conserv$int_end <- sapply(levels(stats_data$class_dist), function(x)
    #  boxplot.stats(stats_data[which(stats_data$class_dist == x),][[var]])[["conf"]][2])
    
    conserv$int_start <- sapply(levels(stats_data$class_dist), function(x)
       t.test(stats_data[which(stats_data$class_dist == x),][[var]])[["conf.int"]][1])
    conserv$int_end <- sapply(levels(stats_data$class_dist), function(x)
       t.test(stats_data[which(stats_data$class_dist == x),][[var]])[["conf.int"]][2])
     
  }

  return(conserv)
}


dist_vs_var <- function(var){
  #pdf(paste(path_evol, "/", ref_sp, "2", target_sp, "_", var, "_to_distance", sample, ".pdf", sep=""), width=8, height=6)
  pdf(paste(path_evol, "/", ref_sp, "_PhastCons_", var, "_to_distance", sample, ".pdf", sep=""), width=8, height=6)
  x = 1 # To change color between each enhancers dataset
  
  for (enh in enhancers){
    ############ Datas #########
    enh_annot <- read.table(paste(path_annot,"/", enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
    obs_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/original_stats.txt", sep=""), header=T)
    simul_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)
    max_dist = 2000000
    
    if (sample != ""){
      obs_stats <- obs_stats[which(obs_stats[[sample]] == 1),]
      simul_stats <- simul_stats[which(simul_stats[[sample]] == 1),]
      max_dist = 1000000
    }
    
    # Take only enhancers with unique BLAT match
    enh_annot <- enh_annot[which(enh_annot$nb_match == 1),] # & (enh_annot$nb_N/enh_annot$length) < 0.2),]
    if (var == "GC_rate"){enh_annot[[var]] <- enh_annot$nb_GC/(enh_annot$length-enh_annot$nb_N)}
    
    ### Distance classes
    obs_stats$class_dist <-cut(obs_stats$med_dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
    simul_stats$class_dist <-cut(simul_stats$med_dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
    class_leg <- c("0", "500Kb", "1Mb", "1.5Mb", "2Mb")
    
    ############ Calcul median (or mean) and confidence intervals of interested variable #########
    if (var == "align_score"){
      conserv <- dist_vs_align_score(obs_stats, enh, enh_annot, align)
      conserv_simul <- dist_vs_align_score(simul_stats, enh, enh_annot, align)
      YLIM = c(0, 0.25)
      ylabel = "Alignment score (mean)"
      
    }else{
      conserv <- dist_vs_other_var(obs_stats, enh, enh_annot, var)
      conserv_simul <- dist_vs_other_var(simul_stats, enh, enh_annot, var)
      if (var == "GC_rate"){YLIM = c(0.35, 0.55)
      ylabel="GC rate"} 
      else if (var == "nb_sample"){YLIM=c(2,12)
      ylabel="Number of sample (mean)"}
      else if (var == "enh_class"){YLIM=c(0.5,0.85)
      ylabel = "Proportion of exonic bp in window"}
      else if (var == "nb_bait"){YLIM=c(1,6)
      ylabel = "Number of contacted baits (mean)"}
      else if (var == "med_score"){YLIM=c(6,7)
      ylabel = "CHiCAGO score (mean)"}
      else {YLIM = c(0, 10)
      ylabel = "Number of contacted genes (mean)"}}
  
    ############ Plot #########
    if (enh == "CAGE"){ # First enhancers dataset
      plot(conserv$result, type="l", col=col[x], xaxt = "n", ylim=YLIM,
           xlab="", ylab=ylabel,
           main=paste(ref_sp, "to", target_sp, "enhancers conservation", sep=" "))
      }else{lines(conserv$result, type="l", col=col[x])} # Add lines of other enhancers datasets
    
    # Add simulated lines
    lines(conserv_simul$result, type="l", lty=2, col=col[x], cex=0.8)
    
    # Add confidences intervals
    for (row in 1:nrow(conserv)){
       segments(x0=row,y0=conserv[row,"int_start"],x1=row,y1=conserv[row,"int_end"], col=col[x], lwd=0.3)}
    for (row in 1:nrow(conserv_simul)){
       segments(x0=row,y0=conserv_simul[row,"int_start"],x1=row,y1=conserv_simul[row,"int_end"], col=col[x], type="l", lty=2, lwd=0.3)}
     
    x = x + 1 # To change color between each enhancers dataset
  }

  ######### Add axis and legends at the end #########
  axis(1, at=seq(1,51,10), labels=F)
  text(seq(1,51,10),par("usr")[3]-0.05, class_leg, xpd = TRUE)
  legend("topleft", fill=c("navy","red","forestgreen", "orange"), legend = c("ENCODE", "CAGE", "RoadMap", "GRO_seq"), bty='n', cex=0.8)
  dev.off()
  
}

sample = ""
################################################# Distance ~ Alignement score   ################################################# 
dist_vs_var("align_score")

################################################# Distance ~ Nb genes   ################################################# 
dist_vs_var("nb_gene")

################################################# Distance ~ Nb baits   ################################################# 
dist_vs_var("nb_bait")

################################################# Distance ~ CHICAGO score   ################################################# 
dist_vs_var("med_score")

################################################# Distance ~ GC rate   ################################################# 
dist_vs_var("GC_rate")

################################################# Distance ~ Nb samples   ################################################# 
dist_vs_var("nb_sample")

################################################# Distance ~ Enhancers classification   ################################################# 
dist_vs_var("enh_class")


pdf(paste(path_evol, "/", ref_sp, "2", target_sp, "_nb_bait.pdf", sep=""), width=6, height=6)
par(mai = c(0.4, 0.7, 0.7, 0.2))
layout(matrix(c(1,2,3,4), nrow = 2, byrow = TRUE))
for (enh in enhancers){
  obs_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/original_stats.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)
  obs <- obs_stats$nb_bait
  obs_data <- rep("Original", length(obs))
  simul <- simul_stats$nb_bait
  simul_data <- rep("Simulated", length(simul))
  test <- data.frame(result = c(obs, simul), data=c(obs_data, simul_data))
  boxplot(test$result~test$data, outline=F, notch=T, border=c("darkgreen", "firebrick3"), xlab="",
          ylab="Contacted baits number", main=enh)
}
dev.off()

  
                  