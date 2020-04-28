ref_sp = "human"
target_sp = "mouse" 
enhancers <- c("CAGE", "ENCODE") #, "RoadMap", "GRO_seq")
col <- c("red", "navy", "forestgreen", "orange")

options(stringsAsFactors = FALSE)

path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, sep="/")
path_local <- "/home/laverre//Documents/Regulatory_Landscape/data/"

#enhancers <- c("CAGE")
conserv_vs_var <- function(var, class, class_leg){
  pdf(paste(path_evol, "/", ref_sp, "2", target_sp, "_enhancers_conserv_PhastCons_by_", var, sample, ".pdf", sep=""), width=8, height=6)
  par(mai = c(0.7, 0.7, 0.7, 0.2))
  layout(matrix(c(1,2,3,4), nrow = 2, byrow = TRUE))
  
  for (enh_data in enhancers){
    
    enh_annot <- read.table(paste(path_annot,"/", enh_data, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")
    #align <- read.table(paste(path_evol,"enhancers_conservation", enh_data, "Alignments_stats_all_species.txt", sep="/"), header=T)
    align <- read.table(paste(path_evol,"/enhancers_conservation/PhastCons/PhastCons_vertebrates_", enh_data, "_MaskedExons_Ensembl94.txt", sep=""), header=T)
    colnames(align) <- c("enh", "chr", "start", "end", "mouse", "coveredlength", "analyzedlength")
    align[which(is.na(align$mouse)),]$mouse <- 0
    
    obs_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh_data, "/original_stats.txt", sep=""), header=T)
    simul_stats <- read.table(paste(path_evol,"/enhancers_conservation/", enh_data, "/simulated_stats.txt", sep=""), header=T)
    if (sample != ""){
      obs_stats <- obs_stats[which(obs_stats[[sample]] == 1),]
      simul_stats <- simul_stats[which(simul_stats[[sample]] == 1),]
    }
    
    # Take only enhancers with unique BLAT match
    enh_annot <- enh_annot[which(enh_annot$Nb_match == 1),] # & (enh_annot$nb_N/enh_annot$length) < 0.2),]
    
    # Select aligned enhancers in specific dataset and with stats informations
    align_obs <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% obs_stats$enh),]
    align_simul <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% simul_stats$enh),]
    obs_stats <- obs_stats[which(obs_stats$enh %in% align_obs$enh),]
    simul_stats <- simul_stats[which(simul_stats$enh %in% align_simul$enh),]

    # Order enhancers
    align_obs <- align_obs[order(align_obs$enh),]
    align_simul <- align_simul[order(align_simul$enh),]
    obs_stats <- obs_stats[order(obs_stats$enh),]
    simul_stats <- simul_stats[order(simul_stats$enh),]
    

    if (var == "GC_rate"){ # Specific because GC rate information is in enh_annot file
      # Define class according to variable
      enh_annot[[var]] <- enh_annot$nb_GC/(enh_annot$length-enh_annot$nb_N)
      enh_annot$class_level <- cut(enh_annot[[var]], breaks=c(class, max(enh_annot[[var]])), include.lowest = T)
      enh_annot <- enh_annot[order(enh_annot[,1]),]
      
      # Original and simulated data in same plot
      boxplot(align_obs$mouse~enh_annot[enh_annot[,1] %in% align_obs$enh,]$class_level,
              boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), names=class_leg,
              main=paste(ref_sp, "to", target_sp, enh_data, "\n sequence conserved", sep=" "),
              ylab="Alignment score", xlab=var, ylim=c(0,0.4))
      
      boxplot(align_simul$mouse~enh_annot[enh_annot[,1] %in% align_simul$enh,]$class_level, add = TRUE, xaxt = "n", yaxt='n',
              border="firebrick3", outline=F, notch=T, boxwex=0.25, at = 1:length(class_leg) - 0.2, outcex=0.2)
      
      boxplot(align_obs$mouse~enh_annot[enh_annot[,1] %in% align_obs$enh,]$class_level, add = TRUE, xaxt = "n", yaxt='n',
              border="darkgreen", outline=F, notch=T, boxwex=0.25, at = 1:length(class_leg) + 0.2, outcex=0.2)
    
      
      }else if (var == "enh_category"){ # Specific because its qualitative variable in enh_class file
        
        enh_class <- read.table(paste(path_local,"/", ref_sp, "/", enh_data, "_classification.txt", sep=""), header=T, sep="\t")
        enh_class_obs <- enh_class[order(which(enh_class$enh %in% align_obs$enh)),]
        enh_class_simul <- enh_class[order(which(enh_class$enh %in% align_simul$enh)),]

        # Original and simulated data in same plot
        boxplot(align_obs$mouse~as.factor(enh_class_obs$type),
                boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), #names=class_leg,
                main=paste(ref_sp, "to", target_sp, enh_data, "\n sequence conserved", sep=" "),
                ylab="Alignment score", xlab=var, ylim=c(0,0.4))
        
        boxplot(align_simul$mouse~as.factor(enh_class_simul$type), add = TRUE, xaxt = "n", yaxt='n',
                border="firebrick3", outline=F, notch=T, boxwex=0.25, at = 1:length(class_leg) - 0.2, outcex=0.2)
        
        boxplot(align_obs$mouse~as.factor(enh_class_obs$type), add = TRUE, xaxt = "n", yaxt='n',
                border="darkgreen", outline=F, notch=T, boxwex=0.25, at = 1:length(class_leg) + 0.2, outcex=0.2)
        
      }else if (var == "pb_exonic" || var == "pb_genic"){
        enh_class <- read.table(paste(path_local,"/", ref_sp, "/", enh_data, "_classification.txt", sep=""), header=T, sep="\t")
        enh_class$window <- (enh_class$end + 100000) - (enh_class$start - 100000)
        
        if (var == "pb_exonic"){enh_class$prop_exonic_bp <- enh_class$X100kb_exonic_bp/enh_class$window
        }else if (var == "pb_genic"){enh_class$prop_exonic_bp <- enh_class$X100kb_genic_bp/enh_class$window}
        
        # Define class according to variable
        enh_class$class_level <-cut(enh_class$prop_exonic_bp, breaks=c(class, max(enh_class$prop_exonic_bp)), include.lowest = T)
        enh_class_obs <- enh_class[which(enh_class$enh %in% align_obs$enh),]
        enh_class_simul <- enh_class[which(enh_class$enh %in% align_simul$enh),]
        
        # Order enhancers
        enh_class_obs <- enh_class_obs[order(enh_class_obs$enh),]
        enh_class_simul <- enh_class_simul[order(enh_class_simul$enh),]
        
        # Original and simulated data in same plot
        boxplot(align_obs$mouse~enh_class_obs$class_level,
                boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), names=class_leg,
                main=paste(ref_sp, "to", target_sp, enh_data, "\n sequence conserved", sep=" "),
                ylab="Alignment score", xlab=var, ylim=c(0,0.4))
        
        boxplot(align_simul$mouse~enh_class_simul$class_level, add = TRUE, xaxt = "n", yaxt='n',
                  border="firebrick3", outline=F, notch=T, boxwex=0.25, at = 1:length(class_leg) - 0.2, outcex=0.2)
        
        boxplot(align_obs$mouse~enh_class_obs$class_level, add = TRUE, xaxt = "n", yaxt='n',
                border="darkgreen", outline=F, notch=T, boxwex=0.25, at = 1:length(class_leg) + 0.2, outcex=0.2)
        
      }else{
        # Define class according to variable
        obs_stats$class_level <-cut(obs_stats[[var]], breaks=c(class, max(obs_stats[[var]])), include.lowest = T)
        
        # Original and simulated data in same plot
        boxplot(align_obs$mouse~obs_stats$class_level,
                boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), names=class_leg,
                main=paste(ref_sp, "to", target_sp, enh_data, "\n sequence conserved", sep=" "),
                ylab="Alignment score", xlab=var, ylim=c(0,0.4))
        
        if (var != "med_score"){
          simul_stats$class_level <-cut(simul_stats[[var]], breaks=c(class, max(obs_stats[[var]])), include.lowest = T)
          
          boxplot(align_simul$mouse~simul_stats$class_level, add = TRUE, xaxt = "n", yaxt='n',
                  border="firebrick3", outline=F, notch=T, boxwex=0.25, at = 1:length(class_leg) - 0.2, outcex=0.2)}
        
        boxplot(align_obs$mouse~obs_stats$class_level, add = TRUE, xaxt = "n", yaxt='n',
                border="darkgreen", outline=F, notch=T, boxwex=0.25, at = 1:length(class_leg) + 0.2, outcex=0.2)
        
      }

  }
  
  dev.off()
  
}

sample = ""
################################################# Conserv ~ Nb genes   ################################################# 
breaks = c(1, 2, 5, 10, 20)
breaks_names = c("1-2","3-5","6-10","11-20", ">20")
conserv_vs_var("nb_gene", breaks, breaks_names)

################################################# Conserv ~ Nb baits   ################################################# 
breaks = c(1, 2, 4, 6)
breaks_names <- c("1-2","2-4","4-6",">6")
conserv_vs_var("nb_bait", breaks, breaks_names)

################################################# Conserv ~ Nb samples   ################################################# 
breaks = c(1, 2, 5, 10, 20)
breaks_names = c("1-2","3-5","6-10","11-20", ">20")
conserv_vs_var("nb_sample", breaks, breaks_names)

################################################# Conserv ~ CHiCAGO score   ################################################# 
breaks = c(5, 6, 10, 15)
breaks_names = c("5-6","6-10","10-15",">15")
conserv_vs_var("med_score", breaks, breaks_names)

################################################# Conserv ~ GC rate   ################################################# 
breaks = c(0, 0.3, 0.45, 0.55, 0.7)
breaks_names <- c("<30","30-45","45-55","55-70", ">70")
conserv_vs_var("GC_rate", breaks, breaks_names)

################################################# Conserv ~ enh category   ################################################# 
breaks = c(0, 0.3, 0.45, 0.55, 0.7)
breaks_names <- c("exonic", "genic", "intergenic")
conserv_vs_var("enh_category", breaks, breaks_names)

################################################# Conserv ~ pb exonic   ################################################# 
breaks = c(0, 0.3, 0.45, 0.55, 0.7)
breaks_names <- c("<30","30-45","45-55","55-70", ">70")
conserv_vs_var("pb_exonic", breaks, breaks_names)

################################################# Conserv ~ pb exonic   ################################################# 
breaks = c(0, 0.05, 0.1, 0.15, 0.2)
breaks_names <- c("<5%","5-10%","10-15","15-20", ">20")
conserv_vs_var("pb_genic", breaks, breaks_names)




