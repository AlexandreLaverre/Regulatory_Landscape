options(stringsAsFactors = FALSE)
ref_sp = "human"

enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")
species <- c("mouse")

path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, sep="/")


par(mfrow=c(2,1))

col <- c("red", "navy", "forestgreen", "orange")
color_n = 1 # To change color between each enhancers dataset

for (enh in enhancers){
  synt_obs <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/human2", sp, "_", enh, "_original_synteny.txt", sep=""), header=T)
  synt_simul <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/human2", sp, "_", enh, "_simulated_synteny.txt", sep=""), header=T)
  
  # Filters
  synt_obs <- synt_obs[which(synt_obs$BLAT_match == 1 & synt_obs$align_score > 0.1),]
  synt_simul <- synt_simul[which(synt_simul$BLAT_match == 1 & synt_simul$align_score > 0.1),]
  synt_obs <- synt_obs %>% replace_with_na(replace = list(target_dist = "trans"))
  synt_simul <- synt_simul %>% replace_with_na(replace = list(target_dist = "trans"))
  
  synt_obs <- synt_obs[which(as.numeric(synt_obs$origin_dist) < 10000000),]
  synt_simul <- synt_simul[which(as.numeric(synt_simul$origin_dist) < 10000000),]
  
  x <- density(log(as.numeric(synt_obs$origin_dist)))
  if (enh == "CAGE"){
    plot(x, xlab="Distance in human (log)", main="Genes-enh interactions human to mouse",  col=col[color_n], ylim=c(0,0.5) )
  }else{
    lines(density(log(as.numeric(synt_obs$origin_dist)), bw=x$bw), col=col[color_n])
  }
  
  lines(density(log(as.numeric(synt_simul$origin_dist)), bw=x$bw), col=col[color_n], lty=2)
  
  color_n = color_n + 1
}

legend("topleft", legend=enhancers, fill=col, bty='n')
