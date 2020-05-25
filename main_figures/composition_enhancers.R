options(stringsAsFactors = FALSE)
ref_sp = "human"
target_sp = "mouse" 
enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")
col <- c("red", "navy", "forestgreen", "orange")

path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, sep="/")

pdf(paste(path_evol, "/", ref_sp, "_enhancers_composition_by_", var, ".pdf", sep=""), width=8, height=6)

colfunc <- colorRampPalette(c("gold", "seagreen"))
colors <- colfunc(4)
colors.alpha <- apply(sapply(colors, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=0.2))  
colors.alpha2 <- apply(sapply(colors, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=0.1))  
width <- .8

list_obs <- list()
list_simul <- list()
for (enh_data in enhancers){
  enh_obs <- read.table(paste(path_annot,"/contacted_", enh_data, "_composition_", ref_sp, "_observed.txt", sep=""), header=T, sep="\t")
  enh_simul <- read.table(paste(path_annot,"/contacted_", enh_data, "_composition_", ref_sp, "_simulated.txt", sep=""), header=T, sep="\t")
  
  enh_obs$prop_repet <- enh_obs$repeat_pb/enh_obs$length
  enh_simul$prop_repet <- enh_simul$repeat_pb/enh_simul$length
  
  enh_obs$prop_exon <- enh_obs$all_exon_pb/(enh_obs$length-enh_obs$repeat_pb)
  enh_simul$prop_exon <- enh_simul$all_exon_pb/(enh_simul$length-enh_simul$repeat_pb)
  
  enh_obs$prop_GC <- enh_obs$GC_pb/(enh_obs$length-enh_obs$repeat_pb)
  enh_simul$prop_GC <- enh_simul$GC_pb/(enh_simul$length-enh_simul$repeat_pb)
  
  # Take only not duplicated enhancers
  list_obs[[enh_data]] <- enh_obs[which(enh_obs$duplication <2),]
  list_simul[[enh_data]] <- enh_simul[which(enh_simul$duplication <2),]
}
  
variables <- colnames(enh_obs[4:13])

for (var in variables){
  boxplot(c(list_obs[[1]][var], list_simul[[1]][var],
            list_obs[[2]][var], list_simul[[2]][var],
            list_obs[[3]][var], list_simul[[3]][var],
            list_obs[[4]][var], list_simul[[4]][var]),
          at=c(0.8,1.2, 1.8,2.2, 2.8,3.2, 3.8,4.2),
          outline=T, notch=T, ylab=var, boxwex=0.2, xaxt='n', main=var)
  
  for(c in c(1:length(list_obs))){
    # obs
    d  <- density(list_obs[[c]][[var]], adjust = 1)
    polygon(c(c-d$y/max(d$y)*width/2, rep(c,length(d$y))), c(d$x,rev(d$x)), border="forestgreen", col=rgb(t(col2rgb('forestgreen')/255), alpha = 0.6), lwd=2)

    #simul
    d  <- density(list_simul[[c]][[var]], adjust = 1)
    polygon(c(c+d$y/max(d$y)*width/2, rep(c,length(d$y))), c(d$x,rev(d$x)), border="firebrick1", col=rgb(t(col2rgb('firebrick1')/255), alpha = 0.6), lwd=2)
  }
  
  axis(1, at=seq(1,4,1), labels=enhancers)
  legend("topleft", fill=c("forestgreen", "firebrick1"), legend=c("Observed", "Simulated"), bty='n')
  
  
}

  
  
