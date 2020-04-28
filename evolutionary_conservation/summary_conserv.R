setwd("/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset6_regulatory_landscape_evolution/")

ref_sp = "human"
target_sp = "mouse" 
enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")

#pdf(paste(ref_sp, "2", target_sp, "_summary_conserv_PhastCons_0.05_by_orthologue_gene.pdf", sep=""))
for (enh in enhancers){
  obs <- read.table(paste(ref_sp,"/", enh, "_original_summary_conserv_0.5.txt", sep=""), header=T)
  #obs <- read.table(paste(ref_sp,"/", enh, "_original_summary_conserv_PhastCons_0.05.txt", sep=""), header=T)
  simul <- read.table(paste(ref_sp,"/", enh, "_simulated_summary_conserv_0.5.txt", sep=""), header=T)
  #simul <- read.table(paste(ref_sp,"/", enh, "_simulated_summary_conserv_PhastCons_0.05.txt", sep=""), header=T)
  data_list <- list(obs, simul)
  data_names <- c("Original", "Simulated")

  par(mai = c(0.7, 0.9, 0.7, 0.5))
  layout(matrix(c(1,2,3,4),nrow = 2, byrow = TRUE))
  
  enh_cover <- sapply(data_list, function(x) x$total_enh_length/x$contacted_coverage)
  B <- boxplot(enh_cover, main=paste(ref_sp, enh, "\n contacted",sep=" "), names=data_names, border=c("darkgreen", "firebrick3"),
               ylab="Proportion \n (enh length / contacted length)", outline=F, notch=T, boxwex=0.5)
  points(B$group, B$out, type="p", pch=1, cex=0.2, col=c("darkgreen", "firebrick3")[B$group])
  
  enh_conserv <- sapply(data_list, function(x) x$nb_seq_conserv/x$nb_total)
  B <- boxplot(enh_conserv, main=paste(ref_sp, enh, "\n conserved in sequence", sep=" "), names=data_names,
               ylab="Proportion \n (align score > 0.8)", outline=F, notch=T, boxwex=0.5, border=c("darkgreen", "firebrick3"))
  points(B$group, B$out, type="p", pch=1, cex=0.2, col=c("darkgreen", "firebrick3")[B$group])
  
  enh_synt2M <- sapply(data_list, function(x) x$nb_synt2M_conserv/x$nb_seq_conserv)
  B <- boxplot(enh_synt2M, main=paste(ref_sp, enh, "\n maintained in synteny (2M)", sep=" "), names=data_names,
               ylab="Proportion", outline=F, notch=T, boxwex=0.5, border=c("darkgreen", "firebrick3"))
  points(B$group, B$out, type="p", pch=1, cex=0.2, col=c("darkgreen", "firebrick3")[B$group])
  
  enh_contact <- sapply(data_list, function(x) x$nb_contact_conserv/x$nb_seq_conserv)
  B <- boxplot(enh_contact, main=paste(ref_sp, enh, "\n maintained in contact", sep=" "), names=data_names,
               ylab="Proportion", outline=F, notch=T, boxwex=0.5, border=c("darkgreen", "firebrick3"))
  points(B$group, B$out, type="p", pch=1, cex=0.2, col=c("darkgreen", "firebrick3")[B$group])

}

#dev.off()



# Fig 2.A : prop overlap enhancer all datasets 
# Plot empty boxplot
pdf(paste(ref_sp, "_sequence_overlap_enhancer_log10.pdf", sep=""), width=7, height=5)
par(mai = c(0.7, 1.1, 0.7, 0.5))
boxplot(seq(0,0.3, 0.1)~enhancers,
        boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), names=enhancers,
        main="", ylab="Proportion \n log10(enh length / contacted length)", xlab="", ylim=c(-4,0))
x=1

for (enh in enhancers){
  obs <- read.table(paste(ref_sp,"/", enh, "_original_summary_conserv_0.8.txt", sep=""), header=T)
  simul <- read.table(paste(ref_sp,"/", enh, "_simulated_summary_conserv_0.8.txt", sep=""), header=T)
  data_list <- list(obs, simul)
  data_names <- c("Original", "Simulated")
  

  enh_cover <- sapply(data_list, function(x) log10(x$total_enh_length/x$contacted_coverage))
  boxplot(enh_cover, names=data_names, add=T, xaxt = "n", yaxt='n', outline=T, notch=T,
          border=c("darkgreen", "firebrick3"),  boxwex=0.25, at = c(x-0.2, x+0.2), outcex=0.1)
  
  x=x+1
  
}

legend("bottomright", border=c("darkgreen", "firebrick3"), fill=c("white","white"), legend = c("original", "simulated"), bty='n')
dev.off()
