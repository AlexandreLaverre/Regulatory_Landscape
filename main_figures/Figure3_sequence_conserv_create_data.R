library(data.table)
options(stringsAsFactors = FALSE)

ref_sp = "human"
target_sp = "mouse"

path <- "/home/laverre/Data/Regulatory_landscape/result/"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, "/", sep="/")

####################################### Restriction fragments datas  ####################################### 
obs <- read.table(paste(path_annot, "/contacted_sequence_composition_", ref_sp, "_observed.txt", sep=""), header=T)
simul <- read.table(paste(path_annot,"/contacted_sequence_composition_", ref_sp, "_simulated.txt", sep=""), header=T)

obs <- obs[which(obs$baited == "unbaited"),]
simul <- simul[which(simul$baited == "unbaited"),]
obs <- obs[which(obs$duplication < 1),]
simul <- simul[which(simul$duplication < 1),]

obs$ID <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
simul$ID <-  do.call(paste,c(simul[c("chr","start","end")],sep=":"))


# Correlation ENCODE ~ align
# x=obs$all_exon_pb
# y=obs$align_target
# R=cor(x, y, method="pearson")
# rho=cor(x, y, method="spearman")
# 
# smoothScatter(x, y)
# mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5)
# abline(lm(y~x), col="red")


################# Alignment score vs all species #################
align <- read.table(paste(path_evol,"sequence_conservation/restriction_fragments/Alignments_stats_all_species_total_ungapped.txt", sep="/"), header=T)

species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")

align_obs <- align[which(align$enh %in% obs$ID), c("enh", species)]
align_obs_enh <- align[which(align$enh %in% obs[which(obs$ENCODE_bp > 0),]$ID), c("enh", species)]
align_simul <- align[which(align$enh %in% simul$ID), c("enh", species) ]


################# Repeat proportion  #################
obs$repeat_part <- obs$repeat_pb/obs$length
simul$repeat_part <- simul$repeat_pb/simul$length

# Repeat prop vs BLAT
# par(mfrow=c(3,1))
# hist(obs$repeat_part, col="red", xlab="Repeat proportion", main="Observed contacted fragment")
# hist(obs[which(obs$duplication == 0),]$repeat_part, add=TRUE, col="white")
# legend("topright", fill="red", legend ="BLAT =! 1", bty='n')
# 
# hist(simul$repeat_part, col="red", xlab="Repeat proportion", main="Simulated contacted fragment")
# hist(simul[which(simul$duplication == 0),]$repeat_part, add=TRUE, col="white")
# 
# hist(obs[which(obs$ENCODE_bp > 0),]$repeat_part, col="red", xlab="Repeat proportion", main="Observed with ENCODE contacted fragment")
# hist(obs[which(obs$ENCODE_bp > 0),]$repeat_part, add=TRUE, col="white")

# Repeat prop according to distance
obs$dist_class <- cut(obs$midist, breaks=seq(from=25000, to=2025000, by=50000), include.lowest = T)
simul$dist_class <-cut(simul$midist, breaks=seq(from=25000, to=2025000, by=50000), include.lowest = T)

obs_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$repeat_part)*100))
obs_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
obs_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)

obs_enh_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)*100))
obs_enh_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)[["conf.int"]][1]*100)
obs_enh_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)[["conf.int"]][2]*100)

simul_repet_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$repeat_part)*100))
simul_repet_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
simul_repet_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)

#################  Exonic prop according to distance #################
obs$exonic_part <- obs$all_exon_pb/obs$length
simul$exonic_part <- simul$all_exon_pb/simul$length

obs_exon_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$exonic_part)*100))
obs_exon_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$exonic_part)[["conf.int"]][1]*100)
obs_exon_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$exonic_part)[["conf.int"]][2]*100)

obs_enh_exon_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$exonic_part)*100))
obs_enh_exon_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$exonic_part)[["conf.int"]][1]*100)
obs_enh_exon_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$exonic_part)[["conf.int"]][2]*100)

simul_exon_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$exonic_part)*100))
simul_exon_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$exonic_part)[["conf.int"]][1]*100)
simul_exon_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$exonic_part)[["conf.int"]][2]*100)


################# Sequence conservation human 2 mouse vs distance to promoters #################
#obs <- obs[which(obs$duplication == 0),]
#simul <- simul[which(simul$duplication == 0),]

obs <- obs[order(obs$ID), ]
simul <- simul[order(simul$ID), ]

obs$align_target <- align_obs[order(align_obs$enh), target_sp]
simul$align_target <- align_simul[order(align_simul$enh), target_sp]

obs_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$align_target)))
obs_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][2])

obs_enh_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)))
obs_enh_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][1])
obs_enh_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$align_target)))
simul_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][2])

################# Sequence conservation human 2 macaque vs distance to promoters #################
obs$align_target <- align_obs[order(align_obs$enh), "macaque"]
simul$align_target <- align_simul[order(align_simul$enh), "macaque"]

obs_dist_mac <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$align_target)))
obs_dist_mac$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][1])
obs_dist_mac$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$align_target)[["conf.int"]][2])

obs_enh_dist_mac <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)))
obs_enh_dist_mac$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][1])
obs_enh_dist_mac$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$align_target)[["conf.int"]][2])

simul_dist_mac <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$align_target)))
simul_dist_mac$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][1])
simul_dist_mac$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$align_target)[["conf.int"]][2])



####################################### Enhancers datas  ####################################### 
enh = "ENCODE"

align <- read.table(paste(path_evol,"sequence_conservation", enh, "Alignments_stats_all_species_total_ungapped.txt", sep="/"), header=T)
obs_stats <- read.table(paste(path_evol,"/sequence_conservation/", enh, "/original_stats.txt", sep=""), header=T)
simul_stats <- read.table(paste(path_evol,"/sequence_conservation/", enh, "/simulated_stats.txt", sep=""), header=T)
enh_annot <- read.table(paste(path_annot,"/", enh, "_BLAT_summary_0.8.txt", sep=""), header=T, sep="\t")

# Filters
species <- c("macaque", target_sp, "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
align <- align[, c("enh", species)]

enh_annot <- enh_annot[which(enh_annot$nb_match < 2),]
align_enhancers_obs <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% obs_stats$enh),]
align_enhancers_simul <- align[which(align$enh %in% enh_annot[,1] & align$enh %in% simul_stats$enh),]

# Output
save(align_obs, align_obs_enh, align_simul,
     obs_repet_dist, obs_enh_repet_dist, simul_repet_dist,
     obs_exon_dist, obs_enh_exon_dist, simul_exon_dist,
     obs_dist, obs_enh_dist, simul_dist,
     obs_dist_mac,obs_enh_dist_mac,simul_dist_mac,
     align_enhancers_obs, align_enhancers_simul,
     file = paste(path, "/Main_figures/Fig3_seq_conserv_", ref_sp, ".Rdata", sep=""))



# Tests
# pdf(paste(path, "/Align_restriction_fragment_distance_all_species_BLAT.pdf", sep=""), width=7, height=5)
# par(mfrow=c(3,3), mai = c(0.3, 0.7, 0.3, 0.2))
# layout( matrix(c(1,2,3,4,5,6,7,8,9),nrow = 3,ncol = 3,byrow = TRUE))


# CEX=1
# CEX_lines=1
# xmin=0 #min(c(simul_dist$int_start, obs_enh$int_start))
# xmax=0.55 #max(obs_enh_dist$inter)
# 
# plot(obs_dist[,"inter"], type="l", col="firebrick1", cex=CEX_lines, main=paste(ref_sp, " to ", target_sp, sep=""),
#      xlab="", ylab="Ungapped Non-exonic Score", xaxt = "n", ylim=c(xmin,xmax), cex.lab=CEX, cex.axis=CEX, cex.main=CEX)
# 
# for (row in 1:nrow(obs_dist)){
#   segments(x0=row,y0=obs_dist[row,"int_start"],x1=row,y1=obs_dist[row,"int_end"], col='firebrick1', lwd=0.3)}
# 
# lines(simul_dist[,"inter"], type="l", col="blue", cex=CEX_lines)
# for (row in 1:nrow(simul_dist)){
#   segments(x0=row,y0=simul_dist[row,"int_start"],x1=row,y1=simul_dist[row,"int_end"], col='blue', lwd=0.3)}
# 
# lines(obs_enh_dist[,"inter"], type="l", col="forestgreen", cex=CEX_lines)
# for (row in 1:nrow(obs_enh_dist)){
#   segments(x0=row,y0=obs_enh_dist[row,"int_start"],x1=row,y1=obs_enh_dist[row,"int_end"], col='forestgreen', lwd=0.3)}
# 
# axis(1, at=seq(1,nrow(obs_dist)+1,10), labels=F)
# text(seq(1,nrow(obs_dist)+1,10),par("usr")[3]-0.02, class_leg, xpd = TRUE, cex=CEX)


#dev.off()

