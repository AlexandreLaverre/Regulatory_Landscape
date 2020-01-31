setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

Align_score <- function(sp_origin, sp_target, data){
  setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/Sequence_conservation/")
  dataframe <- read.table(paste(sp_origin, "/contacted_sequences/PIR_cons_all_overlap_PECAN_", sp_origin,"2",sp_target,data,"_merged.txt_onlyconserv", sep=""), header=T)
  
  dataframe$Allexon_ungapped <- dataframe$exclude_ungapped/dataframe$all_exclude
  dataframe[which(dataframe$Allexon_ungapped == "NaN"),]$Allexon_ungapped <- 0
  
  return(dataframe$Allexon_ungapped)
}

sp_origin = "human"
species <- c("macaque", "mouse", "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
#species <- c("rat", "rabbit", "human", "macaque",  "dog", "cow", "elephant",  "opossum", "chicken")

####### A - Dataframe : sequence conservation target to other species ####### 
# Observed
origin_conserv_seq <- read.table(paste("contacted_sequence_composition_", sp_origin,"_merged.txt", sep=""), header=T)
origin_conserv_seq$length <- origin_conserv_seq$end - origin_conserv_seq$start
for (sp_target in species){origin_conserv_seq[[sp_target]] <- Align_score(sp_origin, sp_target, "")}

# Simulated
setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")
origin_conserv_seq_simul <- read.table(paste("contacted_sequence_composition_", sp_origin,"_simul_merged.txt", sep=""), header=T)
origin_conserv_seq_simul$length <- origin_conserv_seq_simul$end - origin_conserv_seq_simul$start
for (sp_target in species){origin_conserv_seq_simul[[sp_target]] <- Align_score(sp_origin, sp_target, "_simul")}

### Filtres ###
# Length
origin_conserv_seq <- origin_conserv_seq[which(origin_conserv_seq$length > 250 & origin_conserv_seq$length < 20000),]
origin_conserv_seq_simul$length <- origin_conserv_seq_simul$end - origin_conserv_seq_simul$start
origin_conserv_seq_simul <- origin_conserv_seq_simul[which(origin_conserv_seq_simul$length > 250 & origin_conserv_seq_simul$length < 20000),]

# Duplication
origin_conserv_seq <- origin_conserv_seq[which(origin_conserv_seq$duplication == 0),]
origin_conserv_seq_simul <- origin_conserv_seq_simul[which(origin_conserv_seq_simul$duplication == 0),]
# Baited
origin_conserv_seq <- origin_conserv_seq[which(origin_conserv_seq$baited == "unbaited"),]
origin_conserv_seq_simul <- origin_conserv_seq_simul[which(origin_conserv_seq_simul$baited == "unbaited"),]

# Uniq cell type
conserv <- origin_conserv_seq[,species]
obs_conserv <- data.frame(score = apply(conserv, 2, mean))
obs_conserv$int_start <- apply(conserv, 2, function(x) t.test(x)[["conf.int"]][1])
obs_conserv$int_end <- apply(conserv, 2, function(x) t.test(x)[["conf.int"]][2])

#Median
#obs_conserv$int_start <- apply(conserv, 2, function(x) boxplot.stats(x)[["conf"]][1])
#obs_conserv$int_end <- apply(conserv, 2, function(x) boxplot.stats(x)[["conf"]][2])

# Observed with enhancers
conserv_enh <- origin_conserv_seq[which(origin_conserv_seq$CAGE_count >1),species]
obs_conserv_enh <- data.frame(score = apply(conserv_enh, 2, mean))
obs_conserv_enh$int_start <- apply(conserv_enh, 2, function(x) t.test(x)[["conf.int"]][1])
obs_conserv_enh$int_end <- apply(conserv_enh, 2, function(x) t.test(x)[["conf.int"]][2])

conserv_simul <- origin_conserv_seq_simul[,species]
simul_conserv <- data.frame(score = apply(conserv_simul, 2, mean))
simul_conserv$int_start <- apply(conserv_simul, 2, function(x) t.test(x)[["conf.int"]][1])
simul_conserv$int_end <- apply(conserv_simul, 2, function(x) t.test(x)[["conf.int"]][2])

####### B - Conservation ~ Distance to promoters #######
species <- c("chicken", "opossum", "rat", "mouse", "rabbit", "dog", "elephant", "cow", "macaque")
#species <- c("chicken", "opossum", "elephant", "cow", "dog", "macaque", "human", "rabbit", "rat")

conserv_dist_all_sp <- list()

for (sp_target in rev(species)){
  origin_conserv_seq$class <-cut(origin_conserv_seq$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  origin_conserv_seq_simul$class <-cut(origin_conserv_seq_simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
  class_leg <- c("50Kb", "500Kb", "1Mb", "1.5Mb", "2Mb","2.5Mb")
  
  # Observed 
  obs_dist<- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) mean(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])))
  obs_dist$int_start <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])[["conf.int"]][1])
  obs_dist$int_end <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])[["conf.int"]][2])
  
  # Observed with enhancers
  obs_dist_enh <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) mean(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])))
  obs_dist_enh$int_start <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])[["conf.int"]][1])
  obs_dist_enh$int_end <- sapply(levels(origin_conserv_seq$class), function(x) t.test(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])[["conf.int"]][2])
  
  # Simulated
  simul_dist <- data.frame(inter = sapply(levels(origin_conserv_seq_simul$class), function(x) mean(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])))
  simul_dist$int_start <- sapply(levels(origin_conserv_seq_simul$class), function(x) t.test(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])[["conf.int"]][1])
  simul_dist$int_end <- sapply(levels(origin_conserv_seq_simul$class), function(x) t.test(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])[["conf.int"]][2])

  conserv_dist_all_sp[[sp_target]] <- list(obs = obs_dist, obs_enh = obs_dist_enh, simul = simul_dist)
  }


####### B - Conservation ~ Distance to promoters (human/mouse) #######
sp_target = "mouse"
origin_conserv_seq$class <-cut(origin_conserv_seq$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
origin_conserv_seq_simul$class <-cut(origin_conserv_seq_simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
class_leg <- c("50Kb", "500Kb", "1Mb", "1.5Mb", "2Mb","2.5Mb")

# Observed 
obs_dist<- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])))
obs_dist$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])[["conf"]][1])
obs_dist$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x),][[sp_target]])[["conf"]][2])

# Observed 1 cell
obs_dist_1 <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$nb_cell ==  1 & origin_conserv_seq$nb_cell <= 3 & origin_conserv_seq$class == x),][[sp_target]])))
obs_dist_1$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$nb_cell ==  1 & origin_conserv_seq$nb_cell <= 3 & origin_conserv_seq$class == x),][[sp_target]])[["conf"]][1])
obs_dist_1$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$nb_cell ==  1 & origin_conserv_seq$nb_cell <= 3 & origin_conserv_seq$class == x),][[sp_target]])[["conf"]][2])

# 1 - 3 cell
obs_dist_2 <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 3 & origin_conserv_seq$nb_cell <= 10 & origin_conserv_seq$class == x),][[sp_target]])))
obs_dist_2$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 3 & origin_conserv_seq$nb_cell <= 10 & origin_conserv_seq$class == x),][[sp_target]])[["conf"]][1])
obs_dist_2$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 3 & origin_conserv_seq$nb_cell <= 10 & origin_conserv_seq$class == x),][[sp_target]])[["conf"]][2])

# # 1 - 3 cell
# obs_dist_3 <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 10 & origin_conserv_seq$nb_cell <= 20 & origin_conserv_seq$class == x),][[sp_target]])))
# obs_dist_3$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 10 & origin_conserv_seq$nb_cell <= 20 & origin_conserv_seq$class == x),][[sp_target]])[["conf"]][1])
# obs_dist_3$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 10 & origin_conserv_seq$nb_cell <= 20 & origin_conserv_seq$class == x),][[sp_target]])[["conf"]][2])

# > 3 cell
obs_dist_4 <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 10 &  origin_conserv_seq$class == x),][[sp_target]])))
obs_dist_4$int_start <- sapply(levels(origin_conserv_seq$class), function(x) 
  tryCatch(boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 10 &  origin_conserv_seq$class == x),][[sp_target]])[["conf"]][1], error=function(e) NA))
obs_dist_4$int_end <- sapply(levels(origin_conserv_seq$class), function(x) 
  tryCatch(boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$nb_cell > 10 &  origin_conserv_seq$class == x),][[sp_target]])[["conf"]][2], error=function(e) NA))

# Observed with enhancers
obs_dist_enh <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])))
obs_dist_enh$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])[["conf"]][1])
obs_dist_enh$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),][[sp_target]])[["conf"]][2])

# Simulated
simul_dist <- data.frame(inter = sapply(levels(origin_conserv_seq_simul$class), function(x) median(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])))
simul_dist$int_start <- sapply(levels(origin_conserv_seq_simul$class), function(x) boxplot.stats(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])[["conf"]][1])
simul_dist$int_end <- sapply(levels(origin_conserv_seq_simul$class), function(x) boxplot.stats(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),][[sp_target]])[["conf"]][2])

#######  C - Conservation between human and mouse ~ Number of cell  #######
obs_nb_cell <- boxplot(origin_conserv_seq$mouse~origin_conserv_seq$nb_type, plot=F)
obs_nb_cell_enh <- boxplot(origin_conserv_seq[which(origin_conserv_seq$CAGE_count >1),]$mouse~origin_conserv_seq[which(origin_conserv_seq$CAGE_count >1),]$nb_cell, plot=F)


#######  D - Overlap with PhastCons element #######
origin_conserv_seq$phastcons_part <- origin_conserv_seq$phastcons_noexonic250 / (origin_conserv_seq$length-origin_conserv_seq$all_exon250)
origin_conserv_seq_simul$phastcons_part <- origin_conserv_seq_simul$phastcons_noexonic250 / (origin_conserv_seq_simul$length-origin_conserv_seq_simul$all_exon250)

# Observed
obs_phastcons <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) 
  median(na.omit(origin_conserv_seq[which(origin_conserv_seq$class == x),]$phastcons_part))))

obs_phastcons$int_start <- sapply(levels(origin_conserv_seq$class), function(x)
  boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x),]$phastcons_part, na.action=na.omit)[["conf"]][1])

obs_phastcons$int_end <- sapply(levels(origin_conserv_seq$class), function(x)
  boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x),]$phastcons_part)[["conf"]][2])

# Observed with enhancers
obs_phastcons_enh <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x)
  median(na.omit(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$phastcons_part))))

obs_phastcons_enh$int_start <- sapply(levels(origin_conserv_seq$class), function(x)
  boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$phastcons_part)[["conf"]][1])
obs_phastcons_enh$int_end <- sapply(levels(origin_conserv_seq$class), function(x) 
  boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$phastcons_part)[["conf"]][2])

# Simulated
simul_phastcons <- data.frame(inter = sapply(levels(origin_conserv_seq_simul$class), function(x) 
  median(na.omit(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$phastcons_part))))
simul_phastcons$int_start <- sapply(levels(origin_conserv_seq_simul$class), function(x) boxplot.stats(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$phastcons_part)[["conf"]][1])
simul_phastcons$int_end <- sapply(levels(origin_conserv_seq_simul$class), function(x) boxplot.stats(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$phastcons_part)[["conf"]][2])


#######  E - Overlap with Repeat element #######
origin_conserv_seq$repeat_part <- origin_conserv_seq$repeat_pb / origin_conserv_seq$length
origin_conserv_seq_simul$repeat_part <- origin_conserv_seq_simul$repeat_pb / origin_conserv_seq_simul$length

# Observed
obs_repeat <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$class == x),]$repeat_part)))
obs_repeat$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x),]$repeat_part)[["conf"]][1])
obs_repeat$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x),]$repeat_part)[["conf"]][2])

# Observed with enhancers
obs_repeat_enh <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$repeat_part)))
obs_repeat_enh$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$repeat_part)[["conf"]][1])
obs_repeat_enh$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$repeat_part)[["conf"]][2])

# Simulated
simul_repeat <- data.frame(inter = sapply(levels(origin_conserv_seq_simul$class), function(x) median(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$repeat_part)))
simul_repeat$int_start <- sapply(levels(origin_conserv_seq_simul$class), function(x) boxplot.stats(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$repeat_part)[["conf"]][1])
simul_repeat$int_end <- sapply(levels(origin_conserv_seq_simul$class), function(x) boxplot.stats(origin_conserv_seq_simul[which(origin_conserv_seq_simul$class == x),]$repeat_part)[["conf"]][2])

save.image("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/Fig4_human_mean.Rdata")




# ## Nb cell ~ Distance
# 
# # Observed
# obs_cell <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$class == x),]$nb_cell)))
# obs_cell$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x),]$nb_cell)[["conf"]][1])
# obs_cell$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x),]$nb_cell)[["conf"]][2])
# 
# # Observed with enhancers
# obs_cell_enh <- data.frame(inter = sapply(levels(origin_conserv_seq$class), function(x) median(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$nb_cell)))
# obs_cell_enh$int_start <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$nb_cell)[["conf"]][1])
# obs_cell_enh$int_end <- sapply(levels(origin_conserv_seq$class), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq$class == x & origin_conserv_seq$CAGE_count > 0),]$nb_cell)[["conf"]][2])
# 
# plot(obs_cell$inter[1:50], type="l", col="red",main=paste("Cell ~ Distance", sep=""),
#      xlab="", ylab="median nb cell", xaxt = "n", ylim=c(1,4))
# for (row in 1:nrow(obs_cell[1:50,])){
#   segments(x0=row,y0=obs_cell[row,]$int_start,x1=row,y1=obs_cell[row,]$int_end, col='red', lwd=0.3)}
# 
# lines(obs_cell_enh$inter[1:50], type="l", col="forestgreen")
# for (row in 1:nrow(obs_cell_enh[1:50,])){
#   segments(x0=row,y0=obs_cell_enh[row,]$int_start,x1=row,y1=obs_cell_enh[row,]$int_end, col='forestgreen', lwd=0.3)}
# 
# axis(1, at=seq(1,51,10), labels=F)
# text(seq(1,51,10),par("usr")[3]-0.2, class_leg, xpd = TRUE)
# 
# 
# 
# ## Conserv / Cell types
# conserv_cell <- data.frame(inter = sapply(colnames(origin_conserv_seq[,18:31]), function(x) median(origin_conserv_seq[which(origin_conserv_seq[[x]] > 0),]$mouse)))
# conserv_cell$int_start <- sapply(colnames(origin_conserv_seq[,18:31]), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq[[x]] > 0),]$mouse)[["conf"]][1])
# conserv_cell$int_end <- sapply(colnames(origin_conserv_seq[,18:31]), function(x) boxplot.stats(origin_conserv_seq[which(origin_conserv_seq[[x]] > 0),]$mouse)[["conf"]][2])
# 
# conserv_cell <- conserv_cell[order(conserv_cell$inter),]
# 
# plot(conserv_cell$inter, type="p", xaxt="n", ylim=c(0.26, 0.32), xlab="", ylab="median Align Score", main="human Conservation / Sample")
# for (row in 1:nrow(conserv_cell)){
#   segments(x0=row,y0=conserv_cell[row,]$int_start,x1=row,y1=conserv_cell[row,]$int_end)}
# 
# axis(1, at=seq(1,14,1), labels=F)
# text(seq(1,14,1),par("usr")[3]-0.008, rownames(conserv_cell), xpd = TRUE, srt=45)
# 
# r2 <- data.frame(inter = sapply(colnames(origin_conserv_seq[,18:31]), function(x) summary(lm(origin_conserv_seq[[x]]~origin_conserv_seq$mouse))$r.squared))
# 
# 
# 
