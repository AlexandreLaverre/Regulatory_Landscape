##### DATA #####
setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

sp_origin = "human"
species <- c("macaque",  "dog", "cow", "elephant", "rabbit", "rat", "mouse", "opossum", "chicken")

obs_PIR <- read.table(paste("contacted_sequence_composition_", sp_origin,"_merged.txt",sep=""), header=T)
obs_PIR$PIR <- paste(obs_PIR$chr, obs_PIR$start, obs_PIR$end, sep=':')

data <- c()
data_2M <- c()
p_test <- c()
p_test_enh <- c()

for (sp_target in species){
  setwd(paste("/home/laverre/Documents/Regulatory_Landscape/result/conservation/syntenie_conservation/", sp_origin, sep=''))
  obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv.txt", sep=""), header=T)
  simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv_simul.txt", sep=""), header=T)
  
  #### Filtres
  # Length
  obs <- obs[which(obs$bait_length > 250 & obs$bait_length < 20000 & obs$PIR_length >250 & obs$PIR_length < 20000),]
  simul <- simul[which(simul$bait_length > 250 & simul$bait_length < 20000 & simul$PIR_length >250 & simul$PIR_length < 20000),]
  # Duplication
  obs <- obs[which(obs$bait_dupli == 0 & obs$PIR_dupli == 0),]
  simul <- simul[which(simul$bait_dupli == 0 & simul$PIR_dupli == 0),]
  
  obs <- obs[which(!is.na(obs$bait_lift)),]
  obs <- obs[which(!is.na(obs$PIR_lift)),]
  simul <- simul[which(!is.na(simul$bait_lift)),]
  simul <- simul[which(!is.na(simul$PIR_lift)),]
  
  obs$PIR <- as.factor(sub(".*-", "", obs$origin_interaction))
  obs$CAGE <- obs$PIR %in%  obs_PIR[which(obs_PIR$CAGE_count > 0),]$PIR
  
  synt_obs = nrow(obs[which(!is.na(obs$target_dist)),])
  synt_obs_2M = nrow(obs[which(obs$target_dist < 2000000),])
  synt_obs_CAGE = nrow(obs[which(!is.na(obs$target_dist) & obs$CAGE == TRUE),])
  synt_obs_CAGE_2M = nrow(obs[which(obs$target_dist < 2000000 & obs$CAGE == TRUE),])
  synt_simul = nrow(simul[which(!is.na(simul$target_dist)),])
  synt_simul_2M = nrow(simul[which(simul$target_dist < 2000000),])
  
  mat <- matrix(c(synt_obs, synt_simul, nrow(obs)-synt_obs, nrow(simul)-synt_simul),2)
  mat_enh <- matrix(c(synt_obs, synt_obs_CAGE, nrow(obs)-synt_obs, nrow(obs[which(obs$CAGE == TRUE),])-synt_obs_CAGE),2)
  p_test <- append(p_test, prop.test(mat)$p.value)
  p_test_enh <- append(p_test_enh, prop.test(mat_enh)$p.value)
  
  data <- append(data, c( (synt_simul/nrow(simul))*100, (synt_obs/nrow(obs))*100, (synt_obs_CAGE/nrow(obs[which(obs$CAGE == TRUE),]))*100, NA))
  data_2M <- append(data_2M, c( (synt_simul_2M/nrow(simul))*100, (synt_obs_2M/nrow(obs))*100, (synt_obs_CAGE_2M/nrow(obs[which(obs$CAGE == TRUE),]))*100, NA))
  
  message(paste(sp_target,"done !", sep=" "))
}


### C - Syntenie conservation ~~ genomic distance
sp_target = "chicken"
setwd(paste("/home/laverre/Documents/Regulatory_Landscape/result/conservation/syntenie_conservation/", sp_origin, sep=''))
obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv.txt", sep=""), header=T)
simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv_simul.txt", sep=""), header=T)

#### Filtres
# Length
obs <- obs[which(obs$bait_length > 250 & obs$bait_length < 20000 & obs$PIR_length >250 & obs$PIR_length < 20000),]
simul <- simul[which(simul$bait_length > 250 & simul$bait_length < 20000 & simul$PIR_length >250 & simul$PIR_length < 20000),]
# Duplication
obs <- obs[which(obs$bait_dupli == 0 & obs$PIR_dupli == 0),]
simul <- simul[which(simul$bait_dupli == 0 & simul$PIR_dupli == 0),]

obs <- obs[which(!is.na(obs$bait_lift)),]
obs <- obs[which(!is.na(obs$PIR_lift)),]
simul <- simul[which(!is.na(simul$bait_lift)),]
simul <- simul[which(!is.na(simul$PIR_lift)),]

obs$PIR <- as.factor(sub(".*-", "", obs$origin_interaction))
obs$CAGE <- obs$PIR %in%  obs_PIR[which(obs_PIR$ENCODE_count > 0),]$PIR

obs$class <-cut(obs$origin_dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
simul$class <-cut(simul$origin_dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
class_leg <- c("50Kb", "500Kb", "1Mb", "1.5Mb", "2Mb","2.5Mb")

### Cis conserved
# Observed 
obs_dist<- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),])/nrow(obs[which(obs$class == x),]))*100))
obs_dist$int_start <- sapply(levels(obs$class), function(x) (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x),]), p=0.5)$conf.int[1])*100)
obs_dist$int_end <- sapply(levels(obs$class), function(x) (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x),]), p=0.5)$conf.int[2])*100)

# CAGE
obs_dist_enh<- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & obs$CAGE == TRUE & !is.na(obs$target_dist)),])/nrow(obs[which(obs$class == x & obs$CAGE == TRUE),]))*100))
obs_dist_enh$int_start <- sapply(levels(obs$class), function(x) (prop.test(x = nrow(obs[which(obs$class == x & obs$CAGE == TRUE & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x & obs$CAGE == TRUE),]), p=0.5)$conf.int[1])*100)
obs_dist_enh$int_end <- sapply(levels(obs$class), function(x) (prop.test(x = nrow(obs[which(obs$class == x & obs$CAGE == TRUE & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x & obs$CAGE == TRUE),]), p=0.5)$conf.int[2])*100)

# simul
simul_dist<- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),])/nrow(simul[which(simul$class == x),]))*100))
simul_dist$int_start <- sapply(levels(simul$class), function(x) (prop.test(x = nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(simul$class), function(x) (prop.test(x = nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[2])*100)


### Cis + < 2Mb conserv
obs_dist_2M <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x &  obs$target_dist < 2000000),])/nrow(obs[which(obs$class == x),]))*100))
obs_dist_2M$int_start <- sapply(levels(obs$class), function(x) (prop.test(x = nrow(obs[which(obs$class == x & obs$target_dist < 2000000),]), n=nrow(obs[which(obs$class == x),]), p=0.5)$conf.int[1])*100)
obs_dist_2M$int_end <- sapply(levels(obs$class), function(x) (prop.test(x = nrow(obs[which(obs$class == x & obs$target_dist < 2000000),]), n=nrow(obs[which(obs$class == x),]), p=0.5)$conf.int[2])*100)

# CAGE
obs_dist_enh_2M <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & obs$CAGE == TRUE & obs$target_dist < 2000000),])/nrow(obs[which(obs$class == x & obs$CAGE == TRUE),]))*100))
obs_dist_enh_2M$int_start <- sapply(levels(obs$class), function(x) (prop.test(x = nrow(obs[which(obs$class == x & obs$CAGE == TRUE & obs$target_dist < 2000000),]), n=nrow(obs[which(obs$class == x & obs$CAGE == TRUE),]), p=0.5)$conf.int[1])*100)
obs_dist_enh_2M$int_end <- sapply(levels(obs$class), function(x) (prop.test(x = nrow(obs[which(obs$class == x & obs$CAGE == TRUE & obs$target_dist < 2000000),]), n=nrow(obs[which(obs$class == x & obs$CAGE == TRUE),]), p=0.5)$conf.int[2])*100)

# simul
simul_dist_2M <- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & simul$target_dist < 2000000),])/nrow(simul[which(simul$class == x),]))*100))
simul_dist_2M$int_start <- sapply(levels(simul$class), function(x) (prop.test(x = nrow(simul[which(simul$class == x & simul$target_dist < 2000000),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_dist_2M$int_end <- sapply(levels(simul$class), function(x) (prop.test(x = nrow(simul[which(simul$class == x & simul$target_dist < 2000000),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[2])*100)

save.image("/home/laverre/Documents/Regulatory_Landscape/scripts/main_figures/Fig6_human.Rdata")

