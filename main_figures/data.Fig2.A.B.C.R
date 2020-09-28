#########################################################################################################################
library(data.table)
options(stringsAsFactors = FALSE)

source("parameters.R") ## pathFinalData defined based on the user name

ref_sp = "human"
minDistance=25e3
maxDistance=2e6

enhancers <- c("ENCODE", "FANTOM5")
if (ref_sp == "human"){enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")}

obs <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp, "/statistics_contacted_sequence_original.txt", sep=""), header=T)
simul <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp,"/statistics_contacted_sequence_simulated.txt", sep=""), header=T)
obs <- obs[which(obs$baited == "unbaited" & obs$BLAT_match == 1),]
simul <- simul[which(simul$baited == "unbaited" & simul$BLAT_match == 1),]

############################### Fig2-A - Proportions of contacted sequences which overlap with enhancers ############################### 
data <- c()
conf_up <- c()
conf_low <- c()
comp_test <- c()

### Proportion of the sequences that is enhancer
for (enh in enhancers){
  x <- t.test(obs[[paste0(enh, "_bp")]]*100/obs$length)
  data <- c(data, x$estimate)
  conf_up <- c(conf_up, x$conf.int[1])
  conf_low <- c(conf_low, x$conf.int[2])
  
  x <- t.test(simul[[paste0(enh, "_bp")]]*100/simul$length)
  data <- c(data, x$estimate, 0)
  conf_up <- c(conf_up, x$conf.int[1], -1)
  conf_low <- c(conf_low, x$conf.int[2], -1)
  
  comp_test <- c(comp_test, prop.test(x = c(nrow(obs[which(obs[[paste0(enh, "_bp")]] == TRUE),]), nrow(obs[which(simul[[paste0(enh, "_bp")]] == TRUE),])),
                                      n = c(nrow(obs), nrow(simul)))$p.value)
}

enh_prop <- data.frame(data=data, conf_up=conf_up, conf_low=conf_low)

###############################  Fig2-B - According to distance from promoters ############################################################### 
obs$dist_class <-cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
simul$dist_class <- cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)

obs_enh_dist <- data.frame(matrix(vector(), length(levels(obs$dist_class)), 1))
simul_enh_dist <- data.frame(matrix(vector(), length(levels(simul$dist_class)), 1)) 

# Proportion of the sequences that is enh
for (enh in enhancers){
  obs_enh_dist[[paste0(enh)]] <- sapply(levels(obs$dist_class), function(x)
    mean(obs[which(obs$dist_class == x),][[paste0(enh, "_bp")]]/obs[which(obs$dist_class == x),]$length))
  obs_enh_dist[[paste0(enh, "_conflow")]] <- sapply(levels(obs$dist_class), function(x)  
    t.test(obs[which(obs$dist_class == x),][[paste0(enh, "_bp")]]/obs[which(obs$dist_class == x),]$length)[["conf.int"]][1])
  obs_enh_dist[[paste0(enh, "_confup")]] <- sapply(levels(obs$dist_class), function(x)  
    t.test(obs[which(obs$dist_class == x),][[paste0(enh, "_bp")]]/obs[which(obs$dist_class == x),]$length)[["conf.int"]][2])
  
  simul_enh_dist[[paste0(enh)]] <- sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),][[paste0(enh, "_bp")]]/simul[which(simul$dist_class == x),]$length))
  simul_enh_dist[[paste0(enh, "_conflow")]] <- sapply(levels(simul$dist_class), function(x)  
    t.test(simul[which(simul$dist_class == x),][[paste0(enh, "_bp")]]/simul[which(simul$dist_class == x),]$length)[["conf.int"]][1])
  simul_enh_dist[[paste0(enh, "_confup")]] <- sapply(levels(simul$dist_class), function(x)  
    t.test(simul[which(simul$dist_class == x),][[paste0(enh, "_bp")]]/simul[which(simul$dist_class == x),]$length)[["conf.int"]][2])
  
}

prop_dist <- list(obs=obs_enh_dist, simul=simul_enh_dist)

###############################  Fig2-C - According to number of samples ############################################################### 
obs$nb_cell <- as.factor(obs$nb_sample)
obs_enh_cell <- data.frame(matrix(vector(), length(levels(obs$nb_cell)), 1)) 
simul$nb_cell <- as.factor(simul$nb_sample)
simul_enh_cell <- data.frame(matrix(vector(), length(levels(simul$nb_cell)), 1)) 

#### Proportion of the sequences that is enhancer
for (enh in enhancers){
  obs_enh_cell[[paste0(enh)]] <- sapply(levels(obs$nb_cell), function(x) 
    mean(obs[which(obs$nb_cell == x),][[paste0(enh, "_bp")]]/obs[which(obs$nb_cell == x),]$length))
  
  obs_enh_cell[[paste0(enh, "_conflow")]] <- sapply(levels(obs$nb_cell), function(x)
    tryCatch(t.test(obs[which(obs$nb_cell == x),][[paste0(enh, "_bp")]]/obs[which(obs$nb_cell == x),]$length)
             [["conf.int"]][1], error=function(e) 0))
  
  obs_enh_cell[[paste0(enh, "_confup")]] <- sapply(levels(obs$nb_cell), function(x)
    tryCatch(t.test(obs[which(obs$nb_cell == x),][[paste0(enh, "_bp")]]/obs[which(obs$nb_cell == x),]$length)
             [["conf.int"]][2], error=function(e) 0))
  
  simul_enh_cell[[paste0(enh)]] <- sapply(levels(simul$nb_cell), function(x) 
    mean(simul[which(simul$nb_cell == x),][[paste0(enh, "_bp")]]/simul[which(simul$nb_cell == x),]$length))
  
  simul_enh_cell[[paste0(enh, "_conflow")]] <- sapply(levels(simul$nb_cell), function(x)
    tryCatch(t.test(simul[which(simul$nb_cell == x),][[paste0(enh, "_bp")]]/simul[which(simul$nb_cell == x),]$length)
             [["conf.int"]][1], error=function(e) 0))
  
  simul_enh_cell[[paste0(enh, "_confup")]] <- sapply(levels(simul$nb_cell), function(x)
    tryCatch(t.test(simul[which(simul$nb_cell == x),][[paste0(enh, "_bp")]]/simul[which(simul$nb_cell == x),]$length)
             [["conf.int"]][2], error=function(e) 0))
  
}

obs_enh_cell <- obs_enh_cell[-1]
rownames(obs_enh_cell) <- levels(obs$nb_cell)
simul_enh_cell <- simul_enh_cell[-1]
rownames(simul_enh_cell) <- levels(simul$nb_cell)

prop_nb_sample <- list(obs=obs_enh_cell, simul=simul_enh_cell)

################################################# Save RData ################################################# 

save(enh_prop, prop_nb_sample, prop_dist, file = paste(pathFigures, "RData/Fig2_", ref_sp, "_A_B_C.Rdata", sep=""))
