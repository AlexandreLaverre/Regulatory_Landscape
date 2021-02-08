
# Proportion of non exonic repated nucleotide
obs$repeat_part <- obs$repet_noexon_bp/obs$length
simul$repeat_part <- simul$repet_noexon_bp/simul$length

# Repeat prop according to distance
obs$dist_class <- cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
simul$dist_class <-cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)

obs_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$repeat_part)*100))
obs_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
obs_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)

obs_enh_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)*100))
obs_enh_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)[["conf.int"]][1]*100)
obs_enh_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$repeat_part)[["conf.int"]][2]*100)

simul_repet_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$repeat_part)*100))
simul_repet_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
simul_repet_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)

