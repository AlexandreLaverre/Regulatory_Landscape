#########################################################################################################################
library(data.table)
options(stringsAsFactors = FALSE)

source("parameters.R") ## pathFinalDatas are defined based on the user name

ref_sp = "mouse"
minDistance=25e3
maxDistance=2.5e6

enhancers <- c("CAGE", "ENCODE")
if (ref_sp == "human"){enhancers <- c(enhancers, "RoadMap", "GRO_seq")}

############################### Fig2-A - Proportions of contacted sequences which overlap with enhancers ############################### 
data <- c()
conf_up <- c()
conf_low <- c()
comp_test <- c()

### Proportion of the sequences that is enhancer
for (enh in enhancers){
  x <- t.test(obs[[paste0(enh, "_bp")]]/obs$length)
  data <- c(data, x$estimate)
  conf_up <- c(conf_up, x$conf.int[1])
  conf_low <- c(conf_low, x$conf.int[2])
  
  x <- t.test(simul[[paste0(enh, "_bp")]]/simul$length)
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

save(enh_prop, prop_nb_sample, prop_dist, file = paste(pathFigures, "Fig2_", ref_sp, "_A_B_C.Rdata", sep=""))


## OLD ##
# # Proportion of the sequences that is repeat element or phastcons
# stats <- c("repeat_pb", "phastcons_noexonic250")
# for (enh in stats){
#   obs_enh_dist[[paste0(enh)]] <- sapply(levels(obs$dist_class), function(x) 
#     mean(obs[which(obs$dist_class == x),][[paste0(enh)]]/obs[which(obs$dist_class == x ),]$length))
#   obs_enh_dist[[paste0(enh, "_conflow")]] <- sapply(levels(obs$dist_class), function(x)  
#     t.test(obs[which(obs$dist_class == x ),][[paste0(enh)]]/obs[which(obs$dist_class == x ),]$length)[["conf.int"]][1])
#   obs_enh_dist[[paste0(enh, "_confup")]] <- sapply(levels(obs$dist_class), function(x)  
#     t.test(obs[which(obs$dist_class == x ),][[paste0(enh)]]/obs[which(obs$dist_class == x ),]$length)[["conf.int"]][2])
# 
#   obs_enh_dist_with[[paste0(enh)]] <- sapply(levels(obs$dist_class), function(x) 
#     mean(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),][[paste0(enh)]]/obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$length))
#   obs_enh_dist_with[[paste0(enh, "_conflow")]] <- sapply(levels(obs$dist_class), function(x)  
#     tryCatch(t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),][[paste0(enh)]]
#                     /obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$length)[["conf.int"]][1], error=function(e) 0))
#   
#   obs_enh_dist_with[[paste0(enh, "_confup")]] <- sapply(levels(obs$dist_class), function(x)  
#     tryCatch(t.test(obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),][[paste0(enh)]]
#                     /obs[which(obs$dist_class == x & obs$ENCODE_bp > 0),]$length)[["conf.int"]][2],error=function(e) 0))
#   
#   obs_enh_dist_simul[[paste0(enh)]] <- sapply(levels(simul$dist_class), function(x) 
#     mean(simul[which(simul$dist_class == x ),][[paste0(enh)]]/simul[which(simul$dist_class == x ),]$length))
#   obs_enh_dist_simul[[paste0(enh, "_conflow")]] <- sapply(levels(simul$dist_class), function(x)  
#     t.test(simul[which(simul$dist_class == x ),][[paste0(enh)]]/simul[which(simul$dist_class == x ),]$length)[["conf.int"]][1])
#   obs_enh_dist_simul[[paste0(enh, "_confup")]] <- sapply(levels(simul$dist_class), function(x)  
#     t.test(simul[which(simul$dist_class == x ),][[paste0(enh)]]/simul[which(simul$dist_class == x ),]$length)[["conf.int"]][2])
#   
# }
# 
# for (enh in stats){
#   if (enh == "repeat_pb"){YLIM=c(0.30,0.45); YLAB="Repeat proportion"}
#   if (enh == "phastcons_noexonic250"){YLIM=c(0.025,0.1); YLAB="Phastcons proportion"}
#   
#   plot(obs_enh_dist[[enh]][0:50], type="l", col="white", lwd=LWD, ylab=YLAB, ylim=YLIM,
#        xlab="Linear distance to promoters regions (pb)", xaxt = "n", cex.lab=CEX, cex.axis=CEX)
# 
#   points(obs_enh_dist[[paste0(enh)]][0:50], type="l", col="darkgreen", lwd=LWD)
#   for (row in 1:nrow(obs_enh_dist[0:50,])){
#     segments(x0=row,y0=obs_enh_dist[row,paste0(enh, "_conflow")],x1=row,y1=obs_enh_dist[row,paste0(enh, "_confup")], col="darkgreen", lwd=LWD-0.7)}
#   
#   points(obs_enh_dist_with[[paste0(enh)]][0:50], type="l", col="dodgerblue3", lwd=LWD)
#   for (row in 1:nrow(obs_enh_dist[0:50,])){
#     segments(x0=row,y0=obs_enh_dist_with[row,paste0(enh, "_conflow")],x1=row,y1=obs_enh_dist_with[row,paste0(enh, "_confup")], col="dodgerblue3", lwd=LWD-0.7)}
#   
#   points(obs_enh_dist_simul[[paste0(enh)]][0:50], type="l", col="firebrick3", lwd=LWD)
#   for (row in 1:nrow(obs_enh_dist_simul[0:50,])){
#     segments(x0=row,y0=obs_enh_dist_simul[row,paste0(enh, "_conflow")],x1=row,y1=obs_enh_dist_simul[row,paste0(enh, "_confup")], col="firebrick3", lwd=LWD-0.7)}
#   
#   class_leg <- c("25Kb", "500Kb", "1Mb", "1.5Mb", "2Mb", "2.5Mb", "3Mb", "3.5Mb", "4Mb")
#   axis(1, at=seq(1,81,10), labels=F)
#   text(seq(1,81,10), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE, cex=CEX)
#   legend("topright", legend=c("Obs with enh", "Obs", "Simulated"), fill=c("dodgerblue3","forestgreen", "firebrick3"), bty='n', cex=CEX-0.1)
#   
#   
#   }


