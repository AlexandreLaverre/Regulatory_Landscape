################################################################################################################################################
library(data.table)
options(stringsAsFactors = FALSE)

source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name

ref_sp = "human" # to change human or mouse

minDistance=25e3
maxDistance=15e5

enhancers = c("FANTOM5", "ENCODE")

if(ref_sp == "human"){
  enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")
  target_sp = "mouse"
  closest_sp="macaque"
}else{
  target_sp = "human"
  closest_sp="rat"
}

path_evol <- paste(pathFinalData, "SupplementaryDataset7/", ref_sp, "/", sep="")
path_annot <- paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", sep="")


################################################################################################################################################
################################## Fig3-B Alignment score of contacted sequences between all species ###########################################

#### Contacted restriction fragments
obs <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp, "/statistics_contacted_sequence_original.txt", sep=""), header=T)
simul <- read.table(paste(pathFinalData, "SupplementaryDataset5/", ref_sp,"/statistics_contacted_sequence_simulated.txt", sep=""), header=T)

obs$ID <-  do.call(paste,c(obs[c("chr","start","end")],sep=":"))
simul$ID <-  do.call(paste,c(simul[c("chr","start","end")],sep=":"))

# Filters
obs <- obs[which(obs$baited == "unbaited"),]
simul <- simul[which(simul$baited == "unbaited"),]
obs <- obs[which(obs$BLAT_match < 2),]
simul <- simul[which(simul$BLAT_match < 2),]

# Proportion of non exonic repated nucleotide
obs$repeat_part <- obs$repet_noexon_bp/obs$length
simul$repeat_part <- simul$repet_noexon_bp/simul$length

# Repeat prop according to distance
obs$dist_class <- cut(obs$median_dist, breaks=seq(from=minDistance, to=maxDistance+25000, by=50000), include.lowest = T)
simul$dist_class <-cut(simul$median_dist, breaks=seq(from=minDistance, to=maxDistance+25000, by=50000), include.lowest = T)

obs_repet_dist <- data.frame(inter = sapply(levels(obs$dist_class), function(x) mean(obs[which(obs$dist_class == x),]$repeat_part)*100))
obs_repet_dist$int_start <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
obs_repet_dist$int_end <- sapply(levels(obs$dist_class), function(x) t.test(obs[which(obs$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)

simul_repet_dist <- data.frame(inter = sapply(levels(simul$dist_class), function(x) mean(simul[which(simul$dist_class == x),]$repeat_part)*100))
simul_repet_dist$int_start <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][1]*100)
simul_repet_dist$int_end <- sapply(levels(simul$dist_class), function(x) t.test(simul[which(simul$dist_class == x),]$repeat_part)[["conf.int"]][2]*100)


########################  Repeat proportion vs distance to promoters ######################## 
if(ref_sp=="human"){ymin=30; ymax=50}else{ymin=25; ymax=50}
class_leg <- c("0", "0.5", "1", "1.5")
# Original
pdf()
plot(obs_repet_dist[,"inter"], type="l", col=dataset.colors["Original"], main="", axes=F, cex.lab=1.2,
     xlab="", ylab="% length covered by repeat elements", xaxt = "n", ylim=c(ymin,ymax), las=2, , lwd=1.1)

for (row in 1:nrow(obs_repet_dist)){
  segments(x0=row,y0=obs_repet_dist[row,"int_start"],x1=row,y1=obs_repet_dist[row,"int_end"], col=dataset.colors["Original"])}

# Simulated
lines(simul_repet_dist[,"inter"], type="l", col=dataset.colors["Simulated"], lwd=1.1)
for (row in 1:nrow(simul_repet_dist)){
  segments(x0=row,y0=simul_repet_dist[row,"int_start"],x1=row,y1=simul_repet_dist[row,"int_end"], col=dataset.colors["Simulated"])}

# Axis and legend

axis(side=2, mgp=c(3, 0.75, 0), las=2, cex.axis=1.1)
axis(side=1, at=seq(1,nrow(obs_repet_dist)+1,10), mgp=c(3, 0.75, 0), labels=class_leg, cex.axis=1.1)

mtext("distance to promoters (Mb)", side=1, line=2.5, cex=1.2)
