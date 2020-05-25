#########################################################################################################################
sp = "human"
path <- "/home/laverre/Data/Regulatory_landscape/result/"
if (sp == "human"){max_sample = 34}else{max_sample=12}

obs <- read.table(paste(path, "Supplementary_dataset1_original_interactions/", sp, "/all_interactions.txt", sep=""), header=T)
simul <- read.table(paste(path, "Supplementary_dataset2_simulated_interactions/", sp, "/simulated_all_interactions.txt", sep=""), header=T)

obs$length <- obs$end-obs$start
obs$nb_sample <- apply(obs, 1, function(x) sum(!is.na(x[9:max_sample])))
obs$sample_class <- cut(obs$nb_sample, breaks=c(1,2,6,11,16,21,max(obs$nb_sample)), include.lowest = T)
obs$dist_class <- cut(obs$distance, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)

simul$length <- simul$end-simul$start
simul$nb_sample <- apply(simul, 1, function(x) sum(!is.na(x[9:max_sample])))
simul$sample_class <- cut(simul$nb_sample, breaks=c(1,2,6,11,16,21,max(simul$nb_sample)), include.lowest = T)
simul$dist_class <- cut(simul$distance, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)

data_list <- list(obs, simul)
data_names <- c("Original", "Simulated")
data_list <- lapply(data_list, function(x) x[which(x$type == "unbaited" & as.character(x$chr_bait) == as.character(x$chr)),])
data_list <- lapply(data_list, function(x) x[which(x$distance < 10000000 & x$distance > 25000),])
# + duplication

#################### Fig 1.B - Histogram with number of samples in which an interaction is observed #####################
# Nb base pairs
nb_pair <- as.data.frame(sapply(data_list, function(x) sapply(levels(as.factor(obs$sample_class)), function(y) 
  sum(x[which(x$sample_class == y),]$length)/sum(x$length))))

proportion <- as.data.frame(sapply(data_list, function(x) sapply(levels(as.factor(obs$sample_class)), function(y) 
  nrow(x[which(x$sample_class == y),])/nrow(x))))

pdf(paste(path, "/Main_figures/Figure1.B_C.pdf", sep=""), width=8.5, height=4)
par(mai = c(1, 0.8, 0.5, 0.1)) # bottom, left, top, right
layout(matrix(c(1,2), nrow = 1, byrow = TRUE))
cut_names = c("1", "2-5", "6-10", "11-15", "16-20", "21-27")

par(lwd=2)
# barplot(t(as.matrix(nb_pair)), beside=T, main="B.", xlab='Sample frequency', names=cut_names, ylim=c(0,1), space=c(0.4,1),
#         ylab="Base pairs proportion", border=c("forestgreen", "firebrick1"), col="white", lwd=1.5, cex.names=0.8)

barplot(t(as.matrix(proportion)), beside=T, main="B.", xlab='Sample frequency', names=cut_names, ylim=c(0,1), space=c(0.4,1),
        ylab="Nb contact proportion", border=c("forestgreen", "firebrick1"), col="white", lwd=1.5, cex.names=0.8)


#################### Fig 1.C - Distribution of number of samples according to distance #####################
# Mean
dist_sample <- as.data.frame(sapply(data_list, function(x) sapply(levels(obs$dist_class), function(y) 
  mean(x[which(x$dist_class == y),]$nb_sample))))

conf_low <- as.data.frame(sapply(data_list, function(x) sapply(levels(obs$dist_class), function(y) 
  tryCatch(t.test(x[which(x$dist_class == y),]$nb_sample)[["conf.int"]][1],error=function(e) 1))))

conf_up <- as.data.frame(sapply(data_list, function(x) sapply(levels(obs$dist_class), function(y) 
  tryCatch(t.test(x[which(x$dist_class == y),]$nb_sample)[["conf.int"]][2],error=function(e) 1))))

# Median
# dist_sample <- as.data.frame(sapply(data_list, function(x) sapply(levels(obs$dist_class), function(y) 
#   median(x[which(x$dist_class == y),]$nb_sample))))
# 
# conf_low <- as.data.frame(sapply(data_list, function(x) sapply(levels(obs$dist_class), function(y) 
#   boxplot.stats(x[which(x$dist_class == y),]$nb_sample)[["conf"]][1])))
# 
# conf_up <- as.data.frame(sapply(data_list, function(x) sapply(levels(obs$dist_class), function(y) 
#   boxplot.stats(x[which(x$dist_class == y),]$nb_sample)[["conf"]][2])))

par(lwd=1)
legend("topright", legend=c("Original", "Simulated"),border=c("forestgreen", "firebrick1"),fill="white", bty='n', cex=0.8)

plot(dist_sample$V1, type="l", col="forestgreen", xaxt = "n", ylim=c(0,6),
     xlab="Distance to promoter region", ylab="Sample frequency (mean)",main="C.")
lines(dist_sample$V2, type="l", col="firebrick1", lwd=1.5)

for (row in 1:nrow(dist_sample)){
  segments(x0=row,y0=conf_low[row,"V1"],x1=row,y1=conf_up[row,"V1"], col="forestgreen", lwd=0.5)}

for (row in 1:nrow(dist_sample)){
  segments(x0=row,y0=conf_low[row,"V2"],x1=row,y1=conf_up[row,"V2"], col="firebrick1", lwd=1)}

dev.off()

