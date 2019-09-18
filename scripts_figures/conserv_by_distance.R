################################### conservation de séquence par classe de distance ###################################
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/human2mouse/")
obs_pecan <- read.table("PIR_cons_all_overlap_PECAN.txt", header=T)
simul_pecan <- read.table("PIR_cons_all_overlap_PECAN_simul.txt", header=T)

obs_tba <- read.table("PIR_cons_all_overlap_TBA.txt", header=T)
simul_tba <- read.table("PIR_cons_all_overlap_TBA_simul.txt", header=T)

obs <- obs_tba
simul <- simul_tba

#### Filtres ###
# Length
obs <- obs[which(obs$length > 250 & obs$length < 20000),]
simul <- simul[which(simul$length > 250 & simul$length < 20000),]

# Duplication
obs <- obs[which(obs$Duplication == 0),]
simul <- simul[which(simul$Duplication == 0),]


#par(mfrow=c(1,2))
obs <- obs[which(obs$midist_obs < 3000000),]
simul <- simul[which(simul$midist_obs < 3000000),]
obs$class <-cut(obs$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
simul$class <-cut(simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
class_leg <- c("25-75Kb", "475-525Kb", "975Kb-1.02Mb", "1.47-1.52Mb", "1.97-2.02Mb","2.47-2.52Mb")

##### Nb fragment obs par classe de distance #### 
conserv = 0.4
obs_dist <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & obs$Total_ungapped > conserv),])/ nrow(obs[which(obs$class == x),]))*100))
obs_dist$int_start <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & obs$Total_ungapped > conserv),]), n=nrow(obs[which(obs$class == x),])+1, p=0.5)$conf.int[1])*100)
obs_dist$int_end <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & obs$Total_ungapped > conserv),]), n=nrow(obs[which(obs$class == x),])+1, p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & simul$Total_ungapped > conserv),])/ nrow(simul[which(simul$class == x),]))*100))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & simul$Total_ungapped > conserv),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & simul$Total_ungapped > conserv),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[2])*100)

plot(obs_dist$inter[1:50], type="b", col="red", cex =0.7, ylim=c(30,50), main="PECAN not dupli: Total seq > 0.4", ylab="Proportion fragment (%)", xlab="", xaxt = "n")
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10), par("usr")[3] - 1, labels = class_leg, srt = 30, pos = 1, xpd = TRUE, cex=0.8)

#legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n', cex=0.8)

obs_dist <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & obs$Allexon_ungapped > conserv),])/ nrow(obs[which(obs$class == x),]))*100))
obs_dist$int_start <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & obs$Allexon_ungapped > conserv),]), n=nrow(obs[which(obs$class == x),])+1, p=0.5)$conf.int[1])*100)
obs_dist$int_end <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & obs$Allexon_ungapped > conserv),]), n=nrow(obs[which(obs$class == x),])+1, p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & simul$Allexon_ungapped > conserv),])/ nrow(simul[which(simul$class == x),]))*100))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & simul$Allexon_ungapped > conserv),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & simul$Allexon_ungapped > conserv),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[2])*100)

plot(obs_dist$inter[1:50], type="b", col="red", cex =0.7, ylim=c(60,80), main="Human: Extract exon > 0.4", ylab="", xlab="", xaxt = "n")
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10), par("usr")[3] - 1, labels = class_leg, srt = 30, pos = 1, xpd = TRUE, cex=0.8)

#legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n', cex=0.8)






##### Score obs par classe de distance #### 
par(mfrow=c(1,2))

# All
obs_dist <- data.frame(inter = sapply(levels(obs$class), function(x) mean(obs[which(obs$class == x),]$Total_ungapped)))
obs_dist$int_start <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$Total_ungapped)[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$Total_ungapped)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$Total_ungapped)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$Total_ungapped)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$Total_ungapped)[["conf.int"]][2])

plot(obs_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(0.25,0.40), main="TBA : Total seq", xlab="",ylab="Mean score", xaxt = "n")
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10), par("usr")[3]-0.01, labels = class_leg, srt = 30, xpd = TRUE, cex=0.8)

# Extract all exons
obs_dist <- data.frame(inter = sapply(levels(obs$class), function(x) mean(obs[which(obs$class == x),]$Allexon_ungapped)))
obs_dist$int_start <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$Allexon_ungapped)[["conf.int"]][1])
obs_dist$int_end <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$Allexon_ungapped)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$Allexon_ungapped)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$Allexon_ungapped)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$Allexon_ungapped)[["conf.int"]][2])

plot(obs_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(0.25,0.40), main="TBA : Extract exon", xlab="",ylab="Mean score", xaxt = "n")
for (row in 1:nrow(obs_dist[1:50,])){
  segments(x0=row,y0=obs_dist[row,]$int_start,x1=row,y1=obs_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10), par("usr")[3]-0.01, labels = class_leg, srt = 30, xpd = TRUE, cex=0.8)

#legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

