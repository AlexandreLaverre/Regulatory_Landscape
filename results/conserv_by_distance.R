################################### Conservation de séquence par classe de distance ###################################
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/mouse2human/")
conserv <- read.table("PIR_cons_all_overlap.txt", header=T)
simul <- read.table("PIR_cons_all_overlap_simul.txt", header=T)

par(mfrow=c(1,2))
conserv <- conserv[which(conserv$midist_obs < 3000000),]
simul <- simul[which(simul$midist_obs < 3000000),]
conserv$class <-cut(conserv$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
simul$class <-cut(simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
class_leg <- c("25-75Kb", "475-525Kb", "975Kb-1.02Mb", "1.47-1.52Mb", "1.97-2.02Mb","2.47-2.52Mb")

##### Nb fragment conserv par classe de distance #### 
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) (nrow(conserv[which(conserv$class == x & conserv$PIR_score >0),])/ nrow(conserv[which(conserv$class == x),]))*100))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x)  (prop.test(x = nrow(conserv[which(conserv$class == x & conserv$PIR_score >0),]), n=nrow(conserv[which(conserv$class == x),])+1, p=0.5)$conf.int[1])*100)
conserv_dist$int_end <- sapply(levels(conserv$class), function(x)  (prop.test(x = nrow(conserv[which(conserv$class == x & conserv$PIR_score >0),]), n=nrow(conserv[which(conserv$class == x),])+1, p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & simul$PIR_score >0),])/ nrow(simul[which(simul$class == x),]))*100))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & simul$PIR_score >0),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & simul$PIR_score >0),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[2])*100)

plot(conserv_dist$inter[1:50], type="b", col="red", cex =0.7, ylim=c(42,70), main="Mouse2Human", ylab="Proportion fragment conserv (%)", xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10), par("usr")[3] - 1, labels = class_leg, srt = 30, pos = 1, xpd = TRUE, cex=0.8)

legend("topleft", cex=1.2, fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

##### Score conserv par classe de distance #### 
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$PIR_score)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$PIR_score)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$PIR_score)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$PIR_score)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$PIR_score)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$PIR_score)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(0.15,0.28), main="Mouse2Human", xlab="",ylab="Mean score conserv", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10), par("usr")[3]-0.01, labels = class_leg, srt = 60, xpd = TRUE, cex=0.8)

#legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

