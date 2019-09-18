###### MidDistance ~ Conservation #####
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/human2mouse/")
conserv <- read.table("PIR_cons_all_overlap.txt2", header=T)
simul <- read.table("PIR_cons_all_overlap_simul.txt2", header=T)

### Overlap potential enhancers ### 
# Quelle est la proportion de fragments chevauchant avec des potentiels enhancers ?
par(mfrow=c(2,2), mai = c(0.3, 0.7, 0.3, 0.1))
# CAGE
conserv_CAGE <- conserv[which(conserv$CAGE_count < 10 & conserv$CAGE_count > 0),]
conserv_CAGE$CAGE_count <- as.factor(conserv_CAGE$CAGE_count)
conserv_CAGE_dist <- data.frame(inter = sapply(levels(conserv_CAGE$CAGE_count), function(x) (nrow(conserv_CAGE[which(conserv_CAGE$CAGE_count == x),])/nrow(conserv))*100))
conserv_CAGE_dist$int_start <- sapply(levels(conserv_CAGE$CAGE_count), function(x) (prop.test(x = nrow(conserv_CAGE[which(conserv_CAGE$CAGE_count == x),]), n=nrow(conserv), p=0.5)$conf.int[1])*100)
conserv_CAGE_dist$int_end <- sapply(levels(conserv_CAGE$CAGE_count), function(x) (prop.test(x = nrow(conserv_CAGE[which(conserv_CAGE$CAGE_count == x),]), n=nrow(conserv), p=0.5)$conf.int[2])*100)

simul_CAGE <- simul[which(simul$CAGE_count < 10 & simul$CAGE_count > 0),]
simul_CAGE$CAGE_count <- as.factor(simul_CAGE$CAGE_count)
simul_CAGE_dist <- data.frame(inter = sapply(levels(simul_CAGE$CAGE_count), function(x) (nrow(simul_CAGE[which(simul_CAGE$CAGE_count == x),])/nrow(simul))*100))
simul_CAGE_dist$int_start <- sapply(levels(simul_CAGE$CAGE_count), function(x) (prop.test(x = nrow(simul_CAGE[which(simul_CAGE$CAGE_count == x),]), n=nrow(simul), p=0.5)$conf.int[1])*100)
simul_CAGE_dist$int_end <- sapply(levels(simul_CAGE$CAGE_count), function(x) (prop.test(x = nrow(simul_CAGE[which(simul_CAGE$CAGE_count == x),]), n=nrow(simul), p=0.5)$conf.int[2])*100)

plot(conserv_CAGE_dist$inter, type="b", col="red",  main="CAGE peak",ylab="Fragments (%)",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_CAGE_dist[1:50,])){segments(x0=row,y0=conserv_CAGE_dist[row,]$int_start,x1=row,y1=conserv_CAGE_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_CAGE_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_CAGE_dist[1:50,])){segments(x0=row,y0=simul_CAGE_dist[row,]$int_start,x1=row,y1=simul_CAGE_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.2, labels = seq(1,10,1), pos = 1, xpd = TRUE)
legend("topright", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

# ENCODE
conserv_ENCODE <- conserv[which(conserv$ENCODE_count < 10 & conserv$ENCODE_count > 0),]
conserv_ENCODE$ENCODE_count <- as.factor(conserv_ENCODE$ENCODE_count)
conserv_ENCODE_dist <- data.frame(inter = sapply(levels(conserv_ENCODE$ENCODE_count), function(x) (nrow(conserv_ENCODE[which(conserv_ENCODE$ENCODE_count == x),])/nrow(conserv))*100))
conserv_ENCODE_dist$int_start <- sapply(levels(conserv_ENCODE$ENCODE_count), function(x) (prop.test(x = nrow(conserv_ENCODE[which(conserv_ENCODE$ENCODE_count == x),]), n=nrow(conserv), p=0.5)$conf.int[1])*100)
conserv_ENCODE_dist$int_end <- sapply(levels(conserv_ENCODE$ENCODE_count), function(x) (prop.test(x = nrow(conserv_ENCODE[which(conserv_ENCODE$ENCODE_count == x),]), n=nrow(conserv), p=0.5)$conf.int[2])*100)

simul_ENCODE <- simul[which(simul$ENCODE_count < 10 & simul$ENCODE_count > 0),]
simul_ENCODE$ENCODE_count <- as.factor(simul_ENCODE$ENCODE_count)
simul_ENCODE_dist <- data.frame(inter = sapply(levels(simul_ENCODE$ENCODE_count), function(x) (nrow(simul_ENCODE[which(simul_ENCODE$ENCODE_count == x),])/nrow(simul))*100))
simul_ENCODE_dist$int_start <- sapply(levels(simul_ENCODE$ENCODE_count), function(x) (prop.test(x = nrow(simul_ENCODE[which(simul_ENCODE$ENCODE_count == x),]), n=nrow(simul), p=0.5)$conf.int[1])*100)
simul_ENCODE_dist$int_end <- sapply(levels(simul_ENCODE$ENCODE_count), function(x) (prop.test(x = nrow(simul_ENCODE[which(simul_ENCODE$ENCODE_count == x),]), n=nrow(simul), p=0.5)$conf.int[2])*100)

plot(conserv_ENCODE_dist$inter, type="b", col="red",  main="ENCODE peak",ylab="",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_ENCODE_dist[1:50,])){segments(x0=row,y0=conserv_ENCODE_dist[row,]$int_start,x1=row,y1=conserv_ENCODE_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_ENCODE_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_ENCODE_dist[1:50,])){segments(x0=row,y0=simul_ENCODE_dist[row,]$int_start,x1=row,y1=simul_ENCODE_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.2, labels = seq(1,10,1), pos = 1, xpd = TRUE)

# GRO_seq
conserv_GRO_seq <- conserv[which(conserv$GRO_seq_count < 10 & conserv$GRO_seq_count > 0),]
conserv_GRO_seq$GRO_seq_count <- as.factor(conserv_GRO_seq$GRO_seq_count)
conserv_GRO_seq_dist <- data.frame(inter = sapply(levels(conserv_GRO_seq$GRO_seq_count), function(x) (nrow(conserv_GRO_seq[which(conserv_GRO_seq$GRO_seq_count == x),])/nrow(conserv))*100))
conserv_GRO_seq_dist$int_start <- sapply(levels(conserv_GRO_seq$GRO_seq_count), function(x) (prop.test(x = nrow(conserv_GRO_seq[which(conserv_GRO_seq$GRO_seq_count == x),]), n=nrow(conserv), p=0.5)$conf.int[1])*100)
conserv_GRO_seq_dist$int_end <- sapply(levels(conserv_GRO_seq$GRO_seq_count), function(x) (prop.test(x = nrow(conserv_GRO_seq[which(conserv_GRO_seq$GRO_seq_count == x),]), n=nrow(conserv), p=0.5)$conf.int[2])*100)

simul_GRO_seq <- simul[which(simul$GRO_seq_count < 10 & simul$GRO_seq_count > 0),]
simul_GRO_seq$GRO_seq_count <- as.factor(simul_GRO_seq$GRO_seq_count)
simul_GRO_seq_dist <- data.frame(inter = sapply(levels(simul_GRO_seq$GRO_seq_count), function(x) (nrow(simul_GRO_seq[which(simul_GRO_seq$GRO_seq_count == x),])/nrow(simul))*100))
simul_GRO_seq_dist$int_start <- sapply(levels(simul_GRO_seq$GRO_seq_count), function(x) (prop.test(x = nrow(simul_GRO_seq[which(simul_GRO_seq$GRO_seq_count == x),]), n=nrow(simul), p=0.5)$conf.int[1])*100)
simul_GRO_seq_dist$int_end <- sapply(levels(simul_GRO_seq$GRO_seq_count), function(x) (prop.test(x = nrow(simul_GRO_seq[which(simul_GRO_seq$GRO_seq_count == x),]), n=nrow(simul), p=0.5)$conf.int[2])*100)

plot(conserv_GRO_seq_dist$inter, type="b", col="red",  main="GRO_seq peak",ylab="Fragments (%)",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_GRO_seq_dist[1:50,])){segments(x0=row,y0=conserv_GRO_seq_dist[row,]$int_start,x1=row,y1=conserv_GRO_seq_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_GRO_seq_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_GRO_seq_dist[1:50,])){segments(x0=row,y0=simul_GRO_seq_dist[row,]$int_start,x1=row,y1=simul_GRO_seq_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.2, labels = seq(1,10,1), pos = 1, xpd = TRUE)

# RoadMap
conserv_RoadMap <- conserv[which(conserv$RoadMap_count < 10 & conserv$RoadMap_count > 0),]
conserv_RoadMap$RoadMap_count <- as.factor(conserv_RoadMap$RoadMap_count)
conserv_RoadMap_dist <- data.frame(inter = sapply(levels(conserv_RoadMap$RoadMap_count), function(x) (nrow(conserv_RoadMap[which(conserv_RoadMap$RoadMap_count == x),])/nrow(conserv))*100))
conserv_RoadMap_dist$int_start <- sapply(levels(conserv_RoadMap$RoadMap_count), function(x) (prop.test(x = nrow(conserv_RoadMap[which(conserv_RoadMap$RoadMap_count == x),]), n=nrow(conserv), p=0.5)$conf.int[1])*100)
conserv_RoadMap_dist$int_end <- sapply(levels(conserv_RoadMap$RoadMap_count), function(x) (prop.test(x = nrow(conserv_RoadMap[which(conserv_RoadMap$RoadMap_count == x),]), n=nrow(conserv), p=0.5)$conf.int[2])*100)

simul_RoadMap <- simul[which(simul$RoadMap_count < 10 & simul$RoadMap_count > 0),]
simul_RoadMap$RoadMap_count <- as.factor(simul_RoadMap$RoadMap_count)
simul_RoadMap_dist <- data.frame(inter = sapply(levels(simul_RoadMap$RoadMap_count), function(x) (nrow(simul_RoadMap[which(simul_RoadMap$RoadMap_count == x),])/nrow(simul))*100))
simul_RoadMap_dist$int_start <- sapply(levels(simul_RoadMap$RoadMap_count), function(x) (prop.test(x = nrow(simul_RoadMap[which(simul_RoadMap$RoadMap_count == x),]), n=nrow(simul), p=0.5)$conf.int[1])*100)
simul_RoadMap_dist$int_end <- sapply(levels(simul_RoadMap$RoadMap_count), function(x) (prop.test(x = nrow(simul_RoadMap[which(simul_RoadMap$RoadMap_count == x),]), n=nrow(simul), p=0.5)$conf.int[2])*100)

plot(conserv_RoadMap_dist$inter, type="b", col="red",  main="RoadMap peak",ylab="",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_RoadMap_dist[1:50,])){segments(x0=row,y0=conserv_RoadMap_dist[row,]$int_start,x1=row,y1=conserv_RoadMap_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_RoadMap_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_RoadMap_dist[1:50,])){segments(x0=row,y0=simul_RoadMap_dist[row,]$int_start,x1=row,y1=simul_RoadMap_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.2, labels = seq(1,10,1), pos = 1, xpd = TRUE)

####################################### Conserv by potential enhancer #############################################
# Quelle est la conservation des fragments chevauchant ces potentiels enhancers ?
par(mfrow=c(2,2), mai = c(0.3, 0.7, 0.3, 0.1))
#bas, gauche, haut, droite

## CAGE
conserv_CAGE <- conserv[which(conserv$CAGE_count < 10),]
conserv_CAGE$CAGE_count <- as.factor(conserv_CAGE$CAGE_count)
conserv_CAGE_dist <- data.frame(inter = sapply(levels(conserv_CAGE$CAGE_count), function(x) mean(conserv_CAGE[which(conserv_CAGE$CAGE_count == x),]$Score_allexon_ungapped)))
conserv_CAGE_dist$int_start <- sapply(levels(conserv_CAGE$CAGE_count), function(x) t.test(conserv_CAGE[which(conserv_CAGE$CAGE_count == x),]$Score_allexon_ungapped)[["conf.int"]][1])
conserv_CAGE_dist$int_end <- sapply(levels(conserv_CAGE$CAGE_count), function(x) t.test(conserv_CAGE[which(conserv_CAGE$CAGE_count == x),]$Score_allexon_ungapped)[["conf.int"]][2])

simul_CAGE <- simul[which(simul$CAGE_count < 10),]
simul_CAGE$CAGE_count <- as.factor(simul_CAGE$CAGE_count)
simul_CAGE_dist <- data.frame(inter = sapply(levels(simul_CAGE$CAGE_count), function(x) mean(simul_CAGE[which(simul_CAGE$CAGE_count == x),]$Score_allexon_ungapped)))
simul_CAGE_dist$int_start <- sapply(levels(simul_CAGE$CAGE_count), function(x) t.test(simul_CAGE[which(simul_CAGE$CAGE_count == x),]$Score_allexon_ungapped)[["conf.int"]][1])
simul_CAGE_dist$int_end <- sapply(levels(simul_CAGE$CAGE_count), function(x) t.test(simul_CAGE[which(simul_CAGE$CAGE_count == x),]$Score_allexon_ungapped)[["conf.int"]][2])

plot(conserv_CAGE_dist$inter, type="b", col="red", ylim=c(0.2,0.6),main="CAGE peak",ylab="Score conserv",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_CAGE_dist[1:50,])){segments(x0=row,y0=conserv_CAGE_dist[row,]$int_start,x1=row,y1=conserv_CAGE_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_CAGE_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_CAGE_dist[1:50,])){segments(x0=row,y0=simul_CAGE_dist[row,]$int_start,x1=row,y1=simul_CAGE_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.006, labels = seq(0,9,1), pos = 1, xpd = TRUE)
#legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

# ENCODE
conserv_ENCODE <- conserv[which(conserv$ENCODE_count < 21),]
conserv_ENCODE$ENCODE_count <- as.factor(conserv_ENCODE$ENCODE_count)
conserv_ENCODE_dist <- data.frame(inter = sapply(levels(conserv_ENCODE$ENCODE_count), function(x) mean(conserv_ENCODE[which(conserv_ENCODE$ENCODE_count == x),]$Score_allexon_ungapped)))
conserv_ENCODE_dist$int_start <- sapply(levels(conserv_ENCODE$ENCODE_count), function(x) t.test(conserv_ENCODE[which(conserv_ENCODE$ENCODE_count == x),]$Score_allexon_ungapped)[["conf.int"]][1])
conserv_ENCODE_dist$int_end <- sapply(levels(conserv_ENCODE$ENCODE_count), function(x) t.test(conserv_ENCODE[which(conserv_ENCODE$ENCODE_count == x),]$Score_allexon_ungapped)[["conf.int"]][2])

simul_ENCODE <- simul[which(simul$ENCODE_count < 21),]
simul_ENCODE$ENCODE_count <- as.factor(simul_ENCODE$ENCODE_count)
simul_ENCODE_dist <- data.frame(inter = sapply(levels(simul_ENCODE$ENCODE_count), function(x) mean(simul_ENCODE[which(simul_ENCODE$ENCODE_count == x),]$Score_allexon_ungapped)))
simul_ENCODE_dist$int_start <- sapply(levels(simul_ENCODE$ENCODE_count), function(x) t.test(simul_ENCODE[which(simul_ENCODE$ENCODE_count == x),]$Score_allexon_ungapped)[["conf.int"]][1])
simul_ENCODE_dist$int_end <- sapply(levels(simul_ENCODE$ENCODE_count), function(x) t.test(simul_ENCODE[which(simul_ENCODE$ENCODE_count == x),]$Score_allexon_ungapped)[["conf.int"]][2])

plot(conserv_ENCODE_dist$inter, type="b", col="red", ylim=c(0.2,0.6), main="ENCODE peak", ylab="",xaxt = "n")
for (row in 1:nrow(conserv_ENCODE_dist[1:50,])){segments(x0=row,y0=conserv_ENCODE_dist[row,]$int_start,x1=row,y1=conserv_ENCODE_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_ENCODE_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_ENCODE_dist[1:50,])){segments(x0=row,y0=simul_ENCODE_dist[row,]$int_start,x1=row,y1=simul_ENCODE_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,21,2), labels=F)
text(seq(1,21,2), par("usr")[3]-0.006, labels = seq(0,20,2), pos = 1, xpd = TRUE)

## GRO_seq
conserv_GRO <- conserv[which(conserv$GRO_seq_count < 15),]
conserv_GRO$GRO_seq_count <- as.factor(conserv_GRO$GRO_seq_count)
conserv_GRO_dist <- data.frame(inter = sapply(levels(conserv_GRO$GRO_seq_count), function(x) mean(conserv_GRO[which(conserv_GRO$GRO_seq_count == x),]$Score_allexon_ungapped)))
conserv_GRO_dist$int_start <- sapply(levels(conserv_GRO$GRO_seq_count), function(x) t.test(conserv_GRO[which(conserv_GRO$GRO_seq_count == x),]$Score_allexon_ungapped)[["conf.int"]][1])
conserv_GRO_dist$int_end <- sapply(levels(conserv_GRO$GRO_seq_count), function(x) t.test(conserv_GRO[which(conserv_GRO$GRO_seq_count == x),]$Score_allexon_ungapped)[["conf.int"]][2])

simul_GRO <- simul[which(simul$GRO_seq_count < 15),]
simul_GRO$GRO_seq_count <- as.factor(simul_GRO$GRO_seq_count)
simul_GRO_dist <- data.frame(inter = sapply(levels(simul_GRO$GRO_seq_count), function(x) mean(simul_GRO[which(simul_GRO$GRO_seq_count == x),]$Score_allexon_ungapped)))
simul_GRO_dist$int_start <- sapply(levels(simul_GRO$GRO_seq_count), function(x) t.test(simul_GRO[which(simul_GRO$GRO_seq_count == x),]$Score_allexon_ungapped)[["conf.int"]][1])
simul_GRO_dist$int_end <- sapply(levels(simul_GRO$GRO_seq_count), function(x) t.test(simul_GRO[which(simul_GRO$GRO_seq_count == x),]$Score_allexon_ungapped)[["conf.int"]][2])

plot(conserv_GRO_dist$inter, type="b", col="red", ylim=c(0.1,0.5), main="GRO_seq peak", ylab="Score conserv",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_GRO_dist[1:50,])){segments(x0=row,y0=conserv_GRO_dist[row,]$int_start,x1=row,y1=conserv_GRO_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_GRO_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_GRO_dist[1:50,])){segments(x0=row,y0=simul_GRO_dist[row,]$int_start,x1=row,y1=simul_GRO_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,15,2), labels=F)
text(seq(1,15,2), par("usr")[3]-0.006, labels = seq(0,14,2), pos = 1, xpd = TRUE)

## RoadMap
conserv_RoadMap <- conserv[which(conserv$RoadMap_count < 18),]
conserv_RoadMap$RoadMap_count <- as.factor(conserv_RoadMap$RoadMap_count)
conserv_RoadMap_dist <- data.frame(inter = sapply(levels(conserv_RoadMap$RoadMap_count), function(x) mean(conserv_RoadMap[which(conserv_RoadMap$RoadMap_count == x),]$Score_allexon_ungapped)))
conserv_RoadMap_dist$int_start <- sapply(levels(conserv_RoadMap$RoadMap_count), function(x) t.test(conserv_RoadMap[which(conserv_RoadMap$RoadMap_count == x),]$Score_allexon_ungapped)[["conf.int"]][1])
conserv_RoadMap_dist$int_end <- sapply(levels(conserv_RoadMap$RoadMap_count), function(x) t.test(conserv_RoadMap[which(conserv_RoadMap$RoadMap_count == x),]$Score_allexon_ungapped)[["conf.int"]][2])

simul_RoadMap <- simul[which(simul$RoadMap_count < 18),]
simul_RoadMap$RoadMap_count <- as.factor(simul_RoadMap$RoadMap_count)
simul_RoadMap_dist <- data.frame(inter = sapply(levels(simul_RoadMap$RoadMap_count), function(x) mean(simul_RoadMap[which(simul_RoadMap$RoadMap_count == x),]$Score_allexon_ungapped)))
simul_RoadMap_dist$int_start <- sapply(levels(simul_RoadMap$RoadMap_count), function(x) t.test(simul_RoadMap[which(simul_RoadMap$RoadMap_count == x),]$Score_allexon_ungapped)[["conf.int"]][1])
simul_RoadMap_dist$int_end <- sapply(levels(simul_RoadMap$RoadMap_count), function(x) t.test(simul_RoadMap[which(simul_RoadMap$RoadMap_count == x),]$Score_allexon_ungapped)[["conf.int"]][2])

plot(conserv_RoadMap_dist$inter, type="b", col="red", ylim=c(0.3,0.5), main="RoadMap peak", ylab="",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_RoadMap_dist[1:50,])){segments(x0=row,y0=conserv_RoadMap_dist[row,]$int_start,x1=row,y1=conserv_RoadMap_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_RoadMap_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_RoadMap_dist[1:50,])){segments(x0=row,y0=simul_RoadMap_dist[row,]$int_start,x1=row,y1=simul_RoadMap_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,20,2), labels=F)
text(seq(1,18,2), par("usr")[3]-0.006, labels = seq(0,17,2), pos = 1, xpd = TRUE)

################################## Nombre moyen de peaks par classe de distance ################################
par(mfrow=c(2,2), mai = c(0.3, 0.7, 0.3, 0.1))
conserv <- conserv[which(conserv$midist_obs < 3000000),]
simul <- simul[which(simul$midist_obs < 3000000),]
conserv$class <-cut(conserv$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
simul$class <-cut(simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
class_leg <- c("25Kb", "1Mb", "2Mb")

## CAGE
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$CAGE_count)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$CAGE_count)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$CAGE_count)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$CAGE_count)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$CAGE_count)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$CAGE_count)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7,  ylim=c(0.03,0.25), main="CAGE peaks", xlab="",ylab="Mean nb peaks", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.001, labels = class_leg, pos = 1, xpd = TRUE)
#legend("topright", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

## ENCODE
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$ENCODE_count)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$ENCODE_count)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$ENCODE_count)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$ENCODE_count)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$ENCODE_count)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$ENCODE_count)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(0.35, 1.4),main="ENCODE peaks", xlab="",ylab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.03, labels = class_leg, pos = 1, xpd = TRUE)

## GRO_seq
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$GRO_seq_count)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$GRO_seq_count)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$GRO_seq_count)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$GRO_seq_count)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$GRO_seq_count)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$GRO_seq_count)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(0.20, 0.95),main="GRO_seq peaks", xlab="",ylab="Mean nb peaks", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.02, labels = class_leg, pos = 1, xpd = TRUE)

## RoadMap
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$RoadMap_count)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$RoadMap_count)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$RoadMap_count)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$RoadMap_count)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$RoadMap_count)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$RoadMap_count)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(0.35, 1.65),main="RoadMap peaks", xlab="",ylab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.03, labels = class_leg, pos = 1, xpd = TRUE)

