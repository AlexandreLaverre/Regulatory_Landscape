###### MidDistance ~ Conservation #####
library(ggplot2)
library(gridExtra)
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/human2mouse/")
conserv <- read.table("PIR_cons_all_overlap.txt3", header=T)
simul <- read.table("PIR_cons_all_overlap_simul.txt3", header=T)

# Defining ratio of pb according to length of fragments
conserv$part_coding_exon <- (conserv$coding_exon_pb/conserv$length)*100
simul$part_coding_exon <- (simul$coding_exon_pb/simul$length)*100
conserv$part_nocoding_exon <- (conserv$nocoding_exon_pb/conserv$length)*100
simul$part_nocoding_exon <- (simul$nocoding_exon_pb/simul$length)*100
conserv$part_repeat <- (conserv$repeat_pb/conserv$length)*100
simul$part_repeat <- (simul$repeat_pb/simul$length)*100
conserv$part_phastcons <- (conserv$phastcons_noexonic250/conserv$length)*100
simul$part_phastcons <- (simul$phastcons_noexonic250/simul$length)*100
conserv$part_TSS <- (conserv$TSS_count/conserv$length)*100
simul$part_TSS <- (simul$TSS_count/simul$length)*100


# Corrélation exon ~ TSS
plot(log(conserv$coding_exon_pb+1)~conserv$TSS_count, cex=0.1)
abline(lm(log(conserv$coding_exon_pb+1)~conserv$TSS_count), col="red")

plot(conserv$part_coding_exon~conserv$part_TSS, cex=0.1)
abline(lm(conserv$part_coding_exon~conserv$part_TSS), col="red")

plot(conserv$TSS_count~conserv$nb_contact, cex=0.1)
abline(lm(conserv$TSS_count~conserv$nb_contact), col="red")


################################## Nombre moyen de peaks par classe de distance ################################
par(mfrow=c(2,3), mai = c(0.3, 0.7, 0.3, 0.1))
conserv <- conserv[which(conserv$midist_obs < 3000000),]
simul <- simul[which(simul$midist_obs < 3000000),]
conserv$class <-cut(conserv$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
simul$class <-cut(simul$midist_obs, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)
class_leg <- c("25Kb",  "1Mb",  "2Mb")

## TSS
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$part_TSS)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_TSS)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_TSS)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$part_TSS)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$part_TSS)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$part_TSS)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(0,0.02), main="TSS", xlab="",ylab="Mean fragment proportion (%)", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.0003, labels = class_leg, pos = 1, xpd = TRUE)

## coding_exon
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$part_coding_exon)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_coding_exon)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_coding_exon)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$part_coding_exon)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$part_coding_exon)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$part_coding_exon)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(1.5,7), main="Coding.gene exon", xlab="",ylab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.1, labels = class_leg, pos = 1, xpd = TRUE)
#legend("topright", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')


## nocoding_exon
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$part_nocoding_exon)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_nocoding_exon)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_nocoding_exon)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$part_nocoding_exon)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$part_nocoding_exon)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$part_nocoding_exon)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(0.3, 2.8),main="Non-coding.gene exon", xlab="",ylab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.03, labels = class_leg, pos = 1, xpd = TRUE)

## phastcons
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$part_phastcons)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_phastcons)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_phastcons)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$part_phastcons)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$part_phastcons)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$part_phastcons)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", cex=0.7, ylim=c(3, 6),main="Phastcons Element", xlab="",ylab="Mean fragment proportion (%)", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.03, labels = class_leg, pos = 1, xpd = TRUE)

## repeat
conserv_dist <- data.frame(inter = sapply(levels(conserv$class), function(x) mean(conserv[which(conserv$class == x),]$part_repeat)))
conserv_dist$int_start <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_repeat)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$class), function(x) t.test(conserv[which(conserv$class == x),]$part_repeat)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$part_repeat)))
simul_dist$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$part_repeat)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$part_repeat)[["conf.int"]][2])

plot(conserv_dist$inter[1:50], type="b", col="red", ylim=c(38, 55),main="Repeat Element", xlab="",ylab="", xaxt = "n", cex=0.7)
for (row in 1:nrow(conserv_dist[1:50,])){
  segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}

points(simul_dist$inter[1:50], type="b", col="blue", cex=0.7)
for (row in 1:nrow(simul_dist[1:50,])){
  segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,51,20), labels=F)
text(seq(1,51,20), par("usr")[3]-0.03, labels = class_leg, pos = 1, xpd = TRUE)


################### Defining class of proportion (%) ###################
conserv$part_coding_exon <-cut(conserv$part_coding_exon, breaks=seq(from=0, to=100, by=2), include.lowest = T, include.highest=T)
simul$part_coding_exon <-cut(simul$part_coding_exon, breaks=seq(from=0, to=100, by=2), include.lowest = T, include.highest=T)
conserv$part_nocoding_exon <-cut(conserv$part_nocoding_exon, breaks=seq(from=0, to=100, by=10), include.lowest = T, include.highest=T)
simul$part_nocoding_exon <-cut(simul$part_nocoding_exon, breaks=seq(from=0, to=100, by=10), include.lowest = T, include.highest=T)
conserv$part_repeat <-cut(conserv$part_repeat, breaks=seq(from=0, to=100, by=10), include.lowest = T, include.highest=T)
simul$part_repeat <-cut(simul$part_repeat, breaks=seq(from=0, to=100, by=10), include.lowest = T, include.highest=T)
conserv$part_phastcons <-cut(conserv$part_phastcons, breaks=seq(from=0, to=100, by=10), include.lowest = T, include.highest=T)
simul$part_phastcons <-cut(simul$part_phastcons, breaks=seq(from=0, to=100, by=10), include.lowest = T, include.highest=T)
leg <- c("0-10", "10-20", "20-30", "30-40", "40-50","50-60", "60-70", "70-80", "80-90", "90-100")

# Quelle est la proportion de fragments chevauchant avec les différents composants ?
par(mfrow=c(2,2), mai = c(0.3, 0.7, 0.3, 0.1))
# coding_exon
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_coding_exon), function(x) (nrow(conserv[which(conserv$part_coding_exon == x),])/nrow(conserv))*100))
conserv_dist$int_start <- sapply(levels(conserv$part_coding_exon), function(x) (prop.test(x = nrow(conserv[which(conserv$part_coding_exon == x),]), n=nrow(conserv), p=0.5)$conf.int[1])*100)
conserv_dist$int_end <- sapply(levels(conserv$part_coding_exon), function(x) (prop.test(x = nrow(conserv[which(conserv$part_coding_exon == x),]), n=nrow(conserv), p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(simul$part_coding_exon), function(x) (nrow(simul[which(simul$part_coding_exon == x),])/nrow(simul))*100))
simul_dist$int_start <- sapply(levels(simul$part_coding_exon), function(x) (prop.test(x = nrow(simul[which(simul$part_coding_exon == x),]), n=nrow(simul), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(simul$part_coding_exon), function(x) (prop.test(x = nrow(simul[which(simul$part_coding_exon == x),]), n=nrow(simul), p=0.5)$conf.int[2])*100)

plot(conserv_dist$inter, type="b", col="red",  main="Coding.gene exon",ylab="Fragments (%)",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-3, labels = leg, pos = 1, xpd = TRUE)
legend("topright", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

# nocoding_exon
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_nocoding_exon), function(x) (nrow(conserv[which(conserv$part_nocoding_exon == x),])/nrow(conserv))*100))
conserv_dist$int_start <- sapply(levels(conserv$part_nocoding_exon), function(x) (prop.test(x = nrow(conserv[which(conserv$part_nocoding_exon == x),]), n=nrow(conserv), p=0.5)$conf.int[1])*100)
conserv_dist$int_end <- sapply(levels(conserv$part_nocoding_exon), function(x) (prop.test(x = nrow(conserv[which(conserv$part_nocoding_exon == x),]), n=nrow(conserv), p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(simul$part_nocoding_exon), function(x) (nrow(simul[which(simul$part_nocoding_exon == x),])/nrow(simul))*100))
simul_dist$int_start <- sapply(levels(simul$part_nocoding_exon), function(x) (prop.test(x = nrow(simul[which(simul$part_nocoding_exon == x),]), n=nrow(simul), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(simul$part_nocoding_exon), function(x) (prop.test(x = nrow(simul[which(simul$part_nocoding_exon == x),]), n=nrow(simul), p=0.5)$conf.int[2])*100)

plot(conserv_dist$inter, type="b", col="red",  main="Nocoding exon",ylab="",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.2, labels = leg, pos = 1, xpd = TRUE)

# phastcons
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_phastcons), function(x) (nrow(conserv[which(conserv$part_phastcons == x),])/nrow(conserv))*100))
conserv_dist$int_start <- sapply(levels(conserv$part_phastcons), function(x) (prop.test(x = nrow(conserv[which(conserv$part_phastcons == x),]), n=nrow(conserv), p=0.5)$conf.int[1])*100)
conserv_dist$int_end <- sapply(levels(conserv$part_phastcons), function(x) (prop.test(x = nrow(conserv[which(conserv$part_phastcons == x),]), n=nrow(conserv), p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(simul$part_phastcons), function(x) (nrow(simul[which(simul$part_phastcons == x),])/nrow(simul))*100))
simul_dist$int_start <- sapply(levels(simul$part_phastcons), function(x) (prop.test(x = nrow(simul[which(simul$part_phastcons == x),]), n=nrow(simul), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(simul$part_phastcons), function(x) (prop.test(x = nrow(simul[which(simul$part_phastcons == x),]), n=nrow(simul), p=0.5)$conf.int[2])*100)

plot(conserv_dist$inter, type="b", col="red",  main="Phastcons Element",ylab="Fragments (%)",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.2, labels = leg, pos = 1, xpd = TRUE)

# repeat
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_repeat), function(x) (nrow(conserv[which(conserv$part_repeat == x),])/nrow(conserv))*100))
conserv_dist$int_start <- sapply(levels(conserv$part_repeat), function(x) (prop.test(x = nrow(conserv[which(conserv$part_repeat == x),]), n=nrow(conserv), p=0.5)$conf.int[1])*100)
conserv_dist$int_end <- sapply(levels(conserv$part_repeat), function(x) (prop.test(x = nrow(conserv[which(conserv$part_repeat == x),]), n=nrow(conserv), p=0.5)$conf.int[2])*100)

simul_dist <- data.frame(inter = sapply(levels(simul$part_repeat), function(x) (nrow(simul[which(simul$part_repeat == x),])/nrow(simul))*100))
simul_dist$int_start <- sapply(levels(simul$part_repeat), function(x) (prop.test(x = nrow(simul[which(simul$part_repeat == x),]), n=nrow(simul), p=0.5)$conf.int[1])*100)
simul_dist$int_end <- sapply(levels(simul$part_repeat), function(x) (prop.test(x = nrow(simul[which(simul$part_repeat == x),]), n=nrow(simul), p=0.5)$conf.int[2])*100)

plot(conserv_dist$inter, type="b", col="red",  ylim=c(2,20),main="Repeat Element",ylab="",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.2, labels = leg, pos = 1, xpd = TRUE)

####################################### Conserv by potential enhancer #############################################
# Quelle est la conservation des fragments chevauchant ces composants ?
par(mfrow=c(2,2), mai = c(0.3, 0.7, 0.3, 0.1))
#bas, gauche, haut, droite

## coding_exon
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_coding_exon), function(x) mean(conserv[which(conserv$part_coding_exon == x),]$PIR_score)))
conserv_dist$int_start <- sapply(levels(conserv$part_coding_exon), function(x) t.test(conserv[which(conserv$part_coding_exon == x),]$PIR_score)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$part_coding_exon), function(x) t.test(conserv[which(conserv$part_coding_exon == x),]$PIR_score)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$part_coding_exon), function(x) mean(simul[which(simul$part_coding_exon == x),]$PIR_score)))
simul_dist$int_start <- sapply(levels(simul$part_coding_exon), function(x) t.test(simul[which(simul$part_coding_exon == x),]$PIR_score)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$part_coding_exon), function(x) t.test(simul[which(simul$part_coding_exon == x),]$PIR_score)[["conf.int"]][2])

plot(conserv_dist$inter, type="b", col="red", ylim=c(0.15,0.6),main="Coding exon",ylab="Score conserv",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.006, labels = leg, pos = 1, xpd = TRUE)
legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

# nocoding_exon
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_nocoding_exon), function(x) mean(conserv[which(conserv$part_nocoding_exon == x),]$PIR_score)))
conserv_dist$int_start <- sapply(levels(conserv$part_nocoding_exon), function(x) t.test(conserv[which(conserv$part_nocoding_exon == x),]$PIR_score)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$part_nocoding_exon), function(x) t.test(conserv[which(conserv$part_nocoding_exon == x),]$PIR_score)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$part_nocoding_exon), function(x) mean(simul[which(simul$part_nocoding_exon == x),]$PIR_score)))
simul_dist$int_start <- sapply(levels(simul$part_nocoding_exon), function(x) t.test(simul[which(simul$part_nocoding_exon == x),]$PIR_score)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$part_nocoding_exon), function(x) t.test(simul[which(simul$part_nocoding_exon == x),]$PIR_score)[["conf.int"]][2])

plot(conserv_dist$inter, type="b", col="red", ylim=c(0.15,0.34), main="Non coding exon", ylab="",xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,21,2), labels=F)
text(seq(1,10,1), par("usr")[3]-0.006, labels = leg, pos = 1, xpd = TRUE)

## phastcons
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_phastcons), function(x) mean(conserv[which(conserv$part_phastcons == x),]$PIR_score)))
conserv_dist$int_start <- sapply(levels(conserv$part_phastcons), function(x) t.test(conserv[which(conserv$part_phastcons == x),]$PIR_score)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$part_phastcons), function(x) t.test(conserv[which(conserv$part_phastcons == x),]$PIR_score)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$part_phastcons), function(x) mean(simul[which(simul$part_phastcons == x),]$PIR_score)))
simul_dist$int_start <- sapply(levels(simul$part_phastcons), function(x) t.test(simul[which(simul$part_phastcons == x),]$PIR_score)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$part_phastcons), function(x) t.test(simul[which(simul$part_phastcons == x),]$PIR_score)[["conf.int"]][2])

plot(conserv_dist$inter, type="b", col="red", ylim=c(0.10,0.85), main="Phastcons Element", ylab="Score conserv",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,15,2), labels=F)
text(seq(1,10,1), par("usr")[3]-0.006, labels = leg, pos = 1, xpd = TRUE)

## repeat
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_repeat), function(x) mean(conserv[which(conserv$part_repeat == x),]$PIR_score)))
conserv_dist$int_start <- sapply(levels(conserv$part_repeat), function(x) t.test(conserv[which(conserv$part_repeat == x),]$PIR_score)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$part_repeat), function(x) t.test(conserv[which(conserv$part_repeat == x),]$PIR_score)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$part_repeat), function(x) mean(simul[which(simul$part_repeat == x),]$PIR_score)))
simul_dist$int_start <- sapply(levels(simul$part_repeat), function(x) t.test(simul[which(simul$part_repeat == x),]$PIR_score)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$part_repeat), function(x) t.test(simul[which(simul$part_repeat == x),]$PIR_score)[["conf.int"]][2])

plot(conserv_dist$inter, type="b", col="red", ylim=c(0,0.40), main="Repeat Element", ylab="",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,20,2), labels=F)
text(seq(1,10,1), par("usr")[3]-0.006, labels = leg, pos = 1, xpd = TRUE)









###############
## coding_exon
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_coding_exon), function(x) mean(conserv[which(conserv$part_coding_exon == x),]$PIR_score)))
conserv_dist$int_start <- sapply(levels(conserv$part_coding_exon), function(x) t.test(conserv[which(conserv$part_coding_exon == x),]$PIR_score)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$part_coding_exon), function(x) t.test(conserv[which(conserv$part_coding_exon == x),]$PIR_score)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$part_coding_exon), function(x) mean(simul[which(simul$part_coding_exon == x),]$PIR_score)))
simul_dist$int_start <- sapply(levels(simul$part_coding_exon), function(x) t.test(simul[which(simul$part_coding_exon == x),]$PIR_score)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$part_coding_exon), function(x) t.test(simul[which(simul$part_coding_exon == x),]$PIR_score)[["conf.int"]][2])

plot(conserv_dist$inter, type="b", col="red", ylim=c(0.15,0.6),main="Coding exon",ylab="Score conserv",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.006, labels = leg, pos = 1, xpd = TRUE)
legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

## coding_exon
conserv_dist <- data.frame(inter = sapply(levels(conserv$part_coding_exon), function(x) mean(conserv[which(conserv$part_coding_exon == x),]$length)))
conserv_dist$int_start <- sapply(levels(conserv$part_coding_exon), function(x) t.test(conserv[which(conserv$part_coding_exon == x),]$length)[["conf.int"]][1])
conserv_dist$int_end <- sapply(levels(conserv$part_coding_exon), function(x) t.test(conserv[which(conserv$part_coding_exon == x),]$length)[["conf.int"]][2])

simul_dist <- data.frame(inter = sapply(levels(simul$part_coding_exon), function(x) mean(simul[which(simul$part_coding_exon == x),]$length)))
simul_dist$int_start <- sapply(levels(simul$part_coding_exon), function(x) t.test(simul[which(simul$part_coding_exon == x),]$length)[["conf.int"]][1])
simul_dist$int_end <- sapply(levels(simul$part_coding_exon), function(x) t.test(simul[which(simul$part_coding_exon == x),]$length)[["conf.int"]][2])

plot(conserv_dist$inter, type="b", col="red", ylim=c(0,6000), main="Coding exon",ylab="Length (pb)",xlab="", xaxt = "n")
for (row in 1:nrow(conserv_dist[1:50,])){segments(x0=row,y0=conserv_dist[row,]$int_start,x1=row,y1=conserv_dist[row,]$int_end, col='red', lwd=0.3)}
points(simul_dist$inter[1:50], type="b", col="blue")
for (row in 1:nrow(simul_dist[1:50,])){segments(x0=row,y0=simul_dist[row,]$int_start,x1=row,y1=simul_dist[row,]$int_end, col='blue', lwd=0.3)}
axis(1, at=seq(1,10,1), labels=F)
text(seq(1,10,1), par("usr")[3]-0.006, labels = leg, pos = 1, xpd = TRUE)
legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

