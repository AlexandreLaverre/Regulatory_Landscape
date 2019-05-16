setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")
sp_origin = 'mouse'
sp_target = 'human'

obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction.txt2", sep=""), header=T)
simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_interaction_simul.txt2", sep=""), header=T)
 
obs$conserv_class <-cut(obs$PIR_score, breaks=seq(from=0, to=1, by=0.1), include.lowest = T)
simul$conserv_class <-cut(simul$PIR_score, breaks=seq(from=0, to=1, by=0.1), include.lowest = T)
color_pallete_function <- colorRampPalette(colors = c("blue", "orange", "red"), space = "Lab")
conserv_color <- color_pallete_function(nlevels(obs$conserv_class))

#### Filtres ###
# Length
obs <- obs[which(obs$bait_length > 250 & obs$bait_length < 20000 & obs$PIR_length >250 & obs$PIR_length < 20000),]
simul <- simul[which(simul$bait_length > 250 & simul$bait_length < 20000 & simul$PIR_length >250 & simul$PIR_length < 20000),]

# Duplication
obs <- obs[which(obs$bait_dupli == 0 & obs$PIR_dupli == 0),]
simul <- simul[which(simul$bait_dupli == 0 & simul$PIR_dupli == 0),]

# Overlap > 1 potential enh
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/human2mouse/")
obs_enh <- read.table("PIR_cons_all_overlap.txt3", header=T)
simul_enh <- read.table("PIR_cons_all_overlap_simul.txt3", header=T)

obs$PIR <- as.factor(sub(".*-", "", obs$origin_interaction))
simul$PIR <- as.factor(sub(".*-", "", simul$origin_interaction))

obs$CAGE <- obs$PIR %in%  obs_enh[which(obs_enh$CAGE_count > 0),]$PIR
obs$ENCODE <- obs$PIR %in%  obs_enh[which(obs_enh$ENCODE_count > 0),]$PIR
obs$GRO_seq <- obs$PIR %in%  obs_enh[which(obs_enh$GRO_seq_count > 0),]$PIR
obs$RoadMap <- obs$PIR %in%  obs_enh[which(obs_enh$RoadMap_count > 0),]$PIR
simul$CAGE <- simul$PIR %in%  simul_enh[which(simul_enh$CAGE_count > 0),]$PIR
simul$ENCODE <- simul$PIR %in%  simul_enh[which(simul_enh$ENCODE_count > 0),]$PIR
simul$GRO_seq <- simul$PIR %in%  simul_enh[which(simul_enh$GRO_seq_count > 0),]$PIR
simul$RoadMap <- simul$PIR %in%  simul_enh[which(simul_enh$RoadMap_count > 0),]$PIR


total <- nrow(obs[which(!is.na(obs$target_dist)),])/ nrow(obs)
CAGE <- nrow(obs[which(!is.na(obs$target_dist) & obs$CAGE == TRUE),])/ nrow(obs[which(obs$CAGE == TRUE),])
ENCODE <- nrow(obs[which(!is.na(obs$target_dist) & obs$ENCODE == TRUE),])/ nrow(obs[which(obs$ENCODE == TRUE),])
GRO_seq <- nrow(obs[which(!is.na(obs$target_dist) & obs$GRO_seq == TRUE),])/ nrow(obs[which(obs$GRO_seq == TRUE),])
RoadMap <- nrow(obs[which(!is.na(obs$target_dist) & obs$RoadMap == TRUE),])/ nrow(obs[which(obs$RoadMap == TRUE),])

total2 <- nrow(simul[which(!is.na(simul$target_dist)),])/ nrow(simul)
CAGE2 <- nrow(simul[which(!is.na(simul$target_dist) & simul$CAGE == TRUE),])/ nrow(simul[which(simul$CAGE == TRUE),])
ENCODE2 <- nrow(simul[which(!is.na(simul$target_dist) & simul$ENCODE == TRUE),])/ nrow(simul[which(simul$ENCODE == TRUE),])
GRO_seq2 <- nrow(simul[which(!is.na(simul$target_dist) & simul$GRO_seq == TRUE),])/ nrow(simul[which(simul$GRO_seq == TRUE),])
RoadMap2 <- nrow(simul[which(!is.na(simul$target_dist) & simul$RoadMap == TRUE),])/ nrow(simul[which(simul$RoadMap == TRUE),])

test <- c(total, total2, 0, CAGE, CAGE2, 0, ENCODE, ENCODE2, 0, GRO_seq, GRO_seq2, 0, RoadMap, RoadMap2)
col <- rep(c("red", "blue", rgb(1,1,1)),4)

barplot(test, border=col, col=rep(rgb(1,1,1),4), ylim=c(0,0.28), ylab="Proportion interaction conserv (%)")
leg <- c("Total", "CAGE", "ENCODE", "GRO_seq", "RoadMap")
text(c(1.2,5,8.5,12.2,15.8), par("usr")[3]-0.005, labels = leg, pos = 1, xpd = TRUE)
legend("topright", fill=NULL, border=col, legend = c("Paysage observé", "Paysage simulé"), bty='n')


############## Interaction conserv ~ distance ############## 
#obs <- obs[which(obs$CAGE == TRUE),]
#simul <- simul[which(simul$CAGE == TRUE),]

obs$class <-cut(obs$origin_dist, breaks=seq(from=25000, to=8000000, by=50000), include.lowest = T)
simul$class <- cut(simul$origin_dist, breaks=seq(from=25000, to=8000000, by=50000), include.lowest = T)

class_leg <- c("25-75Kb", "1.97-2.02Mb", "3.97-4.02Mb","5.97-6.02Mb","7.97-8.02Mb") #,"9.97-10.02Mb")

obs_conserv <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),])/ nrow(obs[which(obs$class == x),]))*100))
obs_conserv$int_start <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x),])+1, p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$class == x),])+1, p=0.5)$conf.int[2])*100)

obs_conserv <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & !is.na(obs$target_dist) & obs$origin_nb_tissu == 1),])/ nrow(obs[which(obs$class == x & obs$origin_nb_tissu == 1),]))*100))
obs_conserv$int_start <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist) & obs$origin_nb_tissu == 1),]), n=nrow(obs[which(obs$class == x & obs$origin_nb_tissu == 1),])+1, p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist) & obs$origin_nb_tissu == 1),]), n=nrow(obs[which(obs$class == x & obs$origin_nb_tissu == 1),])+1, p=0.5)$conf.int[2])*100)

obs_conserv_consti <- data.frame(inter = sapply(levels(obs$class), function(x) (nrow(obs[which(obs$class == x & !is.na(obs$target_dist) & obs$origin_nb_tissu > 1),])/ nrow(obs[which(obs$class == x & obs$origin_nb_tissu > 1),]))*100))
obs_conserv_consti$int_start <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist) & obs$origin_nb_tissu > 1),]), n=nrow(obs[which(obs$class == x & obs$origin_nb_tissu > 1),])+1, p=0.5)$conf.int[1])*100)
obs_conserv_consti$int_end <- sapply(levels(obs$class), function(x)  (prop.test(x = nrow(obs[which(obs$class == x & !is.na(obs$target_dist) & obs$origin_nb_tissu > 1),]), n=nrow(obs[which(obs$class == x & obs$origin_nb_tissu > 1),])+1, p=0.5)$conf.int[2])*100)

simul_conserv <- data.frame(inter = sapply(levels(simul$class), function(x) (nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),])/ nrow(simul[which(simul$class == x),]))*100))
simul_conserv$int_start <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[1])*100)
simul_conserv$int_end <- sapply(levels(simul$class), function(x)  (prop.test(x = nrow(simul[which(simul$class == x & !is.na(simul$target_dist)),]), n=nrow(simul[which(simul$class == x),]), p=0.5)$conf.int[2])*100)

obs_conserv$nb <- sapply(levels(obs$class), function(x) nrow(obs[which(obs$class == x),]))

plot(obs_conserv$inter, type="b", col="red", cex=0.5, main=paste(sp_origin,"interactions conserved in", sp_target),
     ylab="Proportion interaction conserv (%)", xlab="", xaxt = "n", ylim=c(0,30))
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}

points(obs_conserv_consti$inter, type="b", col="forestgreen", cex=0.5)
for (row in 1:nrow(obs_conserv_consti)){
  segments(x0=row,y0=obs_conserv_consti[row,]$int_start,x1=row,y1=obs_conserv_consti[row,]$int_end, col='forestgreen', lwd=0.3)}

points(simul_conserv$inter, type="b", col="blue", cex=0.5)
for (row in 1:nrow(simul_conserv)){
  segments(x0=row,y0=simul_conserv[row,]$int_start,x1=row,y1=simul_conserv[row,]$int_end, col='blue', lwd=0.3)}

#polygon(x=c(seq(1,199,1),seq(199,1,-1)), y=c(obs_conserv$int_start, rev(obs_conserv$int_end)), col=rgb(1,0,0,0.1), border = NA)
#polygon(x=c(seq(1,199,1),seq(199,1,-1)), y=c(simul_conserv$int_start, rev(simul_conserv$int_end)), col=rgb(0,0,1,0.1), border = NA)

axis(1, at=seq(1,201,40), labels=F)
text(seq(1,201,40), par("usr")[3]-0.3, labels = class_leg, pos = 1, xpd = TRUE, cex=0.8)
legend("topright", fill=c("forestgreen", "red","blue"), legend = c("> 1 cell type", "Cell specific", "Simulated"), bty='n', cex=0.8)

############## Score conserv ~ distance ############## 
obs$class <-cut(obs$origin_dist, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
simul$class <-cut(simul$origin_dist, breaks=seq(from=25000, to=10000000, by=50000), include.lowest = T)
class_leg <- c("25-75Kb", "1.97-2.02Mb", "3.97-4.02Mb","5.97-6.02Mb","7.97-8.02Mb","9.97-10.02Mb")

obs_conserv <- data.frame(inter = sapply(levels(obs$class), function(x) mean(obs[which(obs$class == x),]$PIR_score)))
obs_conserv$int_start <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$PIR_score)[["conf.int"]][1])
obs_conserv$int_end <- sapply(levels(obs$class), function(x) t.test(obs[which(obs$class == x),]$PIR_score)[["conf.int"]][2])

simul_conserv <- data.frame(inter = sapply(levels(simul$class), function(x) mean(simul[which(simul$class == x),]$PIR_score)))
simul_conserv$int_start <- sapply(levels(simul$class), function(x)  t.test(simul[which(simul$class == x),]$PIR_score)[["conf.int"]][1])
simul_conserv$int_end <- sapply(levels(simul$class), function(x) t.test(simul[which(simul$class == x),]$PIR_score)[["conf.int"]][2])

plot(obs_conserv$inter, type="b", col="red", cex=0.5, ylim=c(0.1,0.45),
     main=paste(sp_origin, "contacted regions conserved in", sp_target), ylab="Score conserv", xlab="", xaxt = "n")
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}
points(simul_conserv$inter, type="b", col="blue", cex=0.5)
for (row in 1:nrow(simul_conserv)){
  segments(x0=row,y0=simul_conserv[row,]$int_start,x1=row,y1=simul_conserv[row,]$int_end, col='blue', lwd=0.3)}

axis(1, at=seq(1,201,40), labels=F)
text(seq(1,201,40), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE, cex=0.8)
legend("topleft", fill=c("red","blue"), legend = c("Paysage observé", "Paysage simulé"), bty='n')

### Relation strength ~ frequency
boxplot(obs$origin_strength~obs$origin_nb_tissu, outline=F, notch=T, 
        xlab="Nb cell type", ylab="Contact strength (median)", main=paste(sp_origin, "interactions"))

# Conserv ~ nb cell type
par(mfrow=c(1,1))
barplot(table(obs$origin_nb_tissu), xlab="Nombre d'échantillon", ylab="Fréquence",main=paste(sp_origin, "interactions"))
obs_conserv <- data.frame(inter = sapply(levels(as.factor(obs$origin_nb_tissu)), function(x) (nrow(obs[which(obs$origin_nb_tissu == x & !is.na(obs$target_dist)),])/ nrow(obs[which(obs$origin_nb_tissu == x),]))*100))
obs_conserv$int_start <- sapply(levels(as.factor(obs$origin_nb_tissu)), function(x)  (prop.test(x = nrow(obs[which(obs$origin_nb_tissu == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$origin_nb_tissu == x),]), p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(as.factor(obs$origin_nb_tissu)), function(x)  (prop.test(x = nrow(obs[which(obs$origin_nb_tissu == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$origin_nb_tissu == x),]), p=0.5)$conf.int[2])*100)

plot(obs_conserv$inter, type="b", col="red", main=paste(sp_origin,"interactions conserved in", sp_target),
     ylab="Proportion interaction conserv (%)", xlab="Nombre d'échantillon")
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}

# Conserv ~ strength
obs$origin_strength <- as.integer(obs$origin_strength)
obs[which(obs$origin_strength >= 20),]$origin_strength <- 20

obs$target_strength <- as.integer(obs$target_strength)
obs[which(obs$target_strength >= 20),]$target_strength <- 20

barplot(table(obs$origin_strength), xlab="Force de contact", ylab="Fréquence", main=paste(sp_origin, "interactions"))

obs_conserv <- data.frame(inter = sapply(levels(as.factor(obs$origin_strength)), function(x) (nrow(obs[which(obs$origin_strength == x & !is.na(obs$target_dist)),])/ nrow(obs[which(obs$origin_strength == x),]))*100))
obs_conserv$int_start <- sapply(levels(as.factor(obs$origin_strength)), function(x)  (prop.test(x = nrow(obs[which(obs$origin_strength == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$origin_strength == x),]), p=0.5)$conf.int[1])*100)
obs_conserv$int_end <- sapply(levels(as.factor(obs$origin_strength)), function(x)  (prop.test(x = nrow(obs[which(obs$origin_strength == x & !is.na(obs$target_dist)),]), n=nrow(obs[which(obs$origin_strength == x),]), p=0.5)$conf.int[2])*100)

plot(obs_conserv$inter, type="b", col="red", cex=0.5, main=paste(sp_origin,"interactions conserved in", sp_target),
     ylab="Proportion interaction conserv (%)", xlab="Force du contact (médiane)", xaxt='n', ylim=c(10,35))
for (row in 1:nrow(obs_conserv)){
  segments(x0=row,y0=obs_conserv[row,]$int_start,x1=row,y1=obs_conserv[row,]$int_end, col='red', lwd=0.3)}

axis(1, at=seq(1,16,5), labels=F)
text(seq(1,16,5), par("usr")[3]-0.5, labels = c("5", "10", "15", "20+"), pos = 1, xpd = TRUE)

# Nb tissu target ~ origin
boxplot(obs[which(!is.na(obs$target_nb_tissu)),]$origin_nb_tissu~obs[which(!is.na(obs$target_nb_tissu)),]$target_nb_tissu, 
        outline=F, notch=T, xlab=paste("Nb cell type", sp_target), ylab=paste("Nb cell type", sp_origin),
        main=paste(sp_origin, "interactions conserved in", sp_target))

boxplot(obs[which(!is.na(obs$target_strength)),]$target_strength~obs[which(!is.na(obs$target_strength)),]$origin_strength, 
        outline=F, notch=T, ylab=paste("Force de contact", sp_target), xlab=paste("Force de contact", sp_origin),
        main=paste(sp_origin, "interactions conserved in", sp_target))

########## Chromosome 10 mouse
obs$chr <- as.factor(sub(":.*", "", obs$origin_interaction))

plot(obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]$bait_score~obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]$origin_dist, cex=0.1, 
     ylab="Bait score", xlab="Origin_dist", main="chr10", xlim=c(0,10000000))
plot(obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]$origin_bait_end~obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]$origin_PIR_end, cex=0.1, 
     ylab='Bait position', xlab="PIR position", main="chr10")
abline(a=0,b=1)

data_chr10 <- obs[which(obs$chr == "chr10" & obs$origin_dist > 8000000),]
data_chr10$origin_bait <- as.factor(sub("-.*", "", data_chr10$origin_interaction))
data_chr10$origin_PIR <- as.factor(sub(".*-", "", data_chr10$origin_interaction))
data_chr10 <- data_chr10[,3:4]

write.table(data_chr10, "mouse_chromosome10.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)
