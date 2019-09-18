setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

sp_origin = 'human'
sp_target = 'mouse'

obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_by_bait.txt2", sep=""), header=T)
simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_by_bait_simul.txt2", sep=""), header=T)

obs$type <- "Observed"
simul$type <- "Simulated"
data <- rbind(obs, simul)
data$type <- as.factor(data$type)

### Hist of proportion frequency ###
par(mfrow=c(2,2))
data$seq_class <-cut(data$PIR_lift/data$PIR, breaks=seq(from=0, to=1, by=0.05), include.lowest = T)
barplot(table(data$type, data$seq_class), beside=T, col=c("red", "blue"), main="Conservation sequence")
legend("topright", title="Data", legend=c("Observed", "Simulated"), fill=c("red", "blue"), bty='n')

data$intra_class <-cut(data$PIR_synt/data$PIR_lift, breaks=seq(from=0, to=1, by=0.05), include.lowest = T)
barplot(table(data$type, data$intra_class), beside=T, col=c("red", "blue"), main="Conservation synteny")
legend("topright", title="Data", legend=c("Observed", "Simulated"), fill=c("red", "blue"), bty='n')

data$inter_class <-cut(data$PIR_int/data$PIR, breaks=seq(from=0, to=1, by=0.05), include.lowest = T)
barplot(table(data$type, data$inter_class), beside=T, col=c("red", "blue"), main="Conservation interaction")
legend("topright", title="Data", legend=c("Observed", "Simulated"), fill=c("red", "blue"), bty='n')

#### Calcul class ####
data$seq_class <- (data$PIR_lift/data$PIR)*100
data$intra_class <- (data$PIR_synt/data$PIR_lift)*100
data$inter_class <- (data$PIR_int/data$PIR_lift)*100

### Boxplot conserv by bait ###
par(mfrow=c(1,3))
boxplot(data$seq_class~data$type, notch=T, outline=F, col=c("red", "blue"), boxwex=0.7,
        xlab="", ylab="Conservation (%)", main="Sequence")
boxplot(data$intra_class~data$type, notch=T, outline=F, col=c("red", "blue"), boxwex=0.7,
        xlab="", ylab="Conservation (%)", main="Synteny")
boxplot(data$inter_class~data$type, notch=T, outline=F, col=c("red", "blue"), boxwex=0.7,
        xlab="", ylab="Conservation (%)", main="Interaction")

### Conserv ~ TSS number + LM ### 
par(mfrow=c(2,2))
plot(data[which(data$type == "Observed"),]$seq_class~data[which(data$type == "Observed"),]$TSS, cex=0.1, col=rgb(1,0,0, alpha = 0.2),
     ylab="Conservation (%)", xlab="TSS number", main="Sequence")
abline(lm(data[which(data$type == "Observed"),]$seq_class~data[which(data$type == "Observed"),]$TSS), col="red")
points(data[which(data$type == "Simulated"),]$seq_class~data[which(data$type == "Simulated"),]$TSS, cex=0.1, col=rgb(0,0,1, alpha = 0.2))
abline(lm(data[which(data$type == "Simulated"),]$seq_class~data[which(data$type == "Simulated"),]$TSS), col="blue")

plot(data[which(data$type == "Observed"),]$intra_class~data[which(data$type == "Observed"),]$TSS, cex=0.1, col=rgb(1,0,0, alpha = 0.2),
     ylab="Conservation (%)", xlab="TSS number", main="Synteny")
abline(lm(data[which(data$type == "Observed"),]$intra_class~data[which(data$type == "Observed"),]$TSS), col="red")
points(data[which(data$type == "Simulated"),]$intra_class~data[which(data$type == "Simulated"),]$TSS, cex=0.1, col=rgb(0,0,1, alpha = 0.2))
abline(lm(data[which(data$type == "Simulated"),]$intra_class~data[which(data$type == "Simulated"),]$TSS), col="blue")

plot(data[which(data$type == "Observed"),]$inter_class~data[which(data$type == "Observed"),]$TSS, cex=0.1, col=rgb(1,0,0, alpha = 0.2),
     ylab="Conservation (%)", xlab="TSS number", main="Interaction")
abline(lm(data[which(data$type == "Observed"),]$inter_class~data[which(data$type == "Observed"),]$TSS), col="red")
points(data[which(data$type == "Simulated"),]$inter_class~data[which(data$type == "Simulated"),]$TSS, cex=0.1, col=rgb(0,0,1, alpha = 0.2))
abline(lm(data[which(data$type == "Simulated"),]$inter_class~data[which(data$type == "Simulated"),]$TSS), col="blue")

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("center", legend=c("Observed", "Simulated"), fill=c("red", "blue"), bty='n', cex=1.5,
       title="Proportion of conserved contact\nby bait according to TSS")

### Dvpt process TSS
data$dvpt <- "other"
data[which(data$TSS_dvpt > 0 & data$TSS_immun == 0),]$dvpt <- "dvpt"
data[which(data$TSS_dvpt == 0 & data$TSS_immun > 0),]$dvpt <- "immun"
data$dvpt <- as.factor(data$dvpt)

par(mfrow=c(1,3))
# Sequence
boxplot(data$seq_class~data$type, xlim = c(0.5, length(levels(data$type))+0.6), outline=F, boxfill=rgb(1, 1, 1, alpha=1),
        border=rgb(1, 1, 1, alpha=1), ylim=c(0,100), main="Sequence", ylab='Conservation by bait (%)', xlab="") 

boxplot(data[which(data$dvpt == "dvpt"),]$seq_class~data[which(data$dvpt == "dvpt"),]$type, xaxt = "n", yaxt='n', 
        border=c("red","blue"), add = TRUE, boxfill="green", boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type)) - 0.25) 

boxplot(data[which(data$dvpt == "other"),]$seq_class~data[which(data$dvpt == "other"),]$type, xaxt = "n", yaxt='n', 
        border=c("red", "blue"), add = TRUE, boxfill='grey', boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type))) 

boxplot(data[which(data$dvpt == "immun"),]$seq_class~data[which(data$dvpt == "immun"),]$type, xaxt = "n", yaxt='n',
        border=c("red", "blue"), add = TRUE, boxfill='orange', boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type)) + 0.25) 

# Synteny
boxplot(data$intra_class~data$type, xlim = c(0.5, length(levels(data$type))+0.6), outline=F, boxfill=rgb(1, 1, 1, alpha=1),
        border=rgb(1, 1, 1, alpha=1), ylim=c(75,100), main="Synteny", xlab="", ylab="") 

boxplot(data[which(data$dvpt == "dvpt"),]$intra_class~data[which(data$dvpt == "dvpt"),]$type, xaxt = "n", yaxt='n',
        border=c("red", "blue"), add = TRUE, boxfill="green", boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type)) - 0.25) 

boxplot(data[which(data$dvpt == "other"),]$intra_class~data[which(data$dvpt == "other"),]$type, xaxt = "n", yaxt='n',
        border=c("red", "blue"), add = TRUE, boxfill='grey', boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type))) 

boxplot(data[which(data$dvpt == "immun"),]$intra_class~data[which(data$dvpt == "immun"),]$type, xaxt = "n", yaxt='n',
        border=c("red", "blue"), add = TRUE, boxfill='orange', boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type)) + 0.25) 

# Interaction
boxplot(data$inter_class~data$type, xlim = c(0.5, length(levels(data$type))+0.6), outline=F, boxfill=rgb(1, 1, 1, alpha=1),
        border=rgb(1, 1, 1, alpha=1), ylim=c(0,100), main="Interaction", xlab="", ylab="") 

boxplot(data[which(data$dvpt == "dvpt"),]$inter_class~data[which(data$dvpt == "dvpt"),]$type, xaxt = "n", yaxt='n',
        border=c("red", "blue"), add = TRUE, boxfill="green", boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type)) - 0.25) 

boxplot(data[which(data$dvpt == "other"),]$inter_class~data[which(data$dvpt == "other"),]$type, xaxt = "n", yaxt='n',
        border=c("red", "blue"), add = TRUE, boxfill='grey', boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type))) 

boxplot(data[which(data$dvpt == "immun"),]$inter_class~data[which(data$dvpt == "immun"),]$type, xaxt = "n", yaxt='n',
        border=c("red", "blue"), add = TRUE, boxfill='orange', boxwex=0.2, outline=F, notch=T, at = 1:length(levels(data$type)) + 0.25) 

legend("topright", legend=c("Dvpt", "Other", "Immun"), fill=c("green", "grey", "orange"), bty='n')
