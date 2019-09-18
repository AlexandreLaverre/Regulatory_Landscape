##### DATA #####
setwd("/home/laverre/Documents/Regulatory_Landscape/result/alignments/human2mouse/")
obs_enh <- read.table("PIR_cons_all_overlap_PECAN.txt", header=T)

sp_origin = "human"
species <- c("mouse", "dog", "cow", "elephant", "opossum", "chicken")

data <- c()
test <- c()
test_enh <- c()

for (sp_target in species){
  setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")
  obs <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv.txt2", sep=""), header=T)
  simul <- read.table(paste(sp_origin,"2",sp_target,"_conservation_syntenie_with_notconserv_simul.txt2", sep=""), header=T)
  
  #### Filtres
  # Length
  obs <- obs[which(obs$bait_length > 250 & obs$bait_length < 20000 & obs$PIR_length >250 & obs$PIR_length < 20000),]
  simul <- simul[which(simul$bait_length > 250 & simul$bait_length < 20000 & simul$PIR_length >250 & simul$PIR_length < 20000),]
  # Duplication
  obs <- obs[which(obs$bait_dupli == 0 & obs$PIR_dupli == 0),]
  simul <- simul[which(simul$bait_dupli == 0 & simul$PIR_dupli == 0),]
  
  obs <- obs[which(!is.na(obs$bait_lift)),]
  obs <- obs[which(!is.na(obs$PIR_lift)),]
  simul <- simul[which(!is.na(simul$bait_lift)),]
  simul <- simul[which(!is.na(simul$PIR_lift)),]
  
  obs$PIR <- as.factor(sub(".*-", "", obs$origin_interaction))
  obs$CAGE <- obs$PIR %in%  obs_enh[which(obs_enh$CAGE_count > 0),]$PIR
  
  synt_obs = nrow(obs[which(!is.na(obs$target_dist) & obs$target_dist < 10000000),])
  synt_obs_CAGE = nrow(obs[which(!is.na(obs$target_dist) & obs$target_dist < 10000000 & obs$CAGE == TRUE),])
  synt_simul = nrow(simul[which(!is.na(simul$target_dist) & simul$target_dist < 10000000),])
  
  mat <- matrix(c(synt_obs, synt_simul, nrow(obs)-synt_obs, nrow(simul)-synt_simul),2)
  mat_enh <- matrix(c(synt_obs, synt_obs_CAGE, nrow(obs)-synt_obs, nrow(obs[which(obs$CAGE == TRUE),])-synt_obs_CAGE),2)
  test <- append(test, prop.test(mat)$p.value)
  test_enh <- append(test_enh, prop.test(mat_enh)$p.value)
  
  data <- append(data, c( (synt_simul/nrow(simul))*100, (synt_obs/nrow(obs))*100, (synt_obs_CAGE/nrow(obs[which(obs$CAGE == TRUE),]))*100, NA))

  message(paste(sp_target,"done !", sep=" "))
}

null = rep(c(0,0,0,NA),6) 

plot(null, col = NA, ylim=c(60,103), xaxt='n', ylab='Conserved synteny (%)', xlab="", cex.lab=1.2)
points(data, col=rep(c("blue", "red", "forestgreen", "black"),6), pch=3, cex=2)
axis(1, at=seq(2,24,4), labels=F)
text(seq(2,24,4), par("usr")[3]-0.5, labels = c("Mouse", "Dog", "Cow", "Elephant", "Opossum", "Chicken"), pos = 1, xpd = TRUE, cex=1.2)
legend("topright", fill=c("blue","red", "forestgreen"), legend = c("Simulated", "Oberved", "With enhancers"), bty='n')

segments(1, 97, 2, 97) 
text("***", x= 1.5, y=98, cex=1.3)
segments(2, 99, 3, 99) 
text("*", x= 2.5, y=100, cex=1.3)

segments(5, 99, 6, 99) 
text("***",x=5.5, y=100, cex=1.3)
segments(6, 101, 7, 101) 
text("**", x=6.5, y=102, cex=1.3)

segments(9, 98, 10, 98) 
text("***",x= 9.5, y=99, cex=1.3)
segments(10, 100, 11, 100) 
text("*", x= 10.5, y=101, cex=1.3)

segments(13, 96, 14, 96) 
text("***",x= 13.5, y=97, cex=1.3)
segments(14, 98, 15, 98) 
text("NS", x=14.5, y=99, cex=1)

segments(17, 84, 18, 84)
text("***",x= 17.5, y=85, cex=1.3)
segments(18, 86, 19, 86) 
text("NS", x= 18.5, y=87, cex=1)

segments(21, 73, 22, 73)
text("***",x= 21.5, y=74, cex=1.4)
segments(22, 75, 23, 75) 
text("***", x= 22.5, y=76, cex=1.4)
