setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/interaction_conservation/")

ref_sp = "mouse"
target_sp = "human"

obs_obs <- read.table(paste(ref_sp,"_merged_to_",target_sp,"_by_gene.txt3", sep=""), header=T)
simul_obs <- read.table(paste(ref_sp,"_simul_merged_to_",target_sp,"_by_gene.txt3", sep=""), header=T)
obs_simul <- read.table(paste(ref_sp,"_merged_to_",target_sp,"_simul_by_gene.txt3", sep=""), header=T)
simul_simul <- read.table(paste(ref_sp,"_simul_merged_to_",target_sp,"_simul_by_gene.txt3", sep=""), header=T)

data_list <- list(obs_obs, simul_obs, obs_simul, simul_simul)
data_names <- c("obs2obs", "simul2obs", "obs2simul", "simul2simul")

## Filters
data_list <- lapply(data_list, function(x) x[which(x$duplication == 0),])

#######################################################################################################################
####################################  Contacts conservation proportion  ###############################################

dev.off()
par(mfrow=c(2,2))
### Contact conserv / all contact
# Unbaited
conserv_all_unbaited <- sapply(data_list, function(x) 
  nrow(x[which(x$status == "interaction_conserved" & x$baited_contact == "unbaited"),])
  /nrow(x[which(x$baited_contact == "unbaited"),]))

barplot(conserv_all_unbaited,names=data_names, ylim=c(0,0.5), ylab="Proportion",
        main=paste(ref_sp, " Conserved contact / All contact (Unbaited)", sep=""))

# Baited
conserv_all_baited <- sapply(data_list, function(x) 
  nrow(x[which(x$status == "interaction_conserved" & x$baited_contact == "baited"),])
  /nrow(x[which(x$baited_contact == "baited"),]))

barplot(conserv_all_baited,names=data_names, ylim=c(0,0.5), ylab="Proportion",
        main=paste(ref_sp," Conserved contact / All contact (Baited", sep=""))

### Contact conserv / seq conserv
# Unbaited
conserv_seq_unbaited <- sapply(data_list, function(x) 
  nrow(x[which(x$status == "interaction_conserved" & x$baited_contact == "unbaited"),])
  /nrow(x[which(!is.na(x$conserved_frag) & x$baited_contact == "unbaited"),]))

barplot(conserv_seq_unbaited,names=data_names, ylim=c(0,0.5), ylab="Proportion",
        main=paste(ref_sp, " Conserved contact / Conserved seq (Unbaited)", sep=""))

# Baited
conserv_seq_baited <- sapply(data_list, function(x) 
  nrow(x[which(x$status == "interaction_conserved" & x$baited_contact == "baited"),])
  /nrow(x[which(!is.na(x$conserved_frag) & x$baited_contact == "baited"),]))

barplot(conserv_seq_baited,names=data_names, ylim=c(0,0.5), ylab="Proportion",
        main=paste(ref_sp, " Conserved contact / Conserved seq (Baited)", sep=""))

#######################################################################################################################
####################################  Contacts conservation through distance  #########################################

data_list <- lapply(data_list, function(x) 
  cbind(x, class = cut(x$dist, breaks=seq(from=25000, to=3000000, by=50000), include.lowest = T)))
class_leg <- c("50Kb", "500Kb", "1Mb", "1.5Mb", "2Mb","2.5Mb")

##### Unbaited #####
conserv_contact_unbaited <- as.data.frame(sapply(data_list, function(x) sapply(levels(x$class), function(y) 
  nrow(x[which(x$status == "interaction_conserved" & x$baited_contact == "unbaited" & x$class == y),])
  /nrow(x[which(!is.na(x$conserved_frag) & x$baited_contact == "unbaited" & x$class == y),]))))

colnames(conserv_contact_unbaited) <- data_names

par(mfrow=c(1,1))
plot(conserv_contact_unbaited$obs2obs[1:50], type="l", col="red", 
     ylab="Conserved contact proportion", xlab="Distance (pb)",
     main=paste(ref_sp, "2", target_sp, " Ortho Genes - Conservation Unbaited Contact", sep=""), 
     xaxt = "n", ylim=c(0,0.6))

lines(conserv_contact_unbaited$obs2simul[1:50], type="l", col="navy")
lines(conserv_contact_unbaited$simul2obs[1:50], type="l", col="forestgreen")
lines(conserv_contact_unbaited$simul2simul[1:50], type="l", col="orange")

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.03, class_leg, xpd = TRUE)
legend("topright", fill=c("red","navy","forestgreen", "orange"), 
       legend = c("obs2obs", "obs2simul", "simul2obs", "simul2simul"), bty='n')

# Difference between obs and simul
plot(conserv_contact_unbaited$obs2obs[1:50]-conserv_contact_unbaited$simul2obs[1:50], type="l", col="black",
     ylab="Conserved contact difference proportion", xlab="Distance (pb)",
     main=paste(ref_sp, "2", target_sp, " Ortho gene (Obs - Simul) contacts conservation", sep=""), xaxt = "n")

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.015, class_leg, xpd = TRUE)

##### Baited #####
conserv_contact_baited <- as.data.frame(sapply(data_list, function(x) sapply(levels(x$class), function(y) 
  nrow(x[which(x$status == "interaction_conserved" & x$baited_contact == "baited" & x$class == y),])
  /nrow(x[which(!is.na(x$conserved_frag) & x$baited_contact == "baited" & x$class == y),]))))

colnames(conserv_contact_baited) <- data_names

par(mfrow=c(1,1))
plot(conserv_contact_baited$obs2obs[1:50], type="l", col="red", 
     ylab="Conserved contact proportion", xlab="Distance (pb)",
     main=paste(ref_sp, "2", target_sp, " Ortho Genes - Conservation Baited Contact", sep=""),
     xaxt = "n", ylim=c(0,0.7))

lines(conserv_contact_baited$obs2simul[1:50], type="l", col="navy")
lines(conserv_contact_baited$simul2obs[1:50], type="l", col="forestgreen")
lines(conserv_contact_baited$simul2simul[1:50], type="l", col="orange")

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.03, class_leg, xpd = TRUE)
legend("topright", fill=c("red","navy","forestgreen", "orange"), 
       legend = c("obs2obs", "obs2simul", "simul2obs", "simul2simul"), bty='n')

# Difference between obs and simul
plot(conserv_contact_baited$obs2obs[1:50]-conserv_contact_baited$obs2simul[1:50], type="l", col="black", 
     ylab="Conserved contact difference (%)", xlab="Distance (pb)",
     main=paste(ref_sp, "2", target_sp, " Ortho gene (Obs - Simul) contacts conservation", sep=""), 
     xaxt = "n", ylim=c(0,0.4))

axis(1, at=seq(1,51,10), labels=F)
text(seq(1,51,10),par("usr")[3]-0.02, class_leg, xpd = TRUE)

  
  
  
  