setwd("/home/laverre/Documents/Regulatory_Landscape/result/conservation/")

ref_sp = "mouse"
target_sp = "human"

obs_obs <- read.table(paste(ref_sp,"_merged_to_",target_sp,"_summary_by_gene.txt2", sep=""), header= T, sep="\t", row.names=1)
obs_simul <- read.table(paste(ref_sp,"_merged_to_",target_sp,"_simul_summary_by_gene.txt2", sep=""), header= T, sep="\t", row.names=1)
simul_obs <- read.table(paste(ref_sp,"_simul_merged_to_",target_sp,"_summary_by_gene.txt2", sep=""), header= T, sep="\t", row.names=1)
simul_simul <- read.table(paste(ref_sp,"_simul_merged_to_",target_sp,"_simul_summary_by_gene.txt2", sep=""), header= T, sep="\t", row.names=1)

data_list <- list(obs_obs, simul_obs, obs_simul, simul_simul)
data_names <- c("obs2obs", "simul2obs", "obs2simul", "simul2simul")
data_list_obs_simul <- list(obs_obs, simul_obs)
data_names_obs_simul <- c("obs", "simul")

par(mfrow=c(2,2))
par(mai = c(0.5, 0.8, 0.5, 0.8))
conserv_seq <- as.data.frame(sapply(data_list_obs_simul, function(x) x$nb_seq_conserv/x$nb_contact))
colnames(conserv_seq) <- data_names_obs_simul
boxplot(conserv_seq, main=paste(ref_sp, " Sequence Conservation", sep=""), ylab="Proportion", outline=F, notch=T, boxwex=0.5)


conserv_seq_50 <- as.data.frame(sapply(data_list_obs_simul, function(x) x$seq_conserv_50/x$nb_contact))
colnames(conserv_seq_50) <- data_names_obs_simul
boxplot(conserv_seq_50, main=paste(ref_sp, " Sequence Conservation (overlap 50%)", sep=""), ylab="Proportion", outline=F, notch=T, boxwex=0.5)

conserv_synt <- as.data.frame(sapply(data_list_obs_simul, function(x) x$synt_conserv/x$nb_seq_conserv))
colnames(conserv_synt) <- data_names_obs_simul
boxplot(conserv_synt, main=paste(ref_sp, " Synteny Conservation", sep=""), ylab="Proportion", outline=F, notch=T, boxwex=0.5)

conserv_contact <- as.data.frame(sapply(data_list, function(x) x$nb_int_conserv/x$nb_seq_conserv))
colnames(conserv_contact) <- data_names
boxplot(conserv_contact, main=paste(ref_sp, " Contact Conservation", sep=""), ylab="Proportion", outline=F, notch=T, boxwex=0.5)

####################################################################################################################
############################################# Conserved enhancers ##################################################

obs_obs <- read.table(paste(ref_sp,"_merged_to_",target_sp,"_summary_enhancers_by_gene.txt", sep=""), header= T, sep="\t", row.names=1)
obs_simul <- read.table(paste(ref_sp,"_merged_to_",target_sp,"_simul_summary_enhancers_by_gene.txt", sep=""), header= T, sep="\t", row.names=1)
simul_obs <- read.table(paste(ref_sp,"_simul_merged_to_",target_sp,"_summary_enhancers_by_gene.txt", sep=""), header= T, sep="\t", row.names=1)
simul_simul <- read.table(paste(ref_sp,"_simul_merged_to_",target_sp,"_simul_summary_enhancers_by_gene.txt", sep=""), header= T, sep="\t", row.names=1)

data_list <- list(obs_obs, simul_obs, obs_simul, simul_simul)
data_names <- c("obs2obs", "simul2obs", "obs2simul", "simul2simul")
data_list_obs_simul <- list(obs_obs, simul_obs)
data_names_obs_simul <- c("obs", "simul")

enhancers <- c("CAGE") #, "ENCODE", "RoadMap", "GRO_seq")
for (enh in enhancers){
  par(mai = c(0.5, 0.8, 0.5, 0.8))
  layout(matrix(c(1,2,3,4),nrow = 2,byrow = TRUE))
  CAGE <- as.data.frame(sapply(data_list_obs_simul, function(x) x[[paste0("nb_", enh)]]/x$contacted_length))
  colnames(CAGE) <- data_names_obs_simul
  boxplot(CAGE, main=paste(ref_sp, "contacted", enh, sep=" "), ylab="Nb enhancers / Total contacted length", outline=F, notch=T, boxwex=0.5)
  
  CAGE_conserv <- as.data.frame(sapply(data_list_obs_simul, function(x) x[[paste0("nb_", enh, "_conserv")]]/x[[paste0("nb_", enh)]]))
  colnames(CAGE_conserv) <- data_names_obs_simul
  boxplot(CAGE_conserv, main=paste(ref_sp, "conserved",enh, sep=" "), ylab="Proportion", outline=F, notch=T, boxwex=0.5)
  
  CAGE_contact <- as.data.frame(sapply(data_list_obs_simul, function(x) x[[paste0("nb_", enh, "_synt")]]/x[[paste0("nb_", enh, "_conserv")]]))
  colnames(CAGE_contact) <- data_names_obs_simul
  boxplot(CAGE_contact, main=paste(ref_sp, enh, "maintained in synteny", sep=" "), ylab="Proportion", outline=F, notch=T, boxwex=0.5)
  
  CAGE_contact <- as.data.frame(sapply(data_list, function(x) x[[paste0("nb_", enh, "_contact")]]/x[[paste0("nb_", enh, "_conserv")]]))
  colnames(CAGE_contact) <- data_names
  boxplot(CAGE_contact, main=paste(ref_sp, enh, "maintained in contact", sep=" "), ylab="Proportion", outline=F, notch=T, boxwex=0.5)
}

