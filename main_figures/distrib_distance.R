options(stringsAsFactors = FALSE)

library(data.table)

ref_sp = "human"
target_sp = "mouse" 
path <- "/home/laverre/Data/Regulatory_landscape/result"

path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_contact <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")

#pdf(paste(path_evol, "/", ref_sp, "2", target_sp, "_enhancer_contacted_genes_by_distance.pdf", sep=""))
######################### All interactions ####################
par(mfrow=c(2,1))
par(mai = c(1, 1, 0.5, 1))
### All samples ###
contact_obs <- fread(paste(path, "Supplementary_dataset1_original_interactions/human/all_interactions.txt", sep="/"), header=T)
contact_obs <- contact_obs[which(contact_obs$type == "unbaited" & contact_obs$chr_bait == contact_obs$chr),]
contact_obs <- contact_obs[which(contact_obs$dist < 10000000 & contact_obs$dist > 25000),]
contact_obs$dist <- contact_obs$dist/1000
contact_obs$nb_sample <- apply(contact_obs[,9:34], 1, function(x) sum(!is.na(x)))

contact_simul <- fread(paste(path, "Supplementary_dataset2_simulated_interactions/human/simulated_all_interactions.txt", sep="/"), header=T)
contact_simul <- contact_simul[which(contact_simul$type == "unbaited" & contact_simul$chr_bait == contact_simul$chr),]
contact_simul <- contact_simul[which(contact_simul$distance < 10000000 & contact_simul$distance > 25000),]
contact_simul$distance <- contact_simul$distance/1000
contact_simul$nb_sample <- apply(contact_simul[,9:34], 1, function(x) sum(!is.na(x)))

x <- density(contact_obs$dist)
plot(x, xlab="Distance (kb)", main="Bait-other interactions in all samples (human)",  col="darkgreen", xlim=c(0,1500))
lines(density(contact_simul$distance, bw=x$bw), col="firebrick3")
legend("topright", legend=c("Original", "Simulated"),fill=c("darkgreen", "firebrick3"), bty='n')

### Bcell sample ###
Bcell_obs <- contact_obs[which(!is.na(contact_obs$Bcell)),]
Bcell_simul <- contact_simul[which(!is.na(contact_simul$Bcell)),]

x <- density(Bcell_obs$dist)
plot(x, xlab="Distance (pb)", main="Bait-other interactions in Bcell (human)",  col="darkgreen", xlim=c(0,1500))
lines(density(Bcell_simul$dist, bw=x$bw), col="firebrick3")
legend("topright", legend=c("Original", "Simulated"),fill=c("darkgreen", "firebrick3"), bty='n')

######################### Gene - enhancer interactions ####################
### All samples ###
contact_obs <- fread(paste(path, "Supplementary_dataset4_genes_enhancers_contacts/human/gene_ENCODE_enhancers_original_interactions.txt", sep="/"), header=T)
contact_obs <- contact_obs[which(contact_obs$dist < 10000000 & contact_obs$dist > 25000),]
contact_obs$dist <- contact_obs$dist/1000

contact_simul <- fread(paste(path, "Supplementary_dataset4_genes_enhancers_contacts/human/gene_ENCODE_enhancers_simulated_interactions.txt", sep="/"), header=T)
contact_simul <- contact_simul[which(contact_simul$dist < 10000000 & contact_simul$dist > 25000),]
contact_simul$dist <- contact_simul$dist/1000

x <- density(contact_obs$dist)
plot(x, xlab="Distance (kb)", main="Genes-enh interactions in all samples (human)",  col="darkgreen", xlim=c(0,1500))
lines(density(contact_simul$dist, bw=x$bw), col="firebrick3")
legend("topright", legend=c("Original", "Simulated"),fill=c("darkgreen", "firebrick3"), bty='n')

### Bcell ###
Bcell_obs <- contact_obs[which(!is.na(contact_obs$Bcell)),]
Bcell_simul <- contact_simul[which(!is.na(contact_simul$Bcell)),]

x <- density(contact_obs$dist)
plot(x, xlab="Distance (kb)", main="Genes-enh interactions in all samples (human)",  col="darkgreen", xlim=c(0,1500))
lines(density(contact_simul$dist, bw=x$bw), col="firebrick3")
legend("topright", legend=c("Original", "Simulated"),fill=c("darkgreen", "firebrick3"), bty='n')

