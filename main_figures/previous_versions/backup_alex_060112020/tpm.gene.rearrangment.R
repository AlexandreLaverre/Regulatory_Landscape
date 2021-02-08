################################################################################################################################################
options(stringsAsFactors = FALSE)

ref_sp = "human" # to change human or mouse
target_sp = "mouse"

source("parameters.R") 
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", ref_sp, "/", sep="")
path_annot <- paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", sep="")


minDistance=25e3
maxDistance=2e6

enhancers = c("FANTOM5", "ENCODE")
if(ref_sp == "human"){enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")}

enh="FANTOM5"

synt_obs <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", target_sp, "_original_synteny.txt", sep=""), header=T)
synt_simul <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", target_sp, "_simulated_synteny.txt", sep=""), header=T)

# Filters
synt_obs <- synt_obs[which(synt_obs$BLAT_match < 2 & synt_obs$align_score > 0.8 & synt_obs$origin_dist < maxDistance),]
synt_simul <- synt_simul[which(synt_simul$BLAT_match < 2 & synt_simul$align_score > 0.8 & synt_simul$origin_dist < maxDistance),]

## Gene annotations
annot_human <- read.table(paste(pathFinalData,"SupplementaryDataset3/genes/human_genes_Ensembl94.txt", sep=""), header=T, row.names=1)
annot_mouse <- read.table(paste(pathFinalData,"SupplementaryDataset3/genes/mouse_genes_Ensembl94.txt", sep=""), header=T)
ortho <- read.table(paste(pathFinalData, "SupplementaryDataset3/human2mouse_ortholog_one2one_genes_Ensembl94.txt", sep="/"), h=T, sep="\t")
rownames(ortho) <- ortho$GenestableID

### Nb of contacts by genes
contact_by_gene <- as.data.frame(table(synt_obs$origin_gene))
rownames(contact_by_gene) <- contact_by_gene[,1]
contact_by_gene <- contact_by_gene[which(contact_by_gene$Freq > 2),]


#write.table(contact_by_gene, paste(pathFigures, "synteny_all_genes.txt", sep=""), append = FALSE, sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE, quote=FALSE)


# Nb of rearranged contacts
synt_obs_rearranged <- synt_obs[which(synt_obs$target_dist == "trans"),]
contact_rearranged <- as.data.frame(table(synt_obs_rearranged$origin_gene))

rownames(contact_rearranged) <- contact_rearranged[,1]
colnames(contact_rearranged) <- c("Gene", "TotalEnhancerRearranged")

#add nb contact
contact_by_gene <- contact_by_gene[which(contact_by_gene$Var1 %in% contact_rearranged$Gene),]
contact_rearranged$TotalEnhancerConserved <- contact_by_gene[rownames(contact_rearranged),]$Freq
contact_rearranged$Ratio <- contact_rearranged$TotalEnhancerRearranged/contact_rearranged$TotalEnhancerConserved
contact_rearranged <- contact_rearranged[which(contact_rearranged$TotalEnhancerConserved > 5),]

#add human annot
annot_human <- annot_human[which(rownames(annot_human) %in% contact_rearranged$Gene),]
contact_rearranged$Annotation <- annot_human[rownames(contact_rearranged),]$GeneBiotype

#add mouse annot
ortho <- ortho[which(rownames(ortho) %in% contact_rearranged[,1]),]
contact_rearranged$OrthoGene <- ortho[rownames(contact_rearranged),]$GenestableIDMouse

annot_mouse <- annot_mouse[which(annot_mouse$GeneID %in% contact_rearranged$OrthoGene),]
contact_rearranged$OrthoAnnotation <- annot_human[rownames(contact_rearranged),]$GeneBiotype



# write.table(contact_rearranged[order(-contact_rearranged$rearranged),], paste(pathFigures, "synteny_rearranged_genes.txt", sep=""), append = FALSE, sep = "\t", dec = ".",
#             row.names = FALSE, col.names = FALSE, quote=FALSE)
# 

