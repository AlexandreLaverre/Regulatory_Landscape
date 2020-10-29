options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name

path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))

sp="human"
if (sp == "human"){target_sp = "mouse"}else{target_sp = "human"}

enh="ENCODE"
cells <- c("Bcell", "ESC", "adipo")

#############################################################################################################
################################### Expression divergence ###################################################
expdiv_cells=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CellTypes_MeanTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
colnames(expdiv_cells) = c("IDHuman", "IDMouse", "adipo", "Bcell", "ESC")
expdiv_cells <- expdiv_cells[complete.cases(expdiv_cells), ] # to remove if we find a solution to NA
if (sp == "human"){rownames(expdiv_cells) = expdiv_cells$IDHuman}else{rownames(expdiv_cells) = expdiv_cells$IDMouse}

#select genes: protein-coding genes
annot=gene.annot[[sp]]
coding_genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
genes=intersect(rownames(expdiv_cells), coding_genes)
expdiv_cells = expdiv_cells[genes,]

exp_human = read.table(paste(path_exp, "human/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
exp_mouse = read.table(paste(path_exp, "mouse/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")

if(sp == "human"){rownames(expdiv_cells)=expdiv_cells$IDHuman}else{rownames(expdiv_cells)=expdiv_cells$IDMouse}

for (cell in cells){
  samples = grep(cell, colnames(exp_human), value=T)
  expdiv_cells[[paste0("human_", cell, "_Mean")]] <- apply(exp_human[expdiv_cells[,"IDHuman"], samples], 1, function(x) mean(x))
  
  if (cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intitulate B1, B2...
  samples = grep(cell_name, colnames(exp_mouse), value=T)
  expdiv_cells[[paste0("mouse_", cell, "_Mean")]] <- apply(exp_mouse[expdiv_cells[,"IDMouse"], samples], 1, function(x) mean(x))
  
  expdiv_cells[[paste0(cell, "_Mean")]] <- (expdiv_cells[[paste0("human_", cell, "_Mean")]] + expdiv_cells[[paste0("mouse_", cell, "_Mean")]])/2
  
  lm1 = lm(1-expdiv_cells[,cell]~log2(expdiv_cells[[paste0(cell, "_Mean")]]+1))
  expdiv_cells[[paste0(cell, "_ResidualExpressionConservation")]] = lm1$residuals
}

#########################################################################################################################
############################################# Regulatory Divergence ######################################################

data_cell <- list()
for (cell in cells){
  regland_cell = read.table(paste(path_evol, sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, "_0.4.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  
  regland_cell <- regland_cell[which(regland_cell$nb_total >= 5 & regland_cell$nb_total <= 100),] 
  common=intersect(rownames(expdiv_cells), rownames(regland_cell))
  expdiv_cell=expdiv_cells[common,]
  regland_cell=regland_cell[common,]
  
  # Made decile of nb enhancers
  regland_cell$class_nb_contact <- cut2(regland_cell$nb_total, g=5)
  
  regland_cell$ratio_cons_seq = regland_cell$nb_seq_conserv/regland_cell$nb_total
  regland_cell$ratio_cons_synt = ifelse(regland_cell$nb_seq_conserv > 5 & regland_cell$nb_seq_conserv < 100, regland_cell$nb_synt2M_conserv/regland_cell$nb_seq_conserv, NA)
  regland_cell$ratio_cons_int = ifelse(regland_cell$nb_seq_conserv > 5 & regland_cell$nb_seq_conserv < 100, regland_cell$nb_contact_conserv/regland_cell$nb_seq_conserv, NA)
  
  regland_cell$class_align_score=cut2(regland_cell$med_align_score, g=5, include.lowest=T)
  regland_cell$class_cons_seq=cut(regland_cell$ratio_cons_seq, breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  regland_cell$class_cons_synt=cut(regland_cell$ratio_cons_synt,  breaks=c(0, 0.75, 0.99, 1), include.lowest=T)
  regland_cell$class_cons_int=cut(regland_cell$ratio_cons_int,  breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  
  regland_cell$Conservation <- 1-expdiv_cell[[cell]] 
  regland_cell$ResidualConservation <- expdiv_cell[[paste0(cell, "_ResidualExpressionConservation")]]
  
  data_cell[[cell]] <- regland_cell
}

#########################################################################################################################
############################  Parallel trends among cell types  #######################################################
ortho <- read.table(paste(path_evol, sp, "/gene_orthology/human2mouse_orthologue_dNdS.txt", sep="/"), h=T, sep="\t")
rownames(ortho) <- ortho$GenestableID
ortho$dNdS <- ortho$dN/ortho$dS
ortho <- ortho[which(!is.na(ortho$dNdS) & ortho$dNdS < 50),]

gene_dnds <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
enh_evol <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
contact_conserv <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
correl_expression <- matrix(ncol=2, nrow=3, dimnames=list(cells, c("Pearson","Spearman")))
correl_complexity <- matrix(ncol=2, nrow=3, dimnames=list(cells, c("Pearson","Spearman")))

enhancers_alignment = read.table(paste(path_evol, sp, "sequence_conservation/enhancers", enh, "Alignments_stats_all_species_nonexonic_ungapped.txt", sep="/"), h=T)
enhancers_contact = read.table(paste(pathFinalData, "SupplementaryDataset4", sp, enh,  "statistics_contacted_enhancers_original.txt", sep="/"), h=T, stringsAsFactors = F)
enhancers_contact$enh <-  do.call(paste,c(enhancers_contact[c("chr","start","end")],sep=":"))

for (cell in cells){
  if (cell == "Bcell"){cell_name="Bcell"}
  if (cell == "adipo"){cell_name="pre_adipo"}
  if (cell == "ESC"){cell_name="hESC"}
  
  # Enhancers contacted by gene in cell types
  enhancers_contact_in_cell <- enhancers_contact[which(enhancers_contact[[paste0(cell_name)]] > 0),]
  enh_alignment <- enhancers_alignment[which(enhancers_alignment$enh %in% enhancers_contact_in_cell$enh),]$mouse
  enh_evol[cell,] <- c(mean(enh_alignment), t.test(enh_alignment)[["conf.int"]][1], t.test(enh_alignment)[["conf.int"]][2])
  
  # dN/dS of more expressed genes 
  thirdquantile <- summary(log2(expdiv_cells[[paste(sp, cell, "Mean", sep="_")]]))[5]
  expdiv_top <- expdiv_cells[which(log2(expdiv_cells[[paste(sp, cell, "Mean", sep="_")]]) > thirdquantile ),]
  dNdS <- ortho[which(ortho$GenestableID %in% expdiv_top$IDHuman),]$dNdS
  gene_dnds[cell,] <- c(mean(dNdS), t.test(dNdS)[["conf.int"]][1],  t.test(dNdS)[["conf.int"]][2])
  
  ################## Regulatory landscape conservation ############
  regland_cell = read.table(paste(path_evol, sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  common=intersect(rownames(expdiv_cells), rownames(regland_cell))
  expdiv_cell=expdiv_cells[common,]
  regland_cell=regland_cell[common,]
  
  ###### Correlation of expression ###### 
  correl_expression[cell,] = c(cor(log2(expdiv_cell[[paste(sp, cell, "Mean", sep="_")]]+1), log2(expdiv_cell[[paste(target_sp, cell, "Mean", sep="_")]]+1), method="pearson"),
                               cor(log2(expdiv_cell[[paste(sp, cell, "Mean", sep="_")]]+1), log2(expdiv_cell[[paste(target_sp, cell, "Mean", sep="_")]]+1), method="spearman"))
  
  ###### Calcul complexity zscore ###### 
  regland_cell$zscore <- (regland_cell$nb_total-mean(regland_cell$nb_total)) / sd(regland_cell$nb_total)
  
  # Get complexity zscore of orthologous gene in same cell
  regland_target = read.table(paste(path_evol, target_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  
  regland_target=regland_target[expdiv_cell$IDMouse,]
  regland_target[is.na(regland_target)] <- 0
  regland_cell$zscore_target <- (regland_target$nb_total-mean(regland_target$nb_total)) / sd(regland_target$nb_total)
  
  ## Correlation between regulatory landscapes complexity
  correl_complexity[cell,] = c(cor(regland_cell$zscore, regland_cell$zscore_target, method="pearson"),
                               cor(regland_cell$zscore, regland_cell$zscore_target, method="spearman"))
  
  ######  Ratio conserved contact ###### 
  contact <- regland_cell$nb_contact_conserv/regland_cell$nb_total
  contact_conserv[cell,] <- c(mean(contact), t.test(contact)[["conf.int"]][1], t.test(contact)[["conf.int"]][2])
  
}

#########################################################################################################################
# Output

save(data_cell, gene_dnds, enh_evol, contact_conserv, correl_expression, correl_complexity,
     file = paste(pathFigures, "/RData/Fig6_", sp, "_common_cells.Rdata", sep=""))

