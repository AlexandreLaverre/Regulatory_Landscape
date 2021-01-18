###########################################################################################

options(stringsAsFactors = FALSE)

source("parameters.R") ## pathFinalData are defined based on the user name

path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

###########################################################################################

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.common.cells.expdiv.Rdata", sep=""))

ortho <- read.table(paste(path_evol, "human/gene_orthology/human2mouse_orthologue_dNdS.txt", sep="/"), h=T, sep="\t")

## Defining cells types samples
cells <- c("Bcell", "ESC", "adipo")
ESC.common = list("human"= c("hESC"), "mouse"=c("ESC", "ESC_18", "ESC_wild"))
adipo.common = list("human"= c("pre_adipo"), "mouse"=c("preadip_D0", "preadip_D2", "preadip_4H"))
Bcell.common = list("human"= c("TB", "NB"), "mouse"=c("preB_aged", "preB_young"))
common.cell <- list("ESC"=ESC.common, "adipo"=adipo.common, "Bcell"=Bcell.common)

enh = "ENCODE"


##################################################################################################
############################  Parallel trends among cell types  ##################################
for(sp in c("human", "mouse")){
  target_sp=setdiff(c("human", "mouse"), sp)
  
  if (sp == "human"){rownames(ortho) = ortho$GenestableID}else{rownames(ortho) = ortho$GenestableIDMouse}
  
  ortho$dNdS <- ortho$dN/ortho$dS
  ortho <- ortho[which(!is.na(ortho$dNdS) & ortho$dNdS < 50),]
  
  gene_dnds <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
  enh_evol <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
  seq_conserv <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
  synteny_conserv <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
  contact_conserv <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
  conserv_expression <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
  correl_expression <- matrix(ncol=2, nrow=3, dimnames=list(cells, c("Pearson","Spearman")))
  correl_complexity <- matrix(ncol=2, nrow=3, dimnames=list(cells, c("Pearson","Spearman")))
  
  enhancers_alignment = read.table(paste(path_evol, sp, "sequence_conservation/enhancers", enh, "Alignments_stats_all_species_nonexonic_ungapped.txt", sep="/"), h=T)
  enhancers_contact = enhancer.statistics[[sp]][[enh]][["original"]]
  
  for (cell in cells){
    # Enhancers contacted by gene in cell types
    enhancers_contact_in_cell <- enhancers_contact[which(apply(enhancers_contact[common.cell[[cell]][[sp]]], 1, function(x) any(x>0))),]$enh
    enh_alignment <- enhancers_alignment[which(enhancers_alignment$enh %in% enhancers_contact_in_cell),][[target_sp]]
    enh_evol[cell,] <- c(mean(enh_alignment), t.test(enh_alignment)[["conf.int"]][1], t.test(enh_alignment)[["conf.int"]][2])
    
    # dN/dS of more expressed genes 
    thirdquantile <- summary(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]))[5]
    expdiv_top <- expdiv_cells[which(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]) > thirdquantile ),]
    dNdS <- 1-ortho[which(ortho$GenestableID %in% expdiv_top$IDHuman),]$dNdS
    gene_dnds[cell,] <- c(mean(dNdS), t.test(dNdS)[["conf.int"]][1],  t.test(dNdS)[["conf.int"]][2])
    
    ###### Correlation of expression ###### 
    correl_expression[cell,] = c(cor(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1), log2(expdiv_cells[[paste(cell, target_sp, "MeanRPKM", sep="_")]]+1), method="pearson"),
                                 cor(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1), log2(expdiv_cells[[paste(cell, target_sp, "MeanRPKM", sep="_")]]+1), method="spearman"))
    
    ###### Conservation of expression ###### 
    b = boxplot(expdiv_cells[[paste0(cell, "_ResidualExpressionConservation")]], plot=F)
    conserv_expression[cell,] =  c(b$stats[3], b$conf[1], b$conf[2])
    
    ######  Conservation of regulatory landscape ###### 
    load(paste(pathFigures, "RData/data.", sp, ".common.cells.regland.conservation.RData", sep=""))
    regland_cell =  genes.conservation.cells[["ENCODE"]][[cell]][["obs"]][["all"]]
  
    seq_conserv[cell,] <- c(mean(regland_cell$ratio_cons_seq, na.rm=T), t.test(regland_cell$ratio_cons_seq)[["conf.int"]][1],  t.test(regland_cell$ratio_cons_seq)[["conf.int"]][2])
    contact_conserv[cell,] <- c(mean(regland_cell$ratio_cons_int, na.rm=T), t.test(regland_cell$ratio_cons_int)[["conf.int"]][1],  t.test(regland_cell$ratio_cons_int)[["conf.int"]][2])
    synteny_conserv[cell,] <- c(mean(regland_cell$ratio_cons_synt, na.rm=T), t.test(regland_cell$ratio_cons_synt)[["conf.int"]][1],  t.test(regland_cell$ratio_cons_synt)[["conf.int"]][2])
    
    
    ###### Calcul complexity zscore ###### 
    regland_cell$zscore <- (regland_cell$nb_total-mean(regland_cell$nb_total, na.rm=T)) / sd(regland_cell$nb_total, na.rm=T)
    common=intersect(rownames(ortho), rownames(regland_cell))
    regland_cell=regland_cell[common,]
    ortho_cell=ortho[common,]
    
    # Get complexity zscore of orthologous gene in same cell
    load(paste(pathFigures, "RData/data.", target_sp, ".common.cells.regland.conservation.RData", sep=""))
    regland_target = genes.conservation.cells[["ENCODE"]][[cell]][["obs"]][["all"]]
    
    regland_target$zscore <- (regland_target$nb_total-mean(regland_target$nb_total, na.rm=T)) / sd(regland_target$nb_total, na.rm=T)
    
    if (target_sp == "human"){targetID="GenestableID"}else{targetID="GenestableIDMouse"}
    regland_target=regland_target[ortho_cell[[targetID]],]
    
    ## Correlation between regulatory landscapes complexity
    correl_complexity[cell,] = c(cor(regland_cell$zscore, regland_target$zscore, method="pearson", use="complete.obs"),
                                 cor(regland_cell$zscore, regland_target$zscore, method="spearman", use="complete.obs"))
  }
  
  ### Output
  save(gene_dnds, enh_evol, contact_conserv, synteny_conserv, correl_expression, correl_complexity, conserv_expression,align_score,
       file = paste(pathFigures, "/RData/", sp, ".cells.types.parallel.trends.Rdata", sep=""))
  
}
