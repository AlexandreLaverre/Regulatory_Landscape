###########################################################################################

options(stringsAsFactors = FALSE)

source("../main_figures/parameters.R") ## pathFinalData are defined based on the user name

###########################################################################################

load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.common.cells.expdiv.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.common.cells.regland.conservation.RData", sep=""))

###########################################################################################

## Defining cell types & samples

cells <- c("embryonic stem cells", "pre-adipocytes", "B lymphocytes")
syncells <- c("ESC", "adipo", "Bcell")

common.cells <- lapply(cells, function(x) lapply(c("human", "mouse"), function(y) {z=sampleinfo[[y]]; z[which(z$Broad.cell.type.or.tissue==x),"Sample.ID"]}))

names(common.cells)=syncells

for(c in syncells){
  names(common.cells[[c]])=c("human", "mouse")
}

##################################################################################################

enh = "ENCODE"

##################################################################################################
############################  Parallel trends among cell types  ##################################

for(sp in c("human", "mouse")){
  
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.",sp,".RData", sep=""))

  enh_align=list_align_enh[[enh]][["enh_align_obs"]]
  
  target_sp=setdiff(c("human", "mouse"), sp)
  
  rownames(ortho)=ortho[,sp]
  
  gene_dnds <- matrix(ncol=3, nrow=3, dimnames=list(syncells, c("Mean","Conf_low", "Conf_high")))
  enh_evol <- matrix(ncol=3, nrow=3, dimnames=list(syncells, c("Mean","Conf_low", "Conf_high")))
  seq_conserv <- matrix(ncol=3, nrow=3, dimnames=list(syncells, c("Mean","Conf_low", "Conf_high")))
  synteny_conserv <- matrix(ncol=3, nrow=3, dimnames=list(syncells, c("Mean","Conf_low", "Conf_high")))
  contact_conserv <- matrix(ncol=3, nrow=3, dimnames=list(syncells, c("Mean","Conf_low", "Conf_high")))
  conserv_expression <- matrix(ncol=3, nrow=3, dimnames=list(syncells, c("Mean","Conf_low", "Conf_high")))
  correl_expression <- matrix(ncol=2, nrow=3, dimnames=list(syncells, c("Pearson","Spearman")))
  correl_complexity <- matrix(ncol=2, nrow=3, dimnames=list(syncells, c("Pearson","Spearman")))
  
  enhancers_contact = enhancer.statistics[[sp]][[enh]][["original"]]
  
  for (cell in syncells){
    # Enhancers contacted by gene in cell types
    enhancers_contact_in_cell <- enhancers_contact[which(apply(enhancers_contact[common.cells[[cell]][[sp]]], 1, function(x) any(x>0))),]$enh
    enh_alignment <- enh_align[which(enh_align$ID %in% enhancers_contact_in_cell),][,target_sp]
    enh_evol[cell,] <- c(mean(enh_alignment, na.rm=T), t.test(enh_alignment)[["conf.int"]][1], t.test(enh_alignment)[["conf.int"]][2])
    
    # dN/dS of the top 25% expressed genes 
    thirdquantile <- quantile(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1), p=0.75)
    expdiv_top <- expdiv_cells[which(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1) > thirdquantile),]
    dNdS <- 1-ortho[which(ortho[,"human"] %in% expdiv_top$IDHuman),]$dNdS
    gene_dnds[cell,] <- c(mean(dNdS, na.rm=T), t.test(dNdS)[["conf.int"]][1],  t.test(dNdS)[["conf.int"]][2])
    
    ###### Correlation of expression ###### 
    correl_expression[cell,] = c(cor(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1), log2(expdiv_cells[[paste(cell, target_sp, "MeanRPKM", sep="_")]]+1), method="pearson"),
                                 cor(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1), log2(expdiv_cells[[paste(cell, target_sp, "MeanRPKM", sep="_")]]+1), method="spearman"))
    
    ###### Conservation of expression ###### 
    t = t.test(expdiv_cells[[paste0(cell, "_ResidualExpressionConservation")]])
    conserv_expression[cell,] =  c(mean(expdiv_cells[[paste0(cell, "_ResidualExpressionConservation")]], na.rm=T), t[["conf.int"]][1], t[["conf.int"]][2])
    
    ######  Conservation of regulatory landscape ###### 
   
    regland_cell =  regland.conservation[[sp]][[enh]][[cell]]
  
    seq_conserv[cell,] <- c(mean(regland_cell$mean.aln.score, na.rm=T), t.test(regland_cell$mean.aln.score)[["conf.int"]][1],  t.test(regland_cell$mean.aln.score)[["conf.int"]][2])
    contact_conserv[cell,] <- c(mean(regland_cell$fr.contact.cons, na.rm=T), t.test(regland_cell$fr.contact.cons)[["conf.int"]][1],  t.test(regland_cell$fr.contact.cons)[["conf.int"]][2])
    synteny_conserv[cell,] <- c(mean(regland_cell$fr.synteny.cons, na.rm=T), t.test(regland_cell$fr.synteny.cons)[["conf.int"]][1],  t.test(regland_cell$fr.synteny.cons)[["conf.int"]][2])
    
    
    ###### Calcul complexity zscore ###### 
    regland_cell$zscore <- (regland_cell$nb.contacts-mean(regland_cell$nb.contacts, na.rm=T)) / sd(regland_cell$nb.contacts, na.rm=T)
    common=intersect(rownames(ortho), rownames(regland_cell))
    regland_cell=regland_cell[common,]
    ortho_cell=ortho[common,]
    
    # Get complexity zscore of orthologous gene in same cell
   
    regland_target = regland.conservation[[target_sp]][[enh]][[cell]] 
    
    regland_target$zscore <- (regland_target$nb.contacts-mean(regland_target$nb.contacts, na.rm=T)) / sd(regland_target$nb.contacts, na.rm=T)

    regland_target=regland_target[ortho_cell[,target_sp],]
    
    ## Correlation between regulatory landscapes complexity
    correl_complexity[cell,] = c(cor(regland_cell$zscore, regland_target$zscore, method="pearson", use="complete.obs"),
                                 cor(regland_cell$zscore, regland_target$zscore, method="spearman", use="complete.obs"))
  }
  
  ### Output
  save(gene_dnds, enh_evol, seq_conserv, contact_conserv, synteny_conserv, correl_expression, correl_complexity, conserv_expression,
       file = paste(pathFigures, "/RData/", sp, ".cells.types.parallel.trends.RData", sep=""))
  
}

##################################################################################################################################
