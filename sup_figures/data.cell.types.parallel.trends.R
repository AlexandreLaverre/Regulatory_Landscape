###########################################################################################

options(stringsAsFactors = FALSE)

source("../main_figures/parameters.R") ## pathFinalData are defined based on the user name

###########################################################################################

load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.common.cells.expdiv.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

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
  target_sp=setdiff(c("human", "mouse"), sp)

  ## enhancer alignment statistics
  
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.",sp,".RData", sep=""))
  enh_align=list_align_enh[[enh]][["enh_align_obs"]]

  ## all contacted enhancers
  enhancers_contact = enhancer.statistics[[sp]][[enh]][["original"]]

  ## prepare results

  enh.alignment=list()
  dNdS=list()
  exp.corr=list()
  
  for (cell in syncells){
    ## Enhancers contacted by gene in this cell type
    enhancers_contact_in_cell <- enhancers_contact[which(apply(enhancers_contact[common.cells[[cell]][[sp]]], 1, function(x) any(x>0))),]$enh
    enh.alignment[[cell]] <- enh_align[which(enh_align$ID %in% enhancers_contact_in_cell),][,target_sp]
   
    ## dN/dS of the top 25% expressed genes 
    thirdquantile <- quantile(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1), p=0.75)
    expdiv_top <- expdiv_cells[which(log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1) > thirdquantile),]
    dNdS[[cell]] <- ortho[which(ortho[,"human"] %in% expdiv_top$IDHuman),]$dNdS
   
    ## Correlation of expression

    meanexp.ref=log2(expdiv_cells[[paste(cell, sp, "MeanRPKM", sep="_")]]+1)
    meanexp.tg=log2(expdiv_cells[[paste(cell, target_sp, "MeanRPKM", sep="_")]]+1)
    
    ## bootstrap replicates

    exp.corr[[cell]]=sapply(1:100, function(x) {s=sample(1:length(meanexp.ref), size=length(meanexp.ref), re=T); return(cor(meanexp.ref[s], meanexp.tg[s], method="spearman", use="complete.obs"))})    
  
  }
  
  ### Output
  save(dNdS, exp.corr, enh.alignment, file = paste(pathFigures, "/RData/", sp, ".cells.types.parallel.trends.RData", sep=""))
  
}

##################################################################################################################################
