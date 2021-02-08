###########################################################################
library(stringr)
source("parameters.R")

pathContacts=paste(pathFinalData, "SupplementaryDataset4/", sep="")
pathAlignements=paste(pathFinalData, "SupplementaryDataset7/", sep="")

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

###########################################################################

gene.enhancers.align=list()
cells <- c("B lymphocytes", "embryonic stem cells", "pre-adipocytes") 
names(cells) <- c("Bcell", "ESC", "adipo")
enhancers="ENCODE"

for(sp in c("human")){
  if (sp == "human"){target_sp = "mouse"}else{target_sp="human"}
  
  info=sampleinfo[[sp]]
  rownames(info)=info$Sample.ID
  
  samples=info$Sample.ID 
  celltypes=info$Broad.cell.type.or.tissue
  names(celltypes)=samples
  
  for(enh in enhancers){ #enhancer.datasets[[sp]]
    alignment = read.table(paste(pathAlignements, sp, "/sequence_conservation/enhancers/", enh, "/AlignmentStatistics_Excluding_Exons_", sp, "2", target_sp, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    rownames(alignment) = str_sub(alignment[[paste0("ID.", sp)]], end=-3)
    alignment$score = alignment$FilteredUngappedLength/alignment$FilteredAlignmentLength
    alignment[is.nan(alignment$score),"score"] <- 0
    real=read.table(paste(pathContacts, sp, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

    for (cell in cells){

      selected_cell = names(which(celltypes == cell))
      contacts.in.cell <- which(real[,selected_cell]!="NaN")
      real.subdata <- real[contacts.in.cell, c("gene", "enhancer")]
      
      real.subdata$alignmentscore <- alignment[real.subdata$enhancer, "score"]
      real.subdata[is.na(real.subdata$alignmentscore) & !is.nan(real.subdata$alignmentscore),]$alignmentscore <- 0
      
      # sim=read.table(paste(pathContacts, sp, "/", enh, "/gene_enhancer_contacts_simulated_interactions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
      # contacts.in.cell <- which(sim[,selected_cell]!="NaN", arr.ind=T)
      # sim.subdata <- sim[contacts.in.cell[,1], c("gene", "enhancer")]
      # 
      # sim.subdata$alignmentscore <- alignment[sim.subdata$enhancer, "score"]
      # sim.subdata[is.na(sim.subdata$alignmentscore) & !is.nan(sim.subdata$alignmentscore),]$alignmentscore <- 0
      # 
      cellID = names(which(cells == cell))
      gene.enhancers.align[[sp]][[enh]][[cellID]]=list("real"=real.subdata) #, "simulated"=sim)
    }

  }
}

###########################################################################

save(list=c("gene.enhancers.align"), file=paste(pathFigures, "RData/data.gene.enhancer.align.RData", sep=""))

###########################################################################


