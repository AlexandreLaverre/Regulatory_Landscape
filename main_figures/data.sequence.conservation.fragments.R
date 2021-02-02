###########################################################################

library(data.table)

source("parameters.R") ## paths are defined based on the user name

pathSequenceConservation=paste(pathFinalData, "SupplementaryDataset7/", sep="")

###########################################################################

load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))

###########################################################################

all.species=c("human", "mouse", "macaque", "rat", "rabbit", "dog", "cow", "elephant", "opossum", "chicken")

minAlnLength=10

###########################################################################

for(ref in c("human", "mouse")){

  frag.obs=fragment.statistics[[ref]][["original"]]
  frag.sim=fragment.statistics[[ref]][["simulated"]]

  all.frag=unique(c(frag.obs$ID, frag.sim$ID))

  for(tg in setdiff(all.species, ref)){
    print(paste("fragment sequence conservation", ref, tg))
    
    path=paste(pathSequenceConservation, ref, "/sequence_conservation/restriction_fragments/AlignmentStatistics_Excluding_Exons_", ref,"2", tg,".txt", sep="")

    if(file.exists(path)){

      pcungapped=rep(NA, length(all.frag)) ## we assign NA values to non-lifted elements
      pcidentical=rep(NA, length(all.frag))
      
      names(pcungapped)=all.frag
      names(pcidentical)=all.frag
        
      cons=fread(path, h=T)
      class(cons)="data.frame"

      ## select sequences that are actually aligned (PECAN can return 0 alignments)

      cons=cons[which(cons$TotalAlignmentLength>minAlnLength),]

      ## use only filtered (non-exonic) sequence

      ## ungapped sequence
      this.pcungap=cons$FilteredUngappedLength/cons$FilteredAlignmentLength
      this.pcungap[which(cons$FilteredAlignmentLength < minAlnLength)]=NA
      names(this.pcungap)=cons[,paste("ID", ref, sep=".")]
      
      pcungapped[intersect(names(pcungapped), names(this.pcungap))]=this.pcungap[intersect(names(pcungapped), names(this.pcungap))]

      ## identical sequence
      this.pcid=cons$FilteredIdenticalLength/cons$FilteredUngappedLength
      this.pcid[which(cons$FilteredAlignmentLength < minAlnLength)]=NA
      names(this.pcid)=cons[,paste("ID", ref, sep=".")]
      
      pcidentical[intersect(names(pcidentical), names(this.pcid))]=this.pcid[intersect(names(pcidentical), names(this.pcid))]

      ## save results
      save(list=c("pcidentical", "pcungapped"), file=paste(pathFigures, "RData/data.sequence.conservation.fragments.",ref,"2", tg,".RData", sep=""))
    } else{
      print(paste("cannot find file for ",ref, " and ",tg,sep=""))
    }
  }
}

###########################################################################
