
pathEvolution="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset7/"

#####################################################################################

enhancer.datasets=list()
enhancer.datasets[["human"]]=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")
enhancer.datasets[["mouse"]]=c("ENCODE", "FANTOM5")

exclude=c("length_enh", "repeat_part", "GC_rate", "BLAT_match", "align_score")

#####################################################################################

allspecies=c("human", "macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")

#####################################################################################

for(ref in c("human","mouse")){
  for(tg in setdiff(allspecies, ref)){
    for(enh in enhancer.datasets[[ref]]){
      for(type in c("original", "simulated")){
        print(paste(ref, tg, enh, type))
        
        data=read.table(paste(pathEvolution, ref, "/synteny_conservation/",enh,"/",ref,"2",tg,"_",type,"_synteny.txt",sep=""), h=T, stringsAsFactors=F)
        data=data[,which(!(colnames(data)%in%exclude))]
        
        system(paste("mv ", pathEvolution, ref, "/synteny_conservation/",enh,"/",ref,"2",tg,"_",type,"_synteny.txt ",pathEvolution, ref, "/synteny_conservation/",enh,"/backup_",ref,"2",tg,"_",data,"_synteny.txt ", sep=""))
        
        write.table(data, file=paste(pathEvolution, ref, "/synteny_conservation/",enh,"/",ref,"2",tg,"_",type,"_synteny.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
      }
    }
  }
}

#####################################################################################
