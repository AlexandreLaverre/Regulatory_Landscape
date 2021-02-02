#####################################################################################

pathEvolution="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset7/"

#####################################################################################

enhancer.datasets=list()
enhancer.datasets[["human"]]=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")
enhancer.datasets[["mouse"]]=c("ENCODE", "FANTOM5")

exclude=c("length_enh", "repeat_part", "GC_rate", "BLAT_match", "align_score")

#####################################################################################

for(ref in c("human","mouse")){
  tg=setdiff(c("human", "mouse"), ref)

  for(enh in enhancer.datasets[[ref]]){
     for(refdata in c("original", "simulated")){
       for(tgdata in c("original", "simulated")){

         print(paste(ref, tg, enh, refdata, tgdata))
         
         data=read.table(paste(pathEvolution, ref, "/contact_conservation/", enh,"/",ref,"_",refdata,"2",tg,"_",tgdata,".txt",sep=""), h=T, stringsAsFactors=F)

         data=data[,which(!(colnames(data)%in%exclude))]

         system(paste("mv ",pathEvolution, ref, "/contact_conservation/", enh,"/",ref,"_",refdata,"2",tg,"_",tgdata,".txt ",pathEvolution, ref, "/contact_conservation/", enh,"/backup_",ref,"_",refdata,"2",tg,"_",tgdata,".txt ",sep=""))
         
         write.table(data, file=paste(pathEvolution, ref, "/contact_conservation/", enh,"/",ref,"_",refdata,"2",tg,"_",tgdata,".txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
       }
     }
   }
}

#####################################################################################
