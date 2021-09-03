##########################################################################

pathMappability="/beegfs/data/necsulea/RegulatoryLandscapes/results/mappability/"
pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset1/"

##########################################################################

genomes=c("hg38", "mm10")
names(genomes)=c("human", "mouse")

##########################################################################

lengths=c(50, 51, 101, 69, 75)

for(sp in c("human", "mouse")){
  genome=genomes[sp]

  frag=read.table(paste(pathFinalData, sp, "/frag_coords_",genome, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

  for(rl in lengths){
    system(paste("cut -f 1-7 ",pathMappability, "readlength", rl, "/", sp, "/overlap_restriction_fragments_mapped_regions_margin1000_step10.txt > ",pathMappability, "readlength", rl, "/", sp, "/overlap_restriction_fragments_mapped_regions.txt", sep=""))
    
    map=read.table(paste(pathMappability, "readlength", rl, "/", sp, "/overlap_restriction_fragments_mapped_regions.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    rownames(map)=map$IDFragment
    
    frag[,paste("mappable_length_",rl,sep="")]=map[frag$ID, "TotalMappableLength"]
    frag[,paste("max_mappable_stretch_",rl,sep="")]=map[frag$ID, "MaxMappableStretch"]
  }

  write.table(frag, file=paste(pathFinalData, sp, "/mappability_statistics.txt", sep=""))
}

##########################################################################

