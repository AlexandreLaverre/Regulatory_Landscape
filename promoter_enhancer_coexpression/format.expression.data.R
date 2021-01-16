########################################################################

path="/beegfs/data/necsulea/RegulatoryLandscapes/"
pathFOCS=paste(path, "data/FOCS/", sep="")

options(stringsAsFactors=F)
options(digits=2) ## to make files lighter

########################################################################

datasets=list()
datasets[["human"]]=c("ENCODE", "FANTOM5", "GRO-seq", "RoadmapEpigenomics")
datasets[["mouse"]]=c("ENCODE", "FANTOM5")

genomes=list()
genomes[["human"]]="hg38"
genomes[["mouse"]]="mm10"

########################################################################

for(sp in c("human", "mouse")){
  for(dataset in datasets[[sp]]){
    prom.exp=read.table(paste(pathFOCS, sp, "/", dataset, "/promoter_activity.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    enh.exp=read.table(paste(pathFOCS, sp, "/", dataset, "/enhancer_activity.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    prom.exp=round(prom.exp, digits=2)
    enh.exp=round(enh.exp, digits=2)

    prom.coords=read.table(paste(pathFOCS, sp, "/", dataset, "/promoter_coordinates_",genomes[[sp]],".bed", sep=""), h=F, stringsAsFactors=F, sep="\t")
    enh.coords=read.table(paste(pathFOCS, sp, "/", dataset, "/enhancer_coordinates_",genomes[[sp]],".bed", sep=""), h=F, stringsAsFactors=F, sep="\t")

    colnames(prom.coords)=c("chr", "start", "end", "id")
    colnames(enh.coords)=c("chr", "start", "end", "id")
    
    rownames(prom.coords)=prom.coords$id
    rownames(enh.coords)=enh.coords$id

    prom.exp=cbind(prom.coords[rownames(prom.exp),], prom.exp)
    enh.exp=cbind(enh.coords[rownames(enh.exp),], enh.exp)

    write.table(prom.exp, file=paste(pathFOCS, sp, "/", dataset, "/promoter_coords_activity.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)
    write.table(enh.exp, file=paste(pathFOCS, sp, "/", dataset, "/enhancer_coords_activity.txt", sep=""),  row.names=T, col.names=T, sep="\t", quote=F)
  }
}

########################################################################
