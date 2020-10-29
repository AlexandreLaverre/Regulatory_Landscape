#######################################################################################

set.seed(19) ## we randomly sample colors, we need this to be reproducible

#######################################################################################

options("stringsAsFactors"=FALSE)

#######################################################################################

## define paths

user=as.character(Sys.getenv()["USER"])

if(user=="laverre"){
  pathFinalData="/home/laverre/Manuscript/"
}

if(user=="necsulea"){
  pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
}


if(user=="ubuntu"){
  pathFinalData="/mnt/RegulatoryLandscapesManuscript/"
}


## current directory - we run scripts from the scripts/main_figures folder

dirs=unlist(strsplit(getwd(), split="\\/"))
pathScripts=paste(dirs[1:(length(dirs)-1)], collapse="/")

pathFigures=paste(pathFinalData, "Figures/", sep="")

#######################################################################################

enhancer.datasets=list()
enhancer.datasets[["human"]]=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")
enhancer.datasets[["mouse"]]=c("ENCODE", "FANTOM5")

label.enhancers=c("ENCODE", "FANTOM5", "FOCS GRO-seq", "Roadmap Epigenomics")
names(label.enhancers)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")

enh.syn=c("ENCODE", "FANTOM5", " FOCS\nGRO-seq", " Roadmap\nEpigenomics")
names(enh.syn)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")

#######################################################################################
## colors for the plots

col.enhancers=c("red", "navy", "forestgreen", "orange") ## colors for the datasets
names(col.enhancers)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")


dataset.colors=c("navy", "gray30") ##c("forestgreen", "firebrick1")
dataset.density=c(10,10)
dataset.angle=c(45,-45)

names(dataset.colors)=c("Original", "Simulated")
names(dataset.density)=c("Original", "Simulated")
names(dataset.angle)=c("Original", "Simulated")

col.Shh=dataset.colors["Original"]

#######################################################################################

col.celltypes=c("gray40", "navy", "blue", "steelblue", "deepskyblue1", "slateblue", "purple", "red", "brown", "indianred", "orange", "darkgoldenrod", "darkorange", "forestgreen", "aquamarine3", "black", "indianred", "steelblue", "purple", "darkred", "darkorange")
names(col.celltypes)=c("B lymphocytes", "cardiomyocytes", "embryonic stem cells", "endothelial precursors", "erythroblasts", "fetal thymus", "hematopoietic progenitors", "keratinocytes", "lymphoblastoid cell line", "macrophages", "megakaryocytes","monocytes", "neuroepithelial cells", "neutrophils", "pre-adipocytes", "T lymphocytes", "embryonic stem cells, Nanog KO",  "epiblast stem cells", "ES-derived  hematopoietic progenitors", "fetal liver", "trophoblast stem cells")

syn.celltypes=c("B lymphocytes", "cardiomyocytes", "embryonic stem cells", "endothelial precursors", "erythroblasts", "fetal thymus", "hematopoietic progenitors", "keratinocytes", "lymphoblastoid cells", "macrophages", "megakaryocytes","monocytes", "neuroepithelial cells", "neutrophils", "pre-adipocytes", "T lymphocytes", "embryonic stem cells, Nanog KO",  "epiblast stem cells", "ES-derived\nhematopoietic progenitors", "fetal liver", "trophoblast stem cells")
names(syn.celltypes)=names(col.celltypes)

#######################################################################################

minDistance=25e3
maxDistance=2e6

#######################################################################################
