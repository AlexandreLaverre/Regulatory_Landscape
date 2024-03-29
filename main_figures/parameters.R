#######################################################################################

set.seed(19) ## we randomly sample colors, we need this to be reproducible

#######################################################################################

options("stringsAsFactors"=FALSE)

#######################################################################################

## define paths

user=as.character(Sys.getenv()["USER"])
pathRlibs=NULL

if(user=="laverre"){
  pathFinalData="/home/laverre/Manuscript/"
  pathRepeats="/home/laverre/Regulatory_landscape/result/sequence_composition/"
}


if(user %in% c("necsulea", "alaverre")){
  pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
  pathRlibs="/beegfs/home/necsulea/Rlibs"
  pathRepeats="/beegfs/data/necsulea/RegulatoryLandscapes/results/sequence_composition/"
}

if(user=="ubuntu"){
  pathFinalData="/mnt/mydatalocal/RegulatoryLandscapesManuscript/"
}

## current directory - we run scripts from the scripts/main_figures folder

dirs=unlist(strsplit(getwd(), split="\\/"))
pathScripts=paste(dirs[1:(length(dirs)-1)], collapse="/")

pathFigures=paste(pathFinalData, "tmp_writtable/Figures/", sep="")
pathSuppTables=paste(pathFinalData, "SupplementaryTables/", sep="")

#######################################################################################

enhancer.datasets=list()
enhancer.datasets[["human"]]=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")
enhancer.datasets[["mouse"]]=c("ENCODE", "FANTOM5")

label.enhancers=c("ENCODE", "FANTOM5", "FOCS GRO-seq", "Roadmap Epigenomics")
names(label.enhancers)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")

enh.syn=c("ENCODE", "FANTOM5", "FOCS GRO-seq", "Roadmap Epigenomics")
names(enh.syn)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")

enh.syn.narrow=c("ENCODE", "FANTOM5", "FOCS\nGRO-seq", "Roadmap\nEpigenomics")
names(enh.syn.narrow)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")

#######################################################################################
## colors for the plots

col.enhancers=c("red", "navy", "forestgreen", "orange") ## colors for the datasets
names(col.enhancers)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")


dataset.colors=c("darkorange", "dodgerblue4") ##c("forestgreen", "firebrick1")
names(dataset.colors)=c("Original", "Simulated")
col.Shh=dataset.colors["Original"]

#######################################################################################

col.celltypes=c("gray40", "navy", "blue", "steelblue", "deepskyblue1", "slateblue", "purple", "red", "brown", "indianred", "orange", "darkgoldenrod", "darkorange", "forestgreen", "aquamarine3", "black", "darkblue", "steelblue", "purple", "darkred", "darkorange")
names(col.celltypes)=c("B lymphocytes", "cardiomyocytes", "embryonic stem cells", "endothelial precursors", "erythroblasts", "fetal thymus", "hematopoietic progenitors", "keratinocytes", "lymphoblastoid cell line", "macrophages", "megakaryocytes","monocytes", "neuroepithelial cells", "neutrophils", "pre-adipocytes", "T lymphocytes", "embryonic stem cells, Nanog KO",  "epiblast stem cells", "ES-derived  hematopoietic progenitors", "fetal liver", "trophoblast stem cells")

syn.celltypes=c("B lymphocytes", "cardiomyocytes", "embryonic stem cells", "endothelial precursors", "erythroblasts", "fetal thymus", "hematopoietic progenitors", "keratinocytes", "lymphoblastoid cells", "macrophages", "megakaryocytes","monocytes", "neuroepithelial cells", "neutrophils", "pre-adipocytes", "T lymphocytes", "ESC Nanog KO",  "epiblast stem cells", "ES-derived\nhematopoietic progenitors", "fetal liver", "trophoblast stem cells")
names(syn.celltypes)=names(col.celltypes)

#######################################################################################

minDistance=25e3
maxDistance=2e6

minDistanceSyntenyRef=1e5 ## 100 kb
maxDistanceSyntenyRef=1.5e6 ## 1.5 Mb
maxDistanceSyntenyTarget=2e6 ## 2 Mb

minFragmentSize=150
maxFragmentSize=50000

minSampleTarget = 1
minAlignScore = 0.4 ## for contact conservation

## filters on nb of BLAT hits for fragments and enhancers
minBLAT=-1
maxBLAT=2

## default sequence conservation values, for fragments and enhancers that are not lifted
default.cons = 0

#######################################################################################
# Packages used
# requiredPackages = c('data.table', 'bootBCa', 'ade4', 'plotrix', 'ape', 'vioplot', 'Hmisc', 'imager', 'plyr', 'car', 'MuMIn')
# 
# for(p in requiredPackages){
#   if(!require(p, character.only = TRUE)) install.packages(p)
#   if(!require(p, character.only = TRUE) & p == 'bootBCa') install.packages(p, repos="http://R-Forge.R-project.org")
# }

