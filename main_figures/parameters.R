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

enh.syn=c("ENCODE", "FANTOM5", "FOCS GRO-seq", "Roadmap")
names(enh.syn)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")

#######################################################################################
## colors for the plots

col.enhancers=c("red", "navy", "forestgreen", "orange") ## colors for the datasets
names(col.enhancers)=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")

dataset.colors=c("navy", "gray50") ##c("forestgreen", "firebrick1")
names(dataset.colors)=c("Original", "Simulated")

col.Shh="forestgreen"

#######################################################################################

minDistance=25e3
maxDistance=2.5e6

#######################################################################################
