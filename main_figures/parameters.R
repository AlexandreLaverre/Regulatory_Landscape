#######################################################################################
## define paths

user=as.character(Sys.getenv()["USER"])

if(user=="laverre"){
  pathFinalData="/home/laverre/Manuscript/"
}

if(user=="necsulea"){
  pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
}

## current directory - we run scripts from the scripts/main_figures folder

dirs=unlist(strsplit(getwd(), split="\\/"))
pathScripts=paste(dirs[1:(length(dirs)-1)], collapse="/")

pathFigures=paste(pathFinalData, "Figures/", sep="")

#######################################################################################

enhancer.datasets=list()
enhancer.datasets[["human"]]=c("ENCODE", "FANTOM5", "FOCS_GRO_seq", "RoadmapEpigenomics")
enhancer.datasets[["mouse"]]=c("ENCODE", "FANTOM5")

#######################################################################################

minDistance=25e3
maxDistance=2.5e6

#######################################################################################
