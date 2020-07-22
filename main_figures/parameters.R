#######################################################################################
## define paths

user=as.character(Sys.getenv()["USER"])

if(user=="alaverre"){

}

if(user=="necsulea"){
  pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
}

## current directory - we run scripts from the scripts/main_figures folder

dirs=unlist(strsplit(getwd(), split="\\/"))
pathScripts=paste(dirs[1:(length(dirs)-1)], collapse="/")

#######################################################################################

