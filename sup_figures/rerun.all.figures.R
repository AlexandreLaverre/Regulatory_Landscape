
## all scripts that generate figures

all.files=system("ls", intern=T)
fig.scripts=grep("Figure", all.files, value=T)

fig.scripts.slow=c("SupplementaryFigure3_SupplementaryMaterialFigure22.R","SupplementaryMaterialFigure17.R", "SupplementaryMaterialFigure18.R", "SupplementaryMaterialFigure19.R","SupplementaryMaterialFigure20.R", "SupplementaryMaterialFigure32-33.R","SupplementaryMaterialFigure8.R")

fig.scripts.fast=setdiff(fig.scripts, fig.scripts.slow)

############################################################################

for(file in fig.scripts.fast){
  ## cleanup
  objects=ls()
  objects=setdiff(objects, c("fig.scripts.fast", "fig.scripts.slow", "file"))
  rm(list=objects)

  ## run the script
  print(file)
  source(file)
}

############################################################################

print("Do you want to redo figures with BCa computations?")

consent=scan(what=character(), nmax=1)

if(consent=="yes"){
  
  for(file in fig.scripts.slow){
    ## cleanup
    objects=ls()
    objects=setdiff(objects, c("fig.scripts.slow", "file"))
    rm(list=objects)
    
    ## run the script
    print(file)
    source(file)
  }
}

############################################################################


