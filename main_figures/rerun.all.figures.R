
## all scripts that generate figures

fig.scripts=c("Figure1.R", "Figure2.R", "Figure3.R", "Figure4.R", "Figure5.R") ## "Figure6.R")

############################################################################

for(file in fig.scripts){
  ## cleanup
  objects=ls()
  objects=setdiff(objects, c("fig.scripts", "file"))
  rm(list=objects)

  ## run the script
  print(file)
  source(file)
}

############################################################################


