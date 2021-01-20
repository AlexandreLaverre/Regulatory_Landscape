
## all scripts that generate figures


all.files=system("ls", intern=T)
fig.scripts=c(grep("SupplementaryFigure", all.files, value=T), grep("ExtendedFigure", all.files, value=T))
fig.scripts=grep("SupplementaryFigureX", fig.scripts, value=T, invert=T)

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


