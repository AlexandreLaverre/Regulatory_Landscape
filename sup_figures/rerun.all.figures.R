
## all scripts that generate figures

all.files=system("ls", intern=T)
fig.scripts=grep("Figure", all.files, value=T)

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


