############################################################################

## all scripts that generate .RData objects for the figures

data.scripts=c("data.Shh.figure.R", "Figure2_A_B_C_create_data.R", "Figure2_D_E_F_create_data.R", "Figure3_create_data.R", "Figure4_create_data.R", "Figure5_create_data.R")

############################################################################

for(file in data.scripts){
  ## cleanup
  objects=ls()
  objects=setdiff(objects, "data.scripts", "file")
  rm(list=objects)

  ## run the script
  source(file)
}

############################################################################


