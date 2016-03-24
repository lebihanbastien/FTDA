#--------------------------------------------------------------------------#
# Folders
#--------------------------------------------------------------------------#
homefolder = "~/BackUpBox/PhD/OOFTDA/"

folder = list(
  plot   = paste0(homefolder,"plot/"),
  fprint = paste0(homefolder,"fprint/")
)

#--------------------------------------------------------------------------#
# Routines to build folders name
#--------------------------------------------------------------------------#
printfolder<- function(Model,   #Model
                       Fwrk,
                       Li      #Libration point
                      )    
{
  return(paste0(folder[["fprint"]], Model, "/", Fwrk, "/", Li,"/"))
}


plotfolder<- function(Model,   #Model
                       Fwrk,
                       Li      #Libration point
)    
{
  return(paste0(folder[["plot"]], Model, "/", Fwrk, "/", Li,"/"))
}

