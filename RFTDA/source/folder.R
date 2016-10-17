#--------------------------------------------------------------------------#
# Folders
#--------------------------------------------------------------------------#
ooftdafolder = "~/BackUpBox/PhD/OOFTDA/"
ftincppdafolder = "~/BackUpBox/PhD/FourierTaylorInCpp/"

folder = list(
  plot   = paste0(ooftdafolder,"plot/"),
  fprint = paste0(ooftdafolder,"fprint/")
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

