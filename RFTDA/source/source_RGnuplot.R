#--------------------------------------------------------------------------#
# Routines for RGnuplot interface
#--------------------------------------------------------------------------#
# Get the position of a string in array
#--------------------------------------------------------------------------#
getPos <- function(style, st, max)
{
  n =  which(style == st)[[1]]
  return(n)
}

# Arrays
#--------------------------------------------------------------------------#
#Point Style
pointStyleNumber = 15;
pointStyle = c("+", "x", "*", "s", "fs", "c", "fc", "ut", "fup", "dt", "fdt", "d", "fd", "h", "fh")
#Line Style
lineStyleNumber  = 7;
lineStyle  = c("-", "--", "---", "..", ".-", ".-.-", "...")
#Color
colorStyleNumber  = 9;
colorStyle = c("dark-violet", "#009e73", "#56b4e9", "#e69f00", "#f0e442", "#0072b2", "#e51e10", "black", "gray50")


# Function for set the main cmd
#--------------------------------------------------------------------------#
CGpsetcmd <- function(plotType = "plot", fileName = "test.data", title = "", lt = "fc", ls = "points", lw = 1, lc = 2, replot = F)
{
  #Get the number for lineStyle from arrays
  if(ls == "points"){
    lti = getPos(pointStyle, lt, pointStyleNumber)
  }else{
    lti = getPos(lineStyle,  lt, lineStyleNumber)
  }
  #Get the color
  lci = colorStyle[(lc-1)%%colorStyleNumber+1]
  #Do we need to replot?
  if(replot)
  {
      plotTypei = "replot"
  }else
  {
     plotTypei = plotType
  }
  
  #Set the cmd
  cmd = paste0(plotTypei, " \"", fileName, "\"", " title ", "\"", 
               title, "\"", " with ", ls, " lt ", "\"", 
               lti, "\"", " lw ", lw, " lc rgb ", "\"", lci, "\"")
  print(cmd)
  return(cmd)
}
