#--------------------------------------------------------------------------#
# Misc routines on CRTBP/QBCP/BCP
#--------------------------------------------------------------------------#
ggplot2tikz <- function(plot, width0, height0, file)
{
  #Create a .tex file that will contain the plot
  tikz(file, 
       width = width0, 
       height = height0, 
       standAlone=TRUE, 
       documentDeclaration="\\documentclass{standalone}\n",
       packages = c("\\usepackage[utf8]{inputenc}",
                    "\\usepackage[T1]{fontenc}",
                    "\\usepackage{tikz}", 
                    "\\usepackage{pgf}", 
                    "\\usetikzlibrary{calc}", 
                    "\\usepackage{amssymb}", 
                    "\\usepackage{amsfonts}\n"))
  
  #Plot porb in a trash document
  #plottrash = plot
  #Print it to feed dev
  print(plot)
  #Necessary to close or the tikxDevice .tex file will not be written
  dev.off()
  #return(plottrash)
}

# To stop quietly a code... Used instead of stop() that throw and error.
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()

#--------------------------------------------------------------------------#
# Potential
#--------------------------------------------------------------------------#
cr3bpot <- function(x, y, mu)
{
  r1  = sqrt((x-mu)^2 + y^2);
  r2  = sqrt((x-mu+1)^2 + y^2);
  pot = -0.5*(x^2 + y^2) - (1-mu)/r1 - mu/r2 - 1/2*(1-mu)*mu;
  return (pot)
}

cr3bpot_EM <- function(z)
{
  mu = 0.0121505816234336;
  r1  = sqrt((z[1]-mu)^2 + z[2]^2);
  r2  = sqrt((z[1]-mu+1)^2 + z[2]^2);
  pot = -0.5*(z[1]^2 + z[2]^2) - (1-mu)/r1 - mu/r2 - 1/2*(1-mu)*mu;
  return (pot)
}

#--------------------------------------------------------------------------#
# Return the mass ratio
#--------------------------------------------------------------------------#
muR <- function(FWRK)
{
  if (FWRK == "EM")
  {
    mu = +1.215058162343360e-02;
  }
  
  if (FWRK == "SEM")
  {
    mu = +3.040277404912506e-06;
  }
  
  return(mu)
}

#--------------------------------------------------------------------------#
# Return the value gamma = distance(Li-smallest primary) 
# for L1,2, distance(Li-biggest primary for L3)
#--------------------------------------------------------------------------#
gamma <- function(Li, FWRK)
{
  if (FWRK == "EM")
  {
    if (Li == "L1")
      gamma = 0.150934272990064
    if (Li == "L2")
      gamma = 0.16783273173707
    if (Li == "L3")
      gamma = +9.929120655179929e-01
  }
  
  if (FWRK == "SEM")
  {
    if (Li == "L1")
      gamma = 0.0100108175327472
    if (Li == "L2")
      gamma = 0.0100780785917628
    if (Li == "L3")
      gamma = 0.999998226504847
  }
  
  return(gamma)
}

#--------------------------------------------------------------------------#
# Return the value c1 = f(gamma) for L1,2,3 (see Summary)
#--------------------------------------------------------------------------#
c1 <- function(Li, FWRK)
{
  if (FWRK == "EM")
  {
    if (Li == "L1")
      c1 = -5.5448979798087
    if (Li == "L2")
      c1 = -6.8859163415402
    if (Li == "L3")
      c1 = +1.012237318938304e+00
  }
  
  if (FWRK == "SEM")
  {
    if (Li == "L1")
      c1 = -98.8916378512805
    if (Li == "L2")
      c1 = -100.224961432011
  }
  
  return(c1)
}

#--------------------------------------------------------------------------#
# Return the value L of the distance between the two primaries
#--------------------------------------------------------------------------#
Ldist <- function(FWRK)
{
  if (FWRK == "SEM")
  {
    L = 149.60e6; #Sun-Earth distance in [km]
  }else{
    L = 384400;   #Earth-Moon distance in [km]
  }
  return(L)
}

#--------------------------------------------------------------------------#
# Return the value T of the period associated with the two primaries
#--------------------------------------------------------------------------#
Tcrtbp <- function(FWRK)
{
  if (FWRK == "SEM")
  {
    T = 3.155814950400000e+07;
  }else{
    T = 2.360584684800000e+06;
  }
  return(T)
}

#--------------------------------------------------------------------------#
# Return the value of the SEM period in the FWRK, in normalized units
#--------------------------------------------------------------------------#
SEMperiod <- function(FWRK)
{
  if (FWRK == "SEM")
  {
    T = 5.080085647283307e-01;
  }else{
    T = 6.791193871907917e+00;
  }
  return(T)
}

#--------------------------------------------------------------------------#
# From NC coordinates to C coordinates (Earth-Moon coordinates)
# From NC coordinates to C coordinates (Earth-Moon coordinates, 
# but still centered at Li)
#--------------------------------------------------------------------------#
NCtoC <- function(df, gamma)
{
  df$xC <- -gamma * (df$x);#+6.8859163415402);
  df$yC <- -gamma * (df$y);
  df$zC <- +gamma * (df$z);
  return(df)
}

#--------------------------------------------------------------------------#
# From NC coordinates to SYS coordinates
#--------------------------------------------------------------------------#
NCtoSYS <- function(df, gamma, c1)
{
  df$xEM <- -gamma * (df$x-c1);
  df$yEM <- -gamma * (df$y);
  df$zEM <- +gamma * (df$z);
  return(df)
}


# From EM (or C) to physical units
CtoPH <- function(df, L)
{
  df$xCPH = L * df$xC
  df$yCPH = L * df$yC
  df$zCPH = L * df$zC
  return(df)
}

# From EM to physical units
EMtoPH <- function(df, L)
{
  df$xPH = L * df$xEM
  df$yPH = L * df$yEM
  df$zPH = L * df$zEM
  return(df)
}

# From SYS to physical units
SYStoPH <- function(df, L)
{
  df$xPH = L * df$xEM
  df$yPH = L * df$yEM
  df$zPH = L * df$zEM
  return(df)
}

#--------------------------------------------------------------------------#
# User input
#--------------------------------------------------------------------------#
readInteger <- function(call)
{
  n <- readline(prompt = call)
  n <- as.integer(n)
  if (is.na(n)) {
    n <- readinteger()
  }
  return(n)
}

readString <- function(call)
{
  n <- readline(prompt = call)
  return(n)
}