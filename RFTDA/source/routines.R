#--------------------------------------------------------------------------#
# Misc routines on CRTBP/QBCP/BCP
#--------------------------------------------------------------------------#

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
# From NC coordinates to C coordinates (Earth-Moon coordinates, 
# but still centered at Li)
#--------------------------------------------------------------------------#
NCtoC <- function(df, gamma, c1)
{
  if(missing(c1))
  {
    df$xEM <- -gamma * (df$x);#+6.8859163415402);
  }else
  {
    df$xEM <- -gamma * (df$x-c1);
  }
  df$yEM <- -gamma * (df$y);
  df$zEM <- +gamma * (df$z);
  return(df)
}


# From EM (or C) to physical units
CtoPH <- function(df, L)
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