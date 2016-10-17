# R script to handle a precision map of the QBFBP around EML1,2
#---------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------
# Init
#------------------------------------------------
source("source/init.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li     = "L2"
MODEL  = "QBCP"
FWRK   = "EM"
Type   = "rand" #rand or NF_rand
Energy = 0
OFS_ORDER = 30
vorders = c(5, 10, 15, 20, 25, 30);
currentfolder = paste0(printfolder(MODEL, FWRK, Li))
pS = 2

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
muR = muR(FWRK);
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK);
if(FWRK == "EM")
{
  primaryR    =  1737.10      #m2
}else{
  primaryR    =  6378.10      #m2 
}

#------------------------------------------------
#Additionnal parameters
#------------------------------------------------
isLegendOn = 1;
legendOnly = 0;

#------------------------------------------------
# Building the data.frame of results
#------------------------------------------------
imap = data.frame()
for (i in vorders)  #loop on the orders
{
  if(Energy == 0)
  {
    fileprefix = paste0(currentfolder, "Serv/eIm_", Type, "_ofs_", toString(OFS_ORDER),"_order_",toString(i), "_planar");
  }else
  {
    fileprefix = paste0(currentfolder, "Serv/eIm_", Type, "_ofs_", toString(OFS_ORDER),"_order_",toString(i), "_energy_", toString(Energy));
  }
  filename = paste0(fileprefix, ".bin")
  
  
  # Load csv source
  #------------------------------------------------
  if (file.exists(filename))
  {
    names = c("label", "x", "y", "z", "px", "py", "pz",
              "xEM", "yEM", "zEM", "pxEM", "pyEM", "pzEM",
              "s1", "s2", "s3", "s4", "t", "dHz", "eIm");
    imapc = dffbinary(filename, 20, names);
    imapc$order = i;
    
  }else
  {
    imapc = data.frame()
  }
  
  #rbind in ttm_all
  imap = rbind(imap, imapc);
}

#------------------------------------------------
# Postprocessing
#------------------------------------------------
# Using EIm (EM units) instead of eIm (in NC units) 
imap$EIm = gamma*imap$eIm 
# To centered units
imap = NCtoC(imap, gamma)
# To physical units
imap = CtoPH(imap, L)
# To physical units
imap = EMtoPH(imap, L)
# Compute -log10(precision)
imap$log10eOm = log10(imap$EIm)
# Compute -log10(precision)
imap$flog10eOm = floor(log10(imap$EIm))

#------------------------------------------------
# Get rid of bad precision
#------------------------------------------------
imap = imap[which(imap$log10eOm < -2),]

#------------------------------------------------
# Select only positive some given value in the results
#------------------------------------------------
eps = 1e-7;
isCloseToEPS = abs(imap$EIm - eps) < 0.2*eps;
imapr = imap[which(isCloseToEPS),]
#Ordered version
imapr = imapr[order(imapr$order),]


#---------------------------------------------------------------------------------------------------------------------------
# Plot 
#---------------------------------------------------------------------------------------------------------------------------
#---------------------
# Only a given precision
#---------------------
# Version without background
pplc = plotdf_point(imapr, "xEM", "yEM", "$X$", "$Y$", pointSize = pS, colorCol = "order", colorLabel = "Order", isColorFac = TRUE)

#Version with background (at constant energy)
# pplc = ggplot() + geom_point(data = imap, aes(x = xEM, y = yEM), color = "grey", size = pS, alpha = 0.01)
# pplc = pplc + labs(x = "$X$ $[$-$]$", y = "$Y$ $[$-$]$")
# pplc = pplc + geom_point(data = imapr, aes(x = xEM, y = yEM, color = factor(order)))+ scale_color_discrete(name="Order")
# pplc

#---------------------
#Add m2
#---------------------
primaryPos  =  -(1-muR)  #Add m2
# moon = data.frame(xEM = primaryPos, yEM = 0)
# pplc = pplc + geom_point(data = moon, aes(xEM, yEM), size = 4, colour = "black", fill = "black", pch = 21)
#Add m2 with true apparen radius
moon = circleFun(center = c(primaryPos, 0), diameter = primaryR/L, npoints = 100, start = 0, end = 2, filled = TRUE)
pplc = pplc + geom_polygon(data=moon, aes(x,y), color="black", fill="black")

#---------------------
#Ratio and limits
#---------------------
pplc = pplc + coord_fixed(ratio=1)
if(Li == "L2")
{
  pplc = pplc + scale_y_continuous(limits = c(-0.115,0.115))
  pplc = pplc +legend_inside_theme
}else{
  pplc = pplc + scale_y_continuous(limits = c(-0.115,0.115))
  pplc = pplc +legend_left_theme
}

#---------------------
# Adding the resonnances
#---------------------
if(Li == "L2")
{
  #------------------------------------------------
  # Orbit 24 (rough param)
  #------------------------------------------------
  xc  = c(-1.129, -1.13, -1.136, -1.14, -1.15, -1.16, -1.17, -1.175, -1.177)
  xc2 = c(-1.177, -1.175, -1.17, -1.16, -1.15, -1.14, -1.136, -1.13, -1.129)
  yc  = c(0, 0.031, 0.060, 0.068, 0.069, 0.060, 0.040, 0.020, 0)
  yc2 = c(0, 0.020, 0.040, 0.060, 0.069, 0.068, 0.060, 0.031, 0)
  df24 = data.frame(x = c(xc,  xc2),
                    y = c(yc, -yc2))
  pplc = pplc + geom_path(data=df24, aes(x=x, y=y), size=1, linetype= "dashed")
}else{
  eml1_res = read.table("Single_Orbit/EML1_resonant_orbit.txt", header = F)
  names(eml1_res) = c("x", "y")
  pplc = pplc + geom_path(data=eml1_res, aes(x=x, y=y), size=1, linetype= "dashed")
}

#---------------------
# Other possibilities
#---------------------
ppls = plotdf_smooth(imapr, "xEM", "yEM", "X", "Y", colorCol = "order", colorLabel = "order", isColorFac = TRUE)
ppls

ppll = plotdf_line(imapr, "xEM", "yEM", "X", "Y", colorCol = "order", colorLabel = "order", isColorFac = TRUE)
ppll



#---------------------
# Only a given precision, in physical units
#---------------------
pplcph = plotdf_point(imapr, "xCPH", "yCPH", "X [km]", "Y [km]", pointSize = pS, colorCol = "order", colorLabel = "order", isColorFac = TRUE)
pplcph

pplph = plotdf_point(imapr, "xPH", "yPH", "X [km]", "Y [km]", pointSize = pS, colorCol = "order", colorLabel = "order", isColorFac = TRUE)
#Add m2
primaryPos  =  -(1-muR) * L #Add m2
pplph = addMoon(pplph, x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2
)
pplph

#---------------------
# Only a given order
#---------------------
ppl = plotdf_point(imap[which(imap$order == 5),], "xEM", "yEM", "X", "Y", pointSize = pS, colorCol = "flog10eOm", colorLabel = "log10eOm", isColorFac = TRUE)
ppl

#---------------------
# With respect to the energy
#---------------------
#Select only some given value in the results
eps = 0.015;
isCloseToEPS = abs(imap$dHz - eps) < 1e-4;
imape = imap[which(isCloseToEPS),]
ppe = plotdf_point(imap[which(imap$order == 20),], "xEM", "yEM", "X", "Y", pointSize = pS, colorCol = "flog10eOm", colorLabel = "log10eOm", isColorFac = TRUE)
ppe

#---------------------------------------------------------------------------------------------------------------------------
# Adding the CR3BP zero-velocity curves
#---------------------------------------------------------------------------------------------------------------------------
if(Energy > 0)
{
  #---------------------
  #Build vectors
  #---------------------
  if(Li == "L2")
  {
    x0   = seq(-2,-1,0.001)
  }else{
    x0   = seq(-1,-0.8,0.001)
  }
  y0   = seq(-0.15,0.15,0.001)
  
  z0   = expand.grid(x0, y0)
  potn = apply(z0, 1, cr3bpot_EM)
  
  
  #Build data frame
  df = data.frame(x1 = z0[1], y1 = z0[2], pot = potn)
  names(df) = c('x', 'y', 'pot')
  
  #---------------------
  #Energy associated with the Libration points
  #---------------------
  l1 = c(-0.836915145386502, 0, 0)
  l2 = c(-1.155682150113637, 0, 0)
  l1pot = cr3bpot(l1[1], l1[2], muR)
  l2pot = cr3bpot(l2[1], l2[2], muR)
  
  #---------------------
  #Select a given energy
  #---------------------
  if(Li == "L2")
  {
    eps = l2pot+Energy;
  }else{
    eps = l1pot+Energy;
  }
  isCloseToEPS = abs(df$pot - eps) < 1e-4;
  dfr = df[which(isCloseToEPS),]
  
  #---------------------
  #Plot
  #---------------------
  pplcp = pplc + geom_point(data=dfr, aes(x=x, y=y), size=3)
  pplcp + coord_fixed(ratio=1)
}


#---------------------
# s1-s3 maps
#---------------------
ppls1s3 = plotdf_point(imapr, "s1", "s3", "$s_1$", "$s_3$", pointSize = pS, colorCol = "order", colorLabel = "Order", isColorFac = TRUE)
if(Li == "L2")
{
  ppls1s3 = ppls1s3 + scale_x_continuous(breaks = seq(-50,50,4)) + scale_y_continuous(breaks = seq(-50,50,4))
}else{
  ppls1s3 = ppls1s3 + scale_x_continuous(breaks = seq(-2,2,0.5)) + scale_y_continuous(breaks = seq(-2,2,0.5))
}
ppls1s3

#Save in tikz
#--------------------------------------------------------------------------------------------------------------------------
ggplot2tikz(pplc, xSize, ySize+2, file = paste0(currentfolder, "EIm_", Type, "_ofs_", toString(OFS_ORDER),"_order_",toString(i), "_planar_resonant.tex"))
