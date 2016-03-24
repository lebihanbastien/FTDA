# R script to handle a precision map of the QBFBP around EML1,2
#---------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------
# R options
#------------------------------------------------
options(digits = 15)

#------------------------------------------------
# Load libraries
#------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)

#------------------------------------------------
# Load Source files
#------------------------------------------------
source("source/source_folder.R")
source("source/source_plot.R")
source("source/source_routines.R")
source("source/multiplot.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li    = "L2"
MODEL = "RTBP"
FWRK  = "EM"
Type  = "s1s3" #selection or global
si = "s1"
sj = "s3"
vorders = c(20, 30)
ofs_order = 0
dHz = 0.5
currentfolder = paste0(printfolder(MODEL, FWRK, Li))

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK);

#------------------------------------------------
#Additionnal parameters
#------------------------------------------------
isLegendOn = 1;
legendOnly = 0;

#------------------------------------------------
# Building the data.frame of results
#------------------------------------------------
ttm_sub = data.frame()
ttm_l = list();

#------------------------------------------------
# Loop on the orders
#------------------------------------------------
for (i in c(1,2))  #loop on the orders
{
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0(currentfolder, "Serv/eIm_", Type, "_ofs_",ofs_order,"_order_", toString(vorders[i]), "_hmax_", toString(dHz))
  filename = paste0(fileprefix, ".txt")
  
  # Load csv source
  #------------------------------------------------
  if (file.exists(filename))
  {
    ttm_l[[i]]  = read.csv(filename, header = T, sep = ",")
  }else
  {
    ttm_l[[i]] = data.frame()
  }
  
  # Postprocessing
  #------------------------------------------------
  #Order
  ttm_l[[i]] = ttm_l[[i]][order(ttm_l[[i]]$s1, ttm_l[[i]]$s2, ttm_l[[i]]$s3),]
  # Compute -log10(precision)
  ttm_l[[i]]$log10eOm = log10(ttm_l[[i]]$eOm)
  # From NC to EM units
  ttm_l[[i]] = NCtoC(ttm_l[[i]], gamma)
  # From EM to physical units
  ttm_l[[i]] = CtoPH(ttm_l[[i]], L)
  # Radii from Li
  ttm_l[[i]]$rNC = sqrt(ttm_l[[i]]$x^2+ttm_l[[i]]$y^2+ttm_l[[i]]$z^2)
  ttm_l[[i]]$rPH = sqrt(ttm_l[[i]]$xPH^2+ttm_l[[i]]$yPH^2+ttm_l[[i]]$zPH^2)
}

#------------------------------------------------
#Substracting
#------------------------------------------------
ttm_sub = data.frame(s1=ttm_l[[1]]$s1, 
                     s2=ttm_l[[1]]$s2, 
                     s3=ttm_l[[1]]$s3, 
                     xPH = ttm_l[[1]]$xPH, 
                     yPH = ttm_l[[1]]$yPH, 
                     zPH = ttm_l[[1]]$zPH, 
                     rPH = ttm_l[[1]]$rPH, 
                     log10eOm=ttm_l[[1]]$log10eOm-ttm_l[[2]]$log10eOm)

#Get rid of big values
ttm_sub = ttm_sub[which(ttm_sub$rPH < 1e5),]

#------------------------------------------------
# Selecting one line
#------------------------------------------------
ttm_sub_line = ttm_sub[which(ttm_sub$s1 == ttm_sub$s3),]

#------------------------------------------------
# Plot (Tile)
#------------------------------------------------

#--------------------
#Plot in s space
#--------------------
ppl = plotdf_tile(ttm_sub , si, sj, isLegendOn, 0)

#--------------------
#Plot in z space
#--------------------
xi = "xPH"
yi = "yPH"
xist = "x [km]"
yist = "y [km]"
primaryFactor = 1.0#-gamma*Le

pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "log10eOm", "log10eOm", 0)+scale_colour_gradient2("log10eIm", space="Lab", midpoint = 0, mid = "white", high = muted("blue"))

if(xi == "xPH")
{
  #Earth direction
  pph    = pph+ annotate("text", x = 57000, y = -25000, label = "to Earth \ ", size = 6)
  pph    = pph+ geom_segment(aes(x = 50000, y = -30000, xend = 60000, yend = -30000), 
                                       colour = "black", 
                                       arrow = arrow(length = unit(0.3, "cm"), 
                                                     type = "closed"))
  #Add the moon
  primaryR    = 1737.10/abs(primaryFactor)  #Moon's radius
  primaryPos  = gamma*L/primaryFactor       #Moon position wrt to Li
  pph    = addMoon(pph, x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2)
  pph    = pph+ annotate("text", x = sign(primaryFactor)*gamma*L/primaryFactor, y = 5000, label = "Moon", size = 6)
}

#Add the line
pph = pph + geom_point(data = ttm_sub_line,  aes(x = xPH, y = yPH), color = "black")
