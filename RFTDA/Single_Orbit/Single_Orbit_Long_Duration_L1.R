# Script to plot an orbit along with its associated precisions
# in the parameterization method of the QBFBP/RTBP around L1 of the Earth-Moon system
#
# Used for the computation of the plots in the IOP article (2016).
# BLB 2016
# Working as of 30/08/2016.
#--------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------
# Init
#------------------------------------------------
source("source/init.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li    = "L1"
MODEL = "QBCP"
FWRK  = "EM"

#------------------------------------------------
# Filename
#------------------------------------------------
sizep = "1.06534_test"
filename = "~/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L1/orbit_order_20_size_-0.145455.txt"

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
muR = muR(FWRK);
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK)
Period = ifelse(MODEL=="QBCP", 6.79119387190792, 2*pi)

#------------------------------------------------
# Type of plot (can be changed at will)
#------------------------------------------------
if(MODEL == "QBCP") 
{
  fplot = plotdf_line;
}else
{
  fplot = plotdf_line_dashed;
}
fplot_path = plotdf_path;

#------------------------------------------------
# Strings
#------------------------------------------------
ts = "$t$ (\\%$T$)"

#--------------------------------------------------------------------------------------------------------------------------
#                                             XYZ
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "XYZ_Order_"
# Concatenate the results
orbit = read.table(filename, header = F)
names(orbit) = c("t", "x", "y", "z")

#------------------------------------------------
# Post-processing on coordinates and units
#------------------------------------------------
# From NC coordinates to C coordinates
orbit    = NCtoC(orbit, gamma)
orbit    = CtoPH(orbit, L)
orbit    = NCtoSYS(orbit, gamma, c1)
orbit    = SYStoPH(orbit, L)

#------------------------------------------------
# Kee only up to a certain time
#------------------------------------------------
orbit =  orbit[which(orbit$t < 50*Period),]

#-------------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------------
pxy = fplot_path (orbit,  "xPH", "yPH", "$X$ [km]",  "$Y$ [km]", lineSize = 0.3)
#Scaling
#pxy = pxy + scale_x_continuous(limits = c(-4.75e5, -4.25e5))  #cont scale on the x axis 
#pxy = pxy + scale_y_continuous(limits = c(-4e4,  4e4))   #cont scale on the y axis 
#Theme
pxy = pxy + theme_bw()+custom_theme
#Aspect ratio
#pxy = pxy+ coord_fixed(ratio=1)
pxy

#-------------------------------------------------------------------------------------
# Plot (yz)
#-------------------------------------------------------------------------------------
pyz = fplot_path(orbit,  "yPH", "zPH", "$Y$ [km]", "$Z$ [km]", lineSize = 0.3)
#Scaling
#pyz = pyz + scale_x_continuous(limits = c(-5e4, 5e4))  #cont scale on the x axis 
#pyz = pyz + scale_y_continuous(limits = c(-10000, 13000))  #cont scale on the y axis 
#Theme
pyz = pyz + theme_bw()+custom_theme
#Aspect ratio
pyz = pyz#+ coord_fixed(ratio=1)
pyz


#-------------------------------------------------------------------------------------
# Plot (xz)
#-------------------------------------------------------------------------------------
pxz = fplot_path (orbit,  "xPH", "zPH", "$X$ [km]", "$Z$ [km]", lineSize = 0.3)
#Scaling
#pxz = pxz + scale_x_continuous(limits = c(-50, 40))  #cont scale on the x axis 
#pxz = pxz + scale_y_continuous(limits = c(-3500, 3500))  #cont scale on the y axis 
#Theme
pxz = pxz + theme_bw()+custom_theme
#Aspect ratio
pxz = pxz#+ coord_fixed(ratio=1)
pxz

#--------------------------------------------------------------------------------------------------------------------------
#Save in pdf file
#--------------------------------------------------------------------------------------------------------------------------
ggsave(pxy, file = paste0(currentfolder,  "R_XY_Size_", sizep, ".pdf"))
ggsave(pyz, file = paste0(currentfolder,  "R_YZ_Size_", sizep, ".pdf"))
ggsave(pxz, file = paste0(currentfolder,  "R_XZ_Size_", sizep, ".pdf"))

#--------------------------------------------------------------------------------------------------------------------------
#Save in tex file
#--------------------------------------------------------------------------------------------------------------------------
ggplot2tikz(pxy, xSize, ySize, file = paste0(currentfolder,  "R_XY_Size_", sizep, ".tex"))
ggplot2tikz(pyz, xSize, ySize, file = paste0(currentfolder,  "R_YZ_Size_", sizep, ".tex"))
ggplot2tikz(pxz, xSize, ySize, file = paste0(currentfolder,  "R_XZ_Size_", sizep, ".tex"))

#--------------------------------------------------------------------------------------------------------------------------
#3D
#--------------------------------------------------------------------------------------------------------------------------
lines3D(orbit$x, orbit$y, orbit$z, col = "black")



