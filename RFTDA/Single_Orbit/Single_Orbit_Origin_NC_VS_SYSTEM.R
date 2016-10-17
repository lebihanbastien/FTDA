# Script to plot the dynamical equivalent to the libration points of the Sun-Earth-Moon system.
# This file plots on a single figure the dyneq computed in NCSYS (NCEM or NCSEM) + the dyneq computed in SYS (EM or SEM)
# This file makes use of the date files of type (example) "DYNEQ_QBCP_EM_L1.txt" contains in the folder "OOFTDA/plot/QBCP/DYNEQ/..."
# And the resulting plots are saved in the same folder.
#
# Output files include: eps, pdf, tikz (.tex)
#
# The possible plots are:
#---------------------------------
# DYNEQ_QBCP_EM_L1.txt, using:
# Li    = "L1"
# MODEL = "QBCP"
# FWRK  = "EM"
# FWRK2  = "SEM"
#---------------------------------
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
Li    = "L2"
MODEL = "QBCP"
FWRK  = "SEM"

#Working folder
currentfolder = paste0(ooftdafolder, "plot/QBCP/DYNEQ/")

#Period of the system
Period = ifelse(MODEL=="QBCP", 6.79119387190792, 2*pi)

#------------------------------------------------
#Select the libration point
#------------------------------------------------
LIB = WHICH_LIB_POINT(Li, FWRK);


#------------------------------------------------
# Type of plot (can be changed at will)
#------------------------------------------------
if(MODEL == "QBCP") 
{
  fplot      = plotdf_line;
}else
{
  fplot = plotdf_line_dashed;
}
fplot_path = plotdf_path;

#--------------------------------------------------------------------------------------------------------------------------
#                                             XYZ
#--------------------------------------------------------------------------------------------------------------------------
# Load data.
# Note: the data are already in SYS units (EM or SEM)
# Contrary to the case Single_Orbit_Origin.R
#------------------------------------------------
# Second data: dynamical equivalent computed in NCSYS, back in SYS
name = paste0("DYNEQ_", MODEL, "_", FWRK, "_", Li, "_NC")
nc_orbit = read.table(paste0(currentfolder,name,".txt"), header = F)
colnames(nc_orbit) = c("xEM", "yEM")
nc_orbit$type = name

# First data: dynamical equivalent computed in NCSYS, back in SYS
name = paste0("DYNEQ_", MODEL, "_", FWRK, "_", Li)
sys_orbit = read.table(paste0(currentfolder,name,".txt"), header = F)
colnames(sys_orbit) = c("xEM", "yEM")
sys_orbit$type = paste0(name, "_SYS")

#Rbind the two
orbit = rbind(nc_orbit, sys_orbit)

#------------------------------------------------
#Origin
#------------------------------------------------
Ldf = data.frame(x = 0, y = 0, z =0)
# From NC to EM units
Ldf = NCtoSYS(Ldf, LIB$GAMMA, LIB$C1)
# From EM to physical units
Ldf = SYStoPH(Ldf, LIB$LPRIM)

#------------------------------------------------
#Center manifold
#------------------------------------------------
nc_orbit_start = orbit[1,]

#-------------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------------
#Path
#------------------------------------------------
porb = fplot_path (orbit,  "xEM", "yEM",  "$X$ $[$-$]$",  "$Y$ $[$-$]$", "type", "", TRUE, "type")

#------------------------------------------------
# Starting point
#------------------------------------------------
porb = porb + geom_point(data = nc_orbit_start, aes(xEM, yEM), size = 4, colour = "black", fill = "black", pch = 21)

#------------------------------------------------
# Origin point
#------------------------------------------------
if(FWRK == "EM")
{
  porb = porb + geom_point(data = Ldf, aes(xEM, yEM), size = 4, colour = "black", fill = "white", pch = 21)
}

#------------------------------------------------
#Theme
#------------------------------------------------
porb = porb+custom_theme

#Add an arrow to give the direction of motion along the orbit
ai= 100
porb = porb + geom_segment(aes(x = orbit$xEM[ai], y = orbit$yEM[ai], xend = orbit$xEM[ai+1], yend = orbit$yEM[ai+1]), 
                           colour = "black", 
                           arrow = arrow(length = unit(0.4, "cm"), type = "closed"))

#------------------------------------------------
#Add annotation
#------------------------------------------------
if(Li == "L2" && FWRK == "EM")
{
  porb = porb + annotate("text", x = Ldf$xEM+1e-07, y = Ldf$yEM+3e-07, label = "$L_2$", size = 10)#, parse =TRUE) 
}
if(Li == "L1" && FWRK == "EM"){
  porb = porb + annotate("text", x = Ldf$xEM+1e-07, y = Ldf$yEM+3e-07, label = "$L_1$", size = 10)#, parse =TRUE)
}

#------------------------------------------------
#Color scale
#------------------------------------------------
#Color
#porb = porb+scale_color_discrete(guide = FALSE)
porb = porb+scale_color_manual(guide = FALSE, values=c("black", "#FF6666"))
#Linetype
porb = porb+scale_linetype_manual(guide = FALSE, values=c("solid", "dashed"))

#------------------------------------------------
#Scaling, if necessary
#------------------------------------------------
if(Li == "L1" && FWRK == "EM")
{
  #porb = porb + scale_x_continuous(breaks = c(-0.980075594, -0.980075593))  #cont scale on the x axis
}

if(Li == "L2" && FWRK == "EM"){
  porb = porb + scale_x_continuous(breaks = c(-1.155683, -1.155682))  #cont scale on the x axis
}

if(Li == "L1" && FWRK == "SEM")
{
  #porb = porb + scale_x_continuous(breaks = c(-0.980075594, -0.980075593))  #cont scale on the x axis
}

if(Li == "L2" && FWRK == "SEM"){
  porb = porb + scale_x_continuous(breaks = c(-1.01007730, -1.01007720, -1.01007710))  #cont scale on the x axis
}

#------------------------------------------------
#Y labels
#------------------------------------------------
scientific_10 <- function(x) {
  l = gsub("0e\\+00","0", x);
  l = gsub("e\\-0", "e\\-", l)
  l = gsub("e\\+0", "e\\+", l)
  text=l
}
if(FWRK == "EM")
{
  porb = porb + scale_y_continuous(breaks = seq(-3e-6, 3e-6, 1e-6), labels=scientific_10)
}
if(Li == "L2" && FWRK == "SEM")
{
  porb = porb + scale_y_continuous(breaks = seq(-1e-7, 1e-7, 5e-8), labels=scientific_10)
}

if(Li == "L1" && FWRK == "SEM")
{
  #porb = porb + scale_y_continuous(breaks = seq(-1e-9, 1e-9, 5e-10), labels=scientific_10)
}

#------------------------------------------------
# Aspect ratio, if necessary
#------------------------------------------------
if(FWRK == "SEM")
{
  #Aspect ratio
  porb = porb+ coord_fixed(ratio=1)
  #Change the size
  xSize = 10.4
  ySize = 7.22
}


#--------------------------------------------------------------------------------------------------------------------------
#Save in eps/pdf file
#--------------------------------------------------------------------------------------------------------------------------
ggsave(porb, width = xSize, height = ySize, file = paste0(currentfolder, name, ".eps"))
ggsave(porb, width = xSize, height = ySize, file = paste0(currentfolder, name, ".pdf"))

#--------------------------------------------------------------------------------------------------------------------------
#Save in tikz
#--------------------------------------------------------------------------------------------------------------------------
# Does not currently work
ggplot2tikz(porb, xSize, ySize, file = paste0(currentfolder, name, ".tex"))
