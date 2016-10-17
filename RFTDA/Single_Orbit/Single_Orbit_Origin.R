# Script to plot the dynamical equivalent to the libration points of the Sun-Earth-Moon system.
# This file makes use of the date files of type (example) "DYNEQ of EML1 in EM.txt" contains in the folder "OOFTDA/plot/QBCP/..."
# And the resulting plots are saved in the same folder.
#
# Output files include: eps, pdf, tikz (.tex)
#
# The possible plots are:
#---------------------------------
# DYNEQ of EML1 in EM.txt, using:
# Li    = "L1"
# MODEL = "QBCP"
# FWRK  = "EM"
# FWRK2  = "SEM"
#---------------------------------
# DYNEQ of EML2 in EM.txt, using:
# Li    = "L2"
# MODEL = "QBCP"
# FWRK  = "EM"
# FWRK2  = "SEM"
#---------------------------------
# DYNEQ of SEML1 in SEM.txt, using:
# Li    = "L1"
# MODEL = "QBCP"
# FWRK  = "SEM"
# FWRK2  = "EM"
#---------------------------------
# DYNEQ of SEML1 in SEM.txt, using:
# Li    = "L2"
# MODEL = "QBCP"
# FWRK  = "SEM"
# FWRK2  = "EM"
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
FWRK2  = "EM"

#Working folder
currentfolder = paste0(plotfolder(MODEL, FWRK, Li), "orbits/")

#Period of the system
Period = ifelse(MODEL=="QBCP", 6.79119387190792, 2*pi)

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK);

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
# Load data
#------------------------------------------------
name = paste0("DYNEQ of ", Li, " in ", FWRK)
native_orbit = read.table(paste0(currentfolder,name,".txt"), header = F)
colnames(native_orbit) = c("t", "x", "y", "z")
native_orbit$type = paste0("Dyn. eq. of ", Li, " in ", FWRK, " ")

name = paste0("DYNEQ of ", Li, " prop. in ", FWRK2, " back in ", FWRK)
back_orbit = read.table(paste0(currentfolder, name,".txt"), header = F)
colnames(back_orbit) = c("t", "x", "y", "z")
back_orbit$type = paste0("Dyn. eq. of ", Li, " in ", FWRK2, " ")


#Rbind the two
orbit = rbind(native_orbit)#, back_orbit)


#------------------------------------------------
# Post-processing on coordinates and units
#------------------------------------------------
# From NC to EM units
orbit = NCtoSYS(orbit, gamma, c1)
# From EM to physical units
orbit = SYStoPH(orbit, L)

#------------------------------------------------
#Origin
#------------------------------------------------
Ldf = data.frame(x = 0, y = 0, z =0)
# From NC to EM units
Ldf = NCtoSYS(Ldf, gamma, c1)
# From EM to physical units
Ldf = SYStoPH(Ldf, L)


#------------------------------------------------
#Select half period time
#------------------------------------------------
native_orbit_half = ddply(orbit, ~type, function(x){x[which.min(abs(x$t-0.5*Period)),]})

#------------------------------------------------
#Center manifold
#------------------------------------------------
native_orbit_start = orbit[which(orbit$t == 0.0),]

#-------------------------------------------------------------------------------------
# Plot 1: computed in the native framework (EM or SEM)
#-------------------------------------------------------------------------------------
#Path
#------------------------------------------------
porb = fplot_path (orbit,  "xEM", "yEM",  "$X$ $[$-$]$",  "$Y$ $[$-$]$", "type", "", TRUE, "type")

#------------------------------------------------
# Starting point
#------------------------------------------------
porb = porb + geom_point(data = native_orbit_start, aes(xEM, yEM), size = 4, colour = "black", fill = "black", pch = 21)

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
porb = porb+scale_color_manual(guide = FALSE, values=c("black"))

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
  porb = porb + scale_x_continuous(breaks = c(-0.9899842, -0.9899841, -0.9899840))  #cont scale on the x axis
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
  porb = porb + scale_y_continuous(limits = c(-9.17e-8, 9.17e-8), breaks = seq(-1e-7, 1e-7, 5e-8), labels=scientific_10)
}

if(Li == "L1" && FWRK == "SEM")
{
  porb = porb + scale_y_continuous(limits = c(-9.17e-8, 9.17e-8), breaks = seq(-1e-7, 1e-7, 5e-8), labels=scientific_10)
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
ggsave(porb, width = xSize, height = ySize, file = paste0(currentfolder, "orbit_", Li, ".eps"))
ggsave(porb, width = xSize, height = ySize, file = paste0(currentfolder, "orbit_", Li, ".pdf"))

#--------------------------------------------------------------------------------------------------------------------------
#Save in tikz
#--------------------------------------------------------------------------------------------------------------------------
# Does not currently work
ggplot2tikz(porb, xSize, ySize, file = paste0(currentfolder, "orbit_", Li, ".tex"))

#--------------------------------------------------------------------------------------------------------------------------
# Old code that is now included in ggplot2tikz
#--------------------------------------------------------------------------------------------------------------------------
# #Create a .tex file that will contain your plot as vectors
# #You need to set the size of your plot here, if you do it in LaTeX, font consistency with the rest of the document will be lost
# tikz(paste0(currentfolder, "orbit_", Li, ".tex"), width = xSize, height = ySize, standAlone=TRUE, documentDeclaration="\\documentclass{standalone}\n",
#      packages = c("\\usepackage[utf8]{inputenc}",
#                   "\\usepackage[T1]{fontenc}",
#                   "\\usepackage{tikz}", 
#                   "\\usepackage{pgf}", 
#                   "\\usetikzlibrary{calc}", 
#                   "\\usepackage{amssymb}", 
#                   "\\usepackage{amsfonts}\n"))
# 
# #Simple plot of the dummy data using LaTeX elements
# plot <- porb
# #This line is only necessary if you want to preview the plot right after compiling
# print(plot)
# #Necessary to close or the tikxDevice .tex file will not be written
# dev.off()
