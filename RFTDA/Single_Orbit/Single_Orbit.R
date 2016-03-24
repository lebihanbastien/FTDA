# Script to plot an orbit along with its associated precisions
# in the parameterization method of the QBFBP/RTBP around L1/L2 of the Earth-Moon system
#----------------------------------------------------------------------------------------

#------------------------------------------------
# Load libraries
#------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(png)
library(grImport)

#------------------------------------------------
# Load Source files
#------------------------------------------------
source("source/source_folder.R")
source("source/source_plot.R")
source("source/source_routines.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li    = "L2"
MODEL = "QBCP"
FWRK  = "SEM"
currentfolder = paste0(printfolder(MODEL, FWRK, Li), "orbits/")
size = "30_s1s2" #1.5 7.5 15 30 45
Period = ifelse(MODEL=="QBCP", 6.79119387190792, 2*pi)
maxPrec = 1e-6

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

# Orders
if(MODEL == "QBCP") dfindex = c(10,15,20) else dfindex = c(2,5,20)
dfnumber = 1:length(dfindex)
maxOrder = max(dfindex)

#--------------------------------------------------------------------------------------------------------------------------
#                                             eH
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "eH_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0(currentfolder, dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("x", "y")
  dflist[[i]]$order = dfindex[i]
}
# Concatenate the results
eH = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){eH = rbind(eH, dflist[[i]])}


#Plot
#------------------------------------------------
#Actual plot
pH = fplot(eH,  "x", "y", 
          "\n time (%T)",  "eH [-] \n",
          "order", "Order", 1)
#Cont scale on the x axis 
pH = pH + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), 
                             labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"))  #cont scale on the x axis 
#Logscale on the y axis
if(MODEL == "QBCP") pH = pH + scale_y_log10(limits=c(1e-12, 1e3), breaks = 10^(-8:0*2))else 
  pH = pH + scale_y_log10(limits=c(1e-15, 1e-1), breaks = 10^(-8:0*2)) 
#Display the plot
pH = pH+legend_inside_theme

#--------------------------------------------------------------------------------------------------------------------------
#                                             eI
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "eI_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0(currentfolder, dftype, dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("x", "y")
  dflist[[i]]$order = dfindex[i]
}
# Concatenate the results
eI = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){eI = rbind(eI, dflist[[i]])}


#Plot
#------------------------------------------------
pI = fplot(eI,  "x", "y", 
           "\n time (%T)",  "eI [-] \n",
           "order", "Order", 1)
#Cont scale on the x axis 
pI = pI + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), 
                             labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"))  #cont scale on the x axis 
#Logscale on the y axis
pI = pI + scale_y_log10(limits=c(1e-15, 5e-11), breaks = 10^(-8:0*2))
#Display
pI = pI+legend_inside_theme


#--------------------------------------------------------------------------------------------------------------------------
#                                             eO
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "eO_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0(currentfolder, dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("x", "y")
  dflist[[i]]$order = dfindex[i]
}
# Concatenate the results
eO = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){eO = rbind(eO, dflist[[i]])}

#Plot
#------------------------------------------------
pO = fplot(eO,  "x", "y", 
           "\n time (%T)",  "eO [-] \n",
           "order", "Order", 1)
#Cont scale on the x axis 
pO = pO + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), 
                             labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"))  #cont scale on the x axis 
#Logscale on the y axis
pO = pO + scale_y_log10(limits=c(1e-8,2e1), breaks = 10^(-8:0*2))
#Display
pO = pO+legend_inside_theme


#--------------------------------------------------------------------------------------------------------------------------
#                                             XYZ
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "XYZ_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0(currentfolder, dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("t", "x", "y", "z")
  dflist[[i]]$order = dfindex[i]
}

# Concatenate the results
orbit = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){orbit = rbind(orbit, dflist[[i]])}

#------------------------------------------------
# Load data (PM)
#------------------------------------------------
dflist = list();
dftype = "XYZ_PM_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0(currentfolder, dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("t", "x", "y", "z")
  dflist[[i]]$order = dfindex[i]
}

# Concatenate the results
orbit_pm = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){orbit_pm = rbind(orbit_pm, dflist[[i]])}

#------------------------------------------------
# Post-processing on coordinates and units
#------------------------------------------------
# F rom NC coordinates to C coordinates
orbit    = NCtoC(orbit, gamma)
orbit    = CtoPH(orbit, L)
orbit_pm = NCtoC(orbit_pm, gamma)
orbit_pm = CtoPH(orbit_pm, L)

#------------------------------------------------
#Select half period time
#------------------------------------------------
orbit_half = ddply(orbit, ~order, function(x){x[which.min(abs(x$t-0.25*Period)),]})

#------------------------------------------------
#Select when orbital precision is > maxPrec
#------------------------------------------------
eO_prec = ddply(eO, ~order, function(x){x[which.min(abs(x$y-maxPrec)),]})  #time for which eO ~ maxPrec
orbit_prec = ddply(orbit, ~order, function(orb){orb[which.min(abs(orb$t-eO_prec$x[which(eO_prec$order == maxOrder)])),]}) #select order maxOrder

#------------------------------------------------
#Center manifold
#------------------------------------------------
orbitMaxOrder = orbit_pm[which(orbit_pm$order == maxOrder),]
orbitMaxOrder_prec = orbit_prec[which(orbit_prec$order==maxOrder),]
orbitMaxOrder_start = orbitMaxOrder[which(orbitMaxOrder$t == 0.0),]

#-------------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------------
porb = fplot_path (orbit,  "xPH", "yPH", 
           "\nx [km]",  "y [km]\n",
           "order", "Order", 1 )#, "order")
#Grey scale if needed
#porb = porb + scale_colour_grey(start = 0.9, end = 0.0, name="Order")
#Center manifold plot
porb = porb + geom_path(data  = orbitMaxOrder, aes(xPH, yPH, color = "PM"), size = linesize[["line"]], linetype = "dotted", colour = "black")
# Starting point
porb = porb + geom_point(data = orbitMaxOrder_start, aes(xPH, yPH), size = 4, colour = "black", fill = "white", pch = 21)
# Point at a given precision
porb = porb + geom_point(data = orbitMaxOrder_prec, aes(xPH, yPH), size = 4, pch = 22, colour = "black", fill = "white", solid = T)
#Scaling
#porb = porb + scale_x_continuous(limits = c(-5e4, 10e4))  #cont scale on the x axis 
#porb = porb + scale_y_continuous(limits = c(-5e4, 5e4))  #cont scale on the y axis 
#Theme
porb = porb + theme_bw()+custom_theme
#Annotation
#porb = porb + annotate("text", x = -0.3, y = -0.21, label = "CM \n Order maxOrder", parse = F, size = 8, colour = "limegreen", fontface = "bold")
#porb = porb + geom_segment(aes(x = -0.25, y = -0.2, xend = -0.19, yend = -0.2), arrow = arrow(length = unit(0.5, "cm")), size = psize)
#Add an arrow to give the direction of motion along the orbit
ai= 10
porb = porb + geom_segment(aes(x = orbit$xPH[ai], y = orbit$yPH[ai], xend = orbit$xPH[ai+1], yend = orbit$yPH[ai+1]), 
                           colour = muted("blue"), 
                           arrow = arrow(length = unit(0.4, "cm"), type = "closed"))
porb = porb+legend_inside_theme#+theme(legend.position = c(0.5,0.5))


#Add the Moon
#-------------------
moonR = 1737.10       #Moon's radius
moonPos = gamma*L #Moon position in km wrt to Li
porb = addMoon(porb, x = moonPos, y = 0, moonR, surfSize = 0.4, cratSize = 0.2)+ coord_fixed()


#-------------------------------------------------------------------------------------
# Plot (yz)
#-------------------------------------------------------------------------------------
pyz = fplot_path (orbit,  "yPH", "zPH", 
                   "\ny [km]",  "z [km]\n",
                   "order", "Order", 1 )#, "order")
#Grey scale if needed
#pyz = pyz + scale_colour_grey(start = 0.9, end = 0.0, name="Order")
#Center manifold plot
pyz = pyz + geom_path(data  = orbitMaxOrder, aes(yPH, zPH, color = "PM"), size = linesize[["line"]], linetype = "dotted", colour = "black")
# Starting point
pyz = pyz + geom_point(data = orbitMaxOrder_start, aes(yPH, zPH), size = 4, colour = "black", fill = "white", pch = 21)
# Point at a given precision
pyz = pyz + geom_point(data = orbitMaxOrder_prec, aes(yPH, zPH), size = 4, pch = 22, colour = "black", fill = "white", solid = T)
#Scaling
#pyz = pyz + scale_x_continuous(limits = c(-50, 40))  #cont scale on the x axis 
#pyz = pyz + scale_y_continuous(limits = c(-3500, 3500))  #cont scale on the y axis 
#Theme
pyz = pyz + theme_bw()+custom_theme
#Annotation
#pyz = pyz + annotate("text", x = -0.3, y = -0.21, label = "CM \n Order maxOrder", parse = F, size = 8, colour = "limegreen", fontface = "bold")
#pyz = pyz + geom_segment(aes(x = -0.25, y = -0.2, xend = -0.19, yend = -0.2), arrow = arrow(length = unit(0.5, "cm")), size = psize)
#Display the plot
pyz
pyz = pyz + geom_segment(aes(x = orbit$yPH[ai], y = orbit$zPH[ai], xend = orbit$yPH[ai+1], yend = orbit$zPH[ai+1]), 
                         colour = muted("blue"), 
                         arrow = arrow(length = unit(0.4, "cm"), type = "closed"))
pyz = pyz+legend_inside_theme


#--------------------------------------------------------------------------------------------------------------------------
#Save in eps file
#--------------------------------------------------------------------------------------------------------------------------
ggsave(pH, file = paste0(currentfolder, "R_eH_Size_", size, ".eps"))
ggsave(pI, file = paste0(currentfolder, "R_eI_Size_", size, ".eps"))
ggsave(pO, file = paste0(currentfolder, "R_eO_Size_", size, ".eps"))
ggsave(porb, file = paste0(currentfolder,  "R_XY_Size_", size, ".eps"))
ggsave(pyz,  file = paste0(currentfolder,  "R_YZ_Size_", size, ".eps"))

#--------------------------------------------------------------------------------------------------------------------------
#Initialize gnuplot window
#--------------------------------------------------------------------------------------------------------------------------
# h1<-Gpinit()
# Gpresetplot(h1)
# #Gpcmd(h1,'set size ratio 2')
# #Labels
# Gpcmd(h1, 'set xlabel \"$x$\"')
# Gpcmd(h1, 'set ylabel \"$y$\"')
# Gpcmd(h1, 'set zlabel \"$z$\"')
# #Grid
# Gpcmd(h1, 'set grid')
# Gpcmd(h1, 'set grid xtics')
# Gpcmd(h1, 'set grid ytics')
# Gpcmd(h1, 'set grid ztics')
# #Mesh
# Gpcmd(h1, ' set dgrid3d 60,60 qnorm 2')
# Gpcmd(h1, 'set hidden3d')
# #Range
# # Gpcmd(h1, 'set xrange [-2:42];
# #            set yrange [-2:12];
# #            set zrange [-0:50];')
# #Position of the z-axis
# Gpcmd(h1, 'set ticslevel 0')
# #Borders
# # Gpcmd(h1, 'set border 4095')
# Gpcmd(h1, 'set border 127+256+512')
# #set border 127+256+512 # or set border 1023-128
# #Gnuplot option
# replot = F;
# color = 1
# 
# 
# #------------------------------------------------
# # 3D plot with RGnuplot
# #------------------------------------------------
# #Build temp data
# fileName = paste0("test_orbit", ".data")
# ttm_plot = data.frame(orbit$yPH, orbit$zPH, orbit$zPH);
# write.table(ttm_plot, fileName , sep = " ", col.names = F, row.names = F)
# #Actual plot
# Gpcmd(h1, CGpsetcmd(plotType = "plot", fileName = fileName, lt = "fc", ls = "points", lc = color, replot = F));

