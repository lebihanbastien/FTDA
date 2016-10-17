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
currentfolder = paste0(plotfolder(MODEL, FWRK, Li), "Serv/")
size = "-1.06534_test"
sizep = "1.06534_test"

# Orders
dfindex  = c(8,10,15,20)
dfnumber = 1:length(dfindex)
maxOrder = max(dfindex)


#------------------------------------------------
# Parameters
#------------------------------------------------
maxPrec = 5e-7

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
ts = "$t$ [\\%$T$]"

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
           ts,  "$e_H$ [-]",
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

#Posprocessing
#------------------------------------------------
eI$EI = eI$y*gamma;

#Plot
#------------------------------------------------
pI = fplot(eI,  "x", "EI", ts,  "$E_I$ $[$-$]$", "order", "Order", 1)
#Cont scale on the x axis 
pI = pI + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), 
                             labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"))  #cont scale on the x axis 
#Logscale on the y axis
pI = pI + scale_y_log10(breaks = 10^(-8:0*2))
#Theme
pI = pI+legend_inside_theme
#Display
pI


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
           ts,  "$E_O$ $[$-$]$",
           "order", "Order", 1)
#Cont scale on the x axis 
pO = pO + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"))  #cont scale on the x axis 
#Logscale on the y axis
pO = pO + scale_y_log10(limits=c(1e-10,NaN), breaks = 10^(-8:0*2))
#Theme
pO = pO+legend_inside_theme
#Display
pO

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
# From NC coordinates to C coordinates
orbit    = NCtoC(orbit, gamma)
orbit    = CtoPH(orbit, L)
orbit    = NCtoSYS(orbit, gamma, c1)
orbit    = SYStoPH(orbit, L)
orbit_pm = NCtoC(orbit_pm, gamma)
orbit_pm = CtoPH(orbit_pm, L)
orbit_pm = NCtoSYS(orbit_pm, gamma, c1)
orbit_pm = SYStoPH(orbit_pm, L)

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
orbitMaxOrder_end = orbitMaxOrder[which(orbitMaxOrder$t == max(orbitMaxOrder$t)),]

#-------------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------------
pxy = fplot_path (orbit,  "xPH", "yPH", 
                   "$X$ [km]",  "$Y$ [km]",
                   "order", "Order", 1 )#, "order")
#Grey scale if needed
#pxy = pxy + scale_colour_grey(start = 0.9, end = 0.0, name="Order")
#Center manifold plot
pxy = pxy + geom_path(data  = orbitMaxOrder, aes(xPH, yPH, color = "PM"), size = linesize[["line"]], linetype = "dotted", colour = "black")
# Starting point
pxy = pxy + geom_point(data = orbitMaxOrder_start, aes(xPH, yPH), size = 4, colour = "black", fill = "white", pch = 21)
# Point at a given precision
pxy = pxy + geom_point(data = orbitMaxOrder_end, aes(xPH, yPH), size = 4, pch = 22, colour = "black", fill = "white")
#Scaling
pxy = pxy + scale_x_continuous(limits = c(-4.75e5, -4.25e5))  #cont scale on the x axis 
pxy = pxy + scale_y_continuous(limits = c(-4e4,  4e4))   #cont scale on the y axis 
#Theme
pxy = pxy + theme_bw()+custom_theme
#Annotation
#pxy = pxy + annotate("text", x = -0.3, y = -0.21, label = "CM \n Order maxOrder", parse = F, size = 8, colour = "limegreen", fontface = "bold")
#pxy = pxy + geom_segment(aes(x = -0.25, y = -0.2, xend = -0.19, yend = -0.2), arrow = arrow(length = unit(0.5, "cm")), size = psize)
#Add an arrow to give the direction of motion along the orbit
ai= 10
pxy = pxy + geom_segment(aes(x = orbit$xPH[ai], y = orbit$yPH[ai], xend = orbit$xPH[ai+1], yend = orbit$yPH[ai+1]), 
                           colour = muted("blue"), 
                           arrow = arrow(length = unit(0.4, "cm"), type = "closed"))
pxy = pxy+legend_inside_theme#+theme(legend.position = c(0.5,0.5))


#Add the Moon
#-------------------
moonR = 1737.10       #Moon's radius
moonPos = gamma*L     #Moon position in km wrt to Li
pxy = addMoon(pxy, x = moonPos, y = 0, moonR, surfSize = 0.4, cratSize = 0.2)+ coord_fixed()


#-------------------------------------------------------------------------------------
# Plot (yz)
#-------------------------------------------------------------------------------------
pyz = fplot_path (orbit,  "yPH", "zPH", 
                   "\ny [km]",  "z [km]\n",
                   "order", "Order", 1 )#, "order")
#Grey scale if needed
#pyz = pyz + scale_colour_grey(start = 0.9, end = 0.0, name="Order")
#Center manifold plot
pyz = pyz + geom_path(data  = orbitMaxOrder, aes(yPH, zPH, color = "PM"), size = linesize[["line"]], linetype = "dashed", colour = "black")

# Starting point
pyz = pyz + geom_point(data = orbitMaxOrder_start, aes(yPH, zPH), size = 4, colour = "black", fill = "white", pch = 21)

# Point at a given precision
#pyz = pyz + geom_point(data = orbitMaxOrder_prec, aes(yPH, zPH), size = 4, pch = 22, colour = "black", fill = "white")

#Scaling
pyz = pyz + scale_x_continuous(limits = c(-5e4, 5e4))  #cont scale on the x axis 
pyz = pyz + scale_y_continuous(limits = c(-10000, 13000))  #cont scale on the y axis 
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


#-------------------------------------------------------------------------------------
# Plot (xz)
#-------------------------------------------------------------------------------------
pxz = fplot_path (orbit,  "xPH", "zPH", 
                  "\ny [km]",  "z [km]\n",
                  "order", "Order", 1 )#, "order")
#Grey scale if needed
pxz = pxz + scale_colour_grey(start = 0.9, end = 0.0, name="Order")
#Center manifold plot
pxz = pxz + geom_path(data  = orbitMaxOrder, aes(xPH, zPH, color = "PM"), size = linesize[["line"]], linetype = "dashed", colour = "black")
# Starting point
pxz = pxz + geom_point(data = orbitMaxOrder_start, aes(xPH, zPH), size = 4, colour = "black", fill = "white", pch = 21)
# Point at a given precision
#pxz = pxz + geom_point(data = orbitMaxOrder_prec, aes(xPH, zPH), size = 4, pch = 22, colour = "black", fill = "white")

#Scaling
#pxz = pxz + scale_x_continuous(limits = c(-50, 40))  #cont scale on the x axis 
#pxz = pxz + scale_y_continuous(limits = c(-3500, 3500))  #cont scale on the y axis 
#Theme
pxz = pxz + theme_bw()+custom_theme
#Annotation
#pxz = pxz + annotate("text", x = -0.3, y = -0.21, label = "CM \n Order maxOrder", parse = F, size = 8, colour = "limegreen", fontface = "bold")
#pxz = pxz + geom_segment(aes(x = -0.25, y = -0.2, xend = -0.19, yend = -0.2), arrow = arrow(length = unit(0.5, "cm")), size = psize)
#Display the plot
pxz
pxz = pxz + geom_segment(aes(x = orbit$xPH[ai], y = orbit$zPH[ai], xend = orbit$xPH[ai+1], yend = orbit$zPH[ai+1]), 
                         colour = muted("blue"), 
                         arrow = arrow(length = unit(0.4, "cm"), type = "closed"))
pxz = pxz+legend_inside_theme

#--------------------------------------------------------------------------------------------------------------------------
#Save in pdf file
#--------------------------------------------------------------------------------------------------------------------------
ggsave(pH,   file = paste0(currentfolder, "R_eH_Size_", sizep,  ".pdf"))
ggsave(pI,   file = paste0(currentfolder, "R_eI_Size_", sizep,  ".pdf"))
ggsave(pO,   file = paste0(currentfolder, "R_eO_Size_", sizep,  ".pdf"))
ggsave(pxy,  file = paste0(currentfolder, "R_XY_Size_", sizep, ".pdf"))
ggsave(pyz,  file = paste0(currentfolder, "R_YZ_Size_", sizep, ".pdf"))

#--------------------------------------------------------------------------------------------------------------------------
#Save in tex file
#--------------------------------------------------------------------------------------------------------------------------
ggplot2tikz(pI, xSize, ySize, file = paste0(currentfolder, "R_eI_Size_", sizep,  ".tex"))
ggplot2tikz(pO, xSize, ySize, file = paste0(currentfolder, "R_eO_Size_", sizep,  ".tex"))