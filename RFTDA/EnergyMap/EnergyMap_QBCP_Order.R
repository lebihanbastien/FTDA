# R script to handle an energy map of the QBCP around EML1,2
#===============================================================================

#===============================================================================
# Init
#===============================================================================
# R options
#-------------------------------------------------------------------------------
options(digits = 15)


# Load libraries
#-------------------------------------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)

# Load Source files
#-------------------------------------------------------------------------------
source("source/init.R")


# Select Models & libration point
#-------------------------------------------------------------------------------
Li        = "L2"
MODEL     = "QBCP"
FWRK      = "EM"
Energy    = "0.1"
order     = "1"
ofs_order = "30"
t0        = "0.5";

#Current working folder
currentfolder = paste0(printfolder(MODEL, FWRK, Li))


#Normalized units (gamma, c1)
#-------------------------------------------------------------------------------
gamma = gamma(Li, FWRK);
L     = Ldist(FWRK);
c1    =  c1(Li, FWRK);

#===============================================================================
# Data reading
#===============================================================================

# Filename to check
#-------------------------------------------------------------------------------
fileprefix = paste0(currentfolder, "Serv/Serv_hm_Energy_", Energy, "_order_", order, "_ofs_", ofs_order, "_t0_", t0)
filename   = paste0(fileprefix, ".bin")

# Load bin source
#-------------------------------------------------------------------------------
if (file.exists(filename))
{
  names = c("label", "x", "y", "z", "px", "py", "pz",
            "xEM", "yEM", "zEM", "pxEM", "pyEM", "pzEM",
            "s1", "s2", "s3", "s4", "t", "Hz", "dHz");
  icdf = dffbinary(filename, 20, names);
  icdf$order = order;
  
}else
{
  icdf = data.frame()
}

# Postprocessing
#-------------------------------------------------------------------------------
#Get rid of the origin
isOrigin = icdf$s3 == icdf$s3[length(icdf$s3)] | icdf$s3 == icdf$s3[1] | icdf$s3 == 0 |
           icdf$s1 == icdf$s1[length(icdf$s1)] | icdf$s1 == icdf$s1[1] | icdf$s1 == 0  #& icdf$s2 == 0 & icdf$s3 == 0 & icdf$s4 == 0
icdf = icdf[which(isOrigin),]

# From NC to EM units
icdf = NCtoC(icdf, gamma)
  
# From EM to physical units
icdf = CtoPH(icdf, L)

# Radii from Li
icdf$rNC = sqrt(icdf$x^2+icdf$y^2+icdf$z^2)
icdf$rPH = sqrt(icdf$xCPH^2+icdf$yCPH^2+icdf$zCPH^2)

# Mean value of Hz
mHz = mean(icdf$Hz)

# Deviation wrt the mean value (signed)
icdf$ssdHz = icdf$Hz - mHz


#Select only 0 < dHz <= maxValue
#issdHzSuitable = abs(icdf$dHz) <= maxdHz # & icdf$dHz >= 0.0
#icdf = icdf[which(issdHzSuitable),]

#Select only s2 < 40 (heuristic, getting rid of )
#icdf = icdf[which(abs(icdf$s2) <= 40),]
#icdf = icdf[which(icdf$s3 >= 0),]

#===============================================================================
# Plots
#===============================================================================


# s3 vs s1
#-------------------------------------------------------------------------------
pT = ggplot()+geom_tile(data = icdf, aes(s1, s3, fill = ssdHz, colour = ssdHz))+custom_theme
pT = pT + scale_fill_gradient2(space="Lab", midpoint = 0, mid = "white", high = muted("blue"))
pT = pT + scale_colour_gradient2(space="Lab", midpoint = 0, mid = "white", high = muted("blue"), guide = FALSE)
pT


# x0 vs y0
#-------------------------------------------------------------------------------
pxy = plotdf_point(icdf, "x", "y", "x", "y", "ssdHz", "ssdHz", 0, pointSize = 3)
pxy = pxy + scale_colour_gradient2(space="Lab", midpoint = 0, mid = "white", high = muted("blue"), guide = FALSE)
pxy

#---------------
# dH vs s1
#---------------
plotdf_point(icdf, "s1", "dHz", "s1", "dHz", "s2", "s2", 0)+ scale_colour_gradient2(space="Lab", midpoint = 0, mid = "white", high = muted("blue"))


stop()

#===============================================================================
# Second Postprocessing & Plot with gnuplot (in progress, needs clean up)
#===============================================================================
#-----------------------
#Initialize gnuplot window
#-----------------------
h1<-Gpinit()
Gpresetplot(h1)
#Gpcmd(h1,'set size ratio 2')
#Labels
Gpcmd(h1, 'set xlabel \"$s_1$\"')
Gpcmd(h1, 'set ylabel \"$s_3$\"')
Gpcmd(h1, 'set zlabel \"$s_2$\"')
#Grid
Gpcmd(h1, 'set grid')
Gpcmd(h1, 'set grid xtics')
Gpcmd(h1, 'set grid ytics')
Gpcmd(h1, 'set grid ztics')
#Mesh
Gpcmd(h1, ' set dgrid3d 60,60 qnorm 2')
Gpcmd(h1, 'set hidden3d')
#Range
# Gpcmd(h1, 'set xrange [-2:42];
#            set yrange [-2:12];
#            set zrange [-0:50];')
#Position of the z-axis
Gpcmd(h1, 'set ticslevel 0')
#Borders
# Gpcmd(h1, 'set border 4095')
Gpcmd(h1, 'set border 127+256+512')
#set border 127+256+512 # or set border 1023-128
#Gnuplot option
replot = F;
color = 1

#-----------------------
#Tolerance for the energy
#-----------------------
dHz_vec = c(0.015)
tol = 1e-4;
#Loop on the energy
#---------------------
for(i in dHz_vec)
{
     #Select a given energy value, with a given tolerance
    dHz0 = c(i);
    #Select IC around this value
    icdf_dHz = icdf[which(abs(icdf$dHz - dHz0) <= tol),]
    #Select only positive values for s2
    isPos = icdf_dHz$s2 >= 0
    icdf_dHz = icdf_dHz[which(isPos),]
    
    #------------------------------------------------
    # 3D plot with RGnuplot
    #------------------------------------------------
    #Build temp data
    fileName = paste0("test_", i, ".data")
    ttm_plot = data.frame(icdf_dHz$s1, icdf_dHz$s3, icdf_dHz$s2);
    write.table(ttm_plot, fileName , sep = " ", col.names = F, row.names = F)
    
    #Actual plot
    Gpcmd(h1, CGpsetcmd(plotType = "splot", fileName = fileName, lt = "fc", ls = "points", lc = color, replot = replot));
    
    #For next step
    replot = T;
    color = color+1
}
Gpcmd(h1,'replot')

stop()

#Saving the plot 
#------------------------------------------------
# Gpcmd(h1,'set terminal postscript eps color;set output "testEPS.eps"')
Gpcmd(h1,'set terminal epslatex color;
          set output "testLATEX.eps"')
#Gpcmd(h1,'set terminal latex; set output "testLATEX.tex"')
Gpcmd(h1,'replot')

res = readString("Do you want to close the Gnuplot window? y/n: ")
if(res == "y") h1<-Gpclose(h1)

#Kill temp file
res = readString("Do you want to erase the temp file? y/n: ")
if(res == "y") file.remove(fileName);


# #------------------------------------------------
# # Plotly
# #------------------------------------------------
# py <- plotly(username="b.le-bihan", key="ado73xmcij")  # open plotly connection
# #------------------------------------------------
# # 3D plot
# #------------------------------------------------
# list1 = list(
#   x = icdf$s1,
#   y = icdf$s2,
#   z = icdf$s3,
#   mode = "markers",
#   marker = list(
#     color = icdf$dHz,
#     size = 2,
#     symbol = "circle",
#     line = list(
#       color = "rgb(204, 204, 204)",
#       width = 1
#     ),
#     opacity = 0.9
#   ),
#   type = "scatter3d"
# )
# 
# list2 = list(
#   x = ttm_dHz$s1,
#   y = ttm_dHz$s2,
#   z = ttm_dHz$s3,
#   mode = "markers",
#   marker = list(
#     color = "rgb(127, 127, 127)",
#     size = 2,
#     symbol = "circle",
#     line = list(
#       color = "rgb(256, 256, 256)",
#       width = 1
#     ),
#     opacity = 0.9
#   ),
#   type = "scatter3d"
# )
# 
# data <- list(list1, list2)
# 
# 
# # layout <- list(margin = list(
# #   l = 0,
# #   r = 0,
# #   b = 0,
# #   t = 0))
# 
# response <- py$plotly(data, kwargs=list(filename="random-walk", fileopt="overwrite"))
# url <- response$url

#--------------------------------------------------------------
#Select s2 = s4 = 0
# TO BE CHANGED
#--------------------------------------------------------------
# isPlanar = icdf$s2 == 0 & icdf$s4 == 0
# ttm_planar = icdf[which(isPlanar),]
# ttm_planar_max = ttm_planar[which.max(ttm_planar$dHz),]
# icdf_max = icdf[which.max(icdf$dHz),]
# icdf_inf_planar = icdf[which(icdf$dHz <= max(ttm_planar$dHz)),]
# 
# 
# pt = plotdf_point(icdf_inf_planar , "xPH", "yPH", "xPH", "yPH", "dHz", "dHz", 0)
# pt = pt +geom_point(data = ttm_planar, aes(x = xPH, y = yPH), size = linesize[["point"]],color = "red")
# pt = pt +geom_point(data = ttm_planar_max, aes(x = xPH, y = yPH), color = "green",size = 3)
# pt = pt +geom_point(data = icdf_max, aes(x = xPH, y = yPH), color = "purple",size = 3)
# pt
# 
# plotdf_point(icdf, "xPH", "zPH", "xPH", "zPH", "dHz", "dHz", 0)+geom_point(data = ttm_planar, 
#                                                                             aes(x = xPH, y = zPH), 
#                                                                             size = linesize[["point"]],
#                                                                             color = "red")  
# 
# ttm_dHz = icdf[which(abs(icdf$dHz - ttm_planar$dHz[1]) <= 1e-6),]
# isPlanar = ttm_dHz $s2 == 0 & ttm_dHz $s4 == 0
# ttm_dHz_planar = ttm_dHz[which(isPlanar),]
# plotdf_point(ttm_dHz , "xPH", "zPH", "xPH", "zPH", "dHz", "dHz", 0)
# plotdf_point(ttm_dHz , "xPH", "yPH", "xPH", "yPH", "dHz", "dHz", 0)+geom_point(data = ttm_dHz_planar, aes(x = xPH, y = yPH), size = 3,color = "red")
# plotdf_point(ttm_dHz , "yPH", "zPH", "yPH", "zPH", "dHz", "dHz", 0)



#TODO: 
#Repeat the process above BUT for a given maximum value of dHz (ie check that the planar @ dHz fixed and @ t= 0is first unique and enclose the rest of the curves)
#Uniqueness is okay: we have fixed s2 = s4 = t = 0 + dHz is fixed ==> gives unique relation s1 = f(s3). 

#In fact, dHz is an increasing value of all (s1, s2, s3, s4), when taken POSITIVE.
#But what about the relationship between the si and xPH and yPH? is it monotonous? It depends on the values of the matrix P(t=0)
#But pretty sure that the planar case enclose the other ones. Can be seen afterwards, maybe.

#For example, we can select a planar case of about 1000 km of radius, take its dHz value, 
#impose that dHz ~ this value for all (s1, s2, s3, s4) on a given grid.

#Do we need to use (s1,s2,3) and s4  = f(s1,s2,3)?
#take a look at the CRTBP case to refresh memory about the best choice for s4 (is it better with s1, s2 or s3 instead of s4?)

#Need to implement a rvf-free Poincare maps process, in order to get the z = 0 events without having to compute the whole thing.
#Quite impossible, since we want to keep a 10^-6 error!
#Then, once the map is obtained, check that, actually, all points are enclosed (or not!) inside the planar case.
#Maybe not, since it will corresponds to a weird Poincaré Maps for a non autonomous system (section is not t = 0 mod[T] but z = 0)

