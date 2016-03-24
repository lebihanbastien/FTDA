# Script to plot the precision of the potential and the vector field 
# in the parameterization method of the QBFBP/RTBP around L1/L2 of the Earth-Moon system
#----------------------------------------------------------------------------------------
#Set working directory - uncomment if necessary
#setwd("~/BackUpBox/PhD/OOFTDA/RFTDA/")

# Load libraries
#------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)

# Load Source files
#------------------------------------------------
source("source_folder.R")
source("source_plot.R")
source("source_routines.R")

# Select Models & libration point
#------------------------------------------------
Li = "L2"
MODEL = "RTBP"
currentfolder = printfolder(MODEL, Li)

#Normalized units (gamma, c1)
#------------------------------------------------
gamma = gamma(Li);
c1 =  c1(Li);
L = 384400; #Earth-Moon distance in [km]

# Type of plot (can be changed at will)
#------------------------------------------------
if(MODEL == "QBCP") 
{
  fplot = plotdf_line;
  
}else
{
  fplot = plotdf_line_dashed;
}
# Orders
if(MODEL == "QBCP") dfindex = c(5,10,15,20) else dfindex = c(2,5,10,20)

#--------------------------------------------------------------------------------------------------------------------------
#                                             Potential
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "PotentialvsW_"
for(i in 1:length(dfindex))
{
  dflist[[i]] = read.table(paste0(currentfolder,dftype,dfindex[i],".txt"), header = F)
  colnames(dflist[[i]]) = c("x", "y")
  dflist[[i]]$order = dfindex[i]
}
# Concatenate the results
PvsW = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfindex)){PvsW = rbind(PvsW, dflist[[i]])}

#Euclidian distance in km
PvsW$xPH = gamma*L*PvsW$x

#Plot
#------------------------------------------------
#Actual plot
p = fplot(PvsW,  "xPH", "y", 
                 "\n |W(s,0)| [km]",  "Error on the potential [-]\n",
                 "order", "Order", 1)
#Cont scale on the x axis 
p = p + scale_x_continuous(limits=c(0.0, 19000))
#Logscale on the y axis
if(MODEL == "QBCP") p = p + scale_y_log10(limits=c(1e-12, 1e-3), breaks = 10^(-8:0*2))else 
p = p + scale_y_log10(limits=c(1e-15, 1e-1), breaks = 10^(-8:0*2)) 
#Display the plot
#p = p+theme(legend.justification = c(0,0), legend.position  = c(0.02,0.6) )
p = p + legend_inside_theme



#--------------------------------------------------------------------------------------------------------------------------
#                                             Vector field
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "VFvsW_"
for(i in 1:length(dfindex))
{
  dflist[[i]] = read.table(paste0(currentfolder,dftype,dfindex[i],".txt"), header = F)
  colnames(dflist[[i]]) = c("x", "y")
  dflist[[i]]$order = dfindex[i]
}
# Concatenate the results
VFvsW = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfindex)){VFvsW = rbind(VFvsW, dflist[[i]])}

#Euclidian distance in km
VFvsW$xPH = gamma*L*VFvsW$x

#Plot
#------------------------------------------------
pvf = fplot(VFvsW,  "xPH", "y", 
                    "\n |W(s,0)| [km]",  "Error on the vector field [-]\n",
                    "order", "Order", 1)
#Scaling
pvf = pvf + scale_y_log10(limits=c(1e-15, 1e-1), breaks = 10^(-8:0*2))             #log  scale on the y axis     
pvf = pvf + scale_x_continuous(limits=c(0.0, 19000))
#Display the plot
pvf = pvf+legend_inside_theme

#Save in eps file
#------------------------------------------------
ggsave(p,   file=paste0(currentfolder, "R_PotentialvsW.eps"))
ggsave(pvf, file=paste0(currentfolder, "R_VFvsW.eps"))
