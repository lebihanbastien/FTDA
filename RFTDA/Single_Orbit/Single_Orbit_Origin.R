# Script to plot an native_orbit along with its associated precisions
# in the parameterization method of the QBFBP/RTBP around L1/L2 of the Earth-Moon system
#----------------------------------------------------------------------------------------

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
FWRK2  = "SEM"
currentfolder = paste0(plotfolder(MODEL, FWRK, Li), "orbits/")
maxPrec = 1e-6
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

# Orders
maxOrder = 0

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
orbit = NCtoC(orbit, gamma, c1)
# From EM to physical units
orbit = CtoPH(orbit, L)

#------------------------------------------------
#Origin
#------------------------------------------------
Ldf = data.frame(x = 0, y = 0, z =0)
# From NC to EM units
Ldf = NCtoC(Ldf, gamma, c1)
# From EM to physical units
Ldf = CtoPH(Ldf, L)


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
porb = fplot_path (orbit,  "xEM", "yEM",  "x [-]",  "y [-]", "type", "", TRUE, "type")
# Starting point
porb = porb + geom_point(data = native_orbit_start, aes(xEM, yEM), size = 4, colour = "black", fill = "white", pch = 21)
# Origin point
porb = porb + geom_point(data = Ldf, aes(xEM, yEM), size = 4, colour = "black", fill = "white", pch = 21)
#Theme
porb = porb+custom_theme+theme(legend.position="bottom")

#Add an arrow to give the direction of motion along the orbit
ai= 100
porb = porb + geom_segment(aes(x = orbit$xEM[ai], y = orbit$yEM[ai], xend = orbit$xEM[ai+1], yend = orbit$yEM[ai+1]), 
                           colour = muted("white"), 
                           arrow = arrow(length = unit(0.4, "cm"), type = "closed"))

#Add annotation
if(Li == "L2")
{
  porb = porb + annotate("text", x = Ldf$xEM+1e-07, y = Ldf$yEM+3e-07, label = "L[2]", size = 10, parse =TRUE) 
}else
{
  porb = porb + annotate("text", x = Ldf$xEM+1e-07, y = Ldf$yEM+3e-07, label = "L[1]", size = 10, parse =TRUE)
}

#Scaling, if necessary
#porb = porb + scale_x_continuous(limits = c(-1e-5, 1e-5))  #cont scale on the x axis 
#porb = porb + scale_y_continuous(limits = c(-2e-5, 2e-5))  #cont scale on the y axis
porb = porb+scale_color_manual(guide = FALSE, values=c("black"))

#Save in eps file
#------------------------------------------------
ggsave(porb, file = paste0(currentfolder, "orbit_", Li, ".eps"))



