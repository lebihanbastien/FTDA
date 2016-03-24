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
source("source/folder.R")
source("source/plot.R")
source("source/routines.R")
source("source/multiplot.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li    = "L2"
MODEL = "RTBP"
FWRK  = "EM"
Type  = "planar0" #selection or global
vorders = c(10, 15, 20, 25, 30, 35, 40)
ofs_order = 0
currentfolder = paste0(printfolder(MODEL, FWRK, Li))


#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK);


#------------------------------------------------
# Building the data.frame of results
#------------------------------------------------
ttm_all = data.frame()
for (i in vorders)  #loop on the orders
{
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0(currentfolder, "Serv/eOm_", Type, "_ofs_",ofs_order,"_order_", toString(i))
  filename = paste0(fileprefix, ".txt")
  
  
  # Load csv source
  #------------------------------------------------
  if (file.exists(filename))
  {
    ttm_c  = read.csv(filename, header = T, sep = ",")
  }else
  {
    ttm_c = data.frame()
  }
  
  # Distinguish between order <= 30 and < 30
  #------------------------------------------------
  ttm_c$isG30 = (i > 30)
    
  #rbind in ttm_all
  ttm_all = rbind(ttm_all, ttm_c)
}

#------------------------------------------------
# Postprocessing
#------------------------------------------------
# Compute -log10(precision)
ttm_all$log10eOm = log10(ttm_all$eOm)
#Get rid of the origin
isOrigin = ttm_all$s1 == 0 &
  ttm_all$s2 == 0 & ttm_all$s3 == 0 & ttm_all$s4 == 0
ttm_all = ttm_all[which(!isOrigin),]
# From NC to EM units
ttm_all = NCtoC(ttm_all, gamma)
# From EM to physical units
ttm_all = CtoPH(ttm_all, L)
# Radii from Li
ttm_all$rNC = sqrt(ttm_all$x^2+ttm_all$y^2+ttm_all$z^2)
ttm_all$rPH = sqrt(ttm_all$xPH^2+ttm_all$yPH^2+ttm_all$zPH^2)
#Select only positive s1
#ttm_all = ttm_all[which(ttm_all$s1 >=0),]

#------------------------------------------------
# Plot
#------------------------------------------------
#Plot the precision as a function of rPH
#------------------------------------------------
xlabel  = paste0("Distance from ", Li, " [km]");
ylabel  = expression(log[10](e[O]));
pPH = plotdf_line(ttm_all, "rPH", "log10eOm", xlabel, ylabel, "order", "Order", 1)+legend_inside_theme
pPH = pPH #+ scale_x_continuous(limits=c(0.0, 25000)) 
pPH = pPH #+ scale_y_continuous(limits=c(-9, 0.5))
pPH

#------------------------------------------------
#Plot the precision as a function of s1, with s(0) = (s1, 0, 0, 0) in planar0 case
#------------------------------------------------
xlabel  = expression(s[1]);
ylabel  = expression(log[10](e[O]));
pS1 = plotdf_line(ttm_all, "s1", "log10eOm", xlabel, ylabel, "order", "Order", 1)
pS1 = pS1 #+ scale_y_continuous(limits=c(-7.5, -5))

#------------------------------------------------
#Plot the precision as a function of s1, with s(0) = (s1, 0, 0, 0) in planar0 case
#------------------------------------------------
xlabel  = expression(s[2]);
ylabel  = expression(log[10](e[O]));
pS2 = plotdf_line(ttm_all, "s2", "log10eOm", xlabel, ylabel, "order", "Order", 1)
pS2 = pS2 #+ scale_y_continuous(limits=c(-7.5, -5))

#------------------------------------------------
#Plot the precision as a function of the energy gap
#------------------------------------------------
pE = plotdf_line(ttm    = ttm_all, 
                 colx   = "dHz", 
                 coly   = "log10eOm", 
                 xlabel = expression("Energy gap from the origin"),
                 ylabel = expression(log[10](e[O])), 
                 colorCol = "order", 
                 colorLabel = "Order", 
                 isColorFac = 1)#,
                 #lineTypeCol = "isG30",
                 #lineTypeLabel = "Order > 30")
# pE = pE+scale_x_continuous(limits=c(0.0, 0.02))
pE = pE+scale_y_continuous(limits=c(-12, 0))
pE


#------------------------------------------------pE
#Plot the precision as a function of yEM
#------------------------------------------------
xlabel  = expression("yEM");
ylabel  = expression(log[10](e[O]));
pY = plotdf_line(ttm_all, "yEM", "log10eOm", xlabel, ylabel, "order", "Order", 1)

#------------------------------------------------
#Plot the energy gap as a function of s1
#------------------------------------------------
plotdf_line(ttm_all, "s1", "dHz", expression(s[1]), "Energy gap from the origin", "order", "Order", 1)
plotdf_line(ttm_all, "s2", "dHz", expression(s[2]), "Energy gap from the origin", "order", "Order", 1)

#------------------------------------------------
#dHz as a function of rPH
#------------------------------------------------
plotdf_line(ttm_all, "rPH", "dHz", "rPH", "dHz", "order", "Order", 1)

#------------------------------------------------
#Save in eps file
#------------------------------------------------
#ggsave(pS1, file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_RCM.eps"))
#ggsave(pY, file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_EM.eps"))
#ggsave(pPH, file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_NC.eps"))
ggsave(pE, file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order,  "_", Type, "_Energy_EML2_RTBP.pdf"))
ggsave(pE, file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order,  "_", Type, "_Energy_EML2_RTBP.eps"))