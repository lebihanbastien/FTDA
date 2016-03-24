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
MODEL = "QBCP"
FWRK  = "EM"
#FWRK  = "SEM"
#Type  = "3d_complete" #selection or global
Type  = "3d_s1s2" #selection or global
vorders = c(10)
ofs_order = 30
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

#------------------------------------------------
# Building the data.frame of results
#------------------------------------------------
ttm_all = data.frame()
for (i in vorders)  #loop on the orders
{
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0(currentfolder, "Serv/eOm_Serv_", Type, "_ofs_",ofs_order,"_order_", toString(i))
  filename = paste0(fileprefix, ".txt")
  
  
  # Load csv source
  #------------------------------------------------
  if (file.exists(filename))
  {
    ttm_c1  = read.csv(filename, header = T, sep = ",")
  }else
  {
    ttm_c1 = data.frame()
  }
  #Style
  ttm_c1$Style = "GS"
  #Order
  ttm_c1 = ttm_c1[order(ttm_c1$s1, ttm_c1$s2, ttm_c1$s3),]
  # Compute -log10(precision)
  ttm_c1$log10eOm = log10(ttm_c1$eOm)
  #rbind in ttm_all
  ttm_all = rbind(ttm_all, ttm_c1)
  
  
  
  
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0(currentfolder, "Serv/eOm_Serv_NF_", Type, "_ofs_",ofs_order,"_order_", toString(i))
  filename = paste0(fileprefix, ".txt")
  
  
  # Load csv source
  #------------------------------------------------
  if (file.exists(filename))
  {
    ttm_c2  = read.csv(filename, header = T, sep = ",")
  }else
  {
    ttm_c2 = data.frame()
  }
  #Style
  ttm_c2$Style = "NF"
  #Order
  ttm_c2 = ttm_c2[order(ttm_c2$s1, ttm_c2$s2, ttm_c2$s3),]
  # Compute -log10(precision)
  ttm_c2$log10eOm = log10(ttm_c2$eOm)
  
  #rbind in ttm_all
  ttm_all = rbind(ttm_all, ttm_c2)
}


#------------------------------------------------
#Substracting
#------------------------------------------------
ttm_sub = data.frame(s1=ttm_c1$s1, s2=ttm_c1$s2, s3=ttm_c1$s3, log10eOm=ttm_c1$log10eOm-ttm_c2$log10eOm)




#------------------------------------------------
# Plot
#------------------------------------------------
#------------------------------------------------
# Complet Plot (eOm = f(s1, s2))
#------------------------------------------------
Ri     = 15000;
dRi    = 400;
Rcolor = muted("green")
Rsize  = 0.2
dHzi   = 0.0025
ddHzi  = 1e-4

#Multiplot for all orders
#---------------------
ttm_l   = list();
ttm_lr  = list();
ttm_le  = list();
ttm_str = list();
ppl     = list();
for(i in 1:length(vorders))
{
  order = vorders[i];
  ttm_l[[i]]  = ttm_sub[which(ttm_all$order == order),]      #select order
  ppl[[i]]    = plotdf_tile(ttm_l[[i]], "s1", "s2", isLegendOn, 0)
}
#Actual plot in multiplot format
pMult = multiplot(plotlist = ppl, cols = 2)

#Legend only (don't forget to display the legend in plotfd_tile!)
#--------------------
#legend <- g_legend(ppl[[1]]) 
#pLegend = grid.draw(legend) 

#Save in eps file
#------------------------------------------------
ggsave(ppl[[1]], file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", "Comparison", ".eps"))

#TODO: METTRE AU POINT LES PLOTS DANS X-Y AVEG GGPOINT

