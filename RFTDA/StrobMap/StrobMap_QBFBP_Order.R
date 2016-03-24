# R script to handle a poincare map of the QBFBP around EML1,2
#---------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
# Init
#------------------------------------------------------------------------------------
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
source("source/source_folder.R")
source("source/source_plot.R")
source("source/source_routines.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
MODEL     = "QBCP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
Energy    = "0.0025"
order     = "15" #available for L2: 20 and 15
vorders   = c(10,15,20,30)
ofs_order = "30"

#Current working folder
currentfolder = paste0(printfolder(MODEL, FWRK, Li))

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
gamma = gamma(Li, FWRK);
c1 =  c1(Li, FWRK);
if(FWRK == "SEM")
{
  L = 149.60e6; #Sun-Earth distance in [km]
}else{
  L = 384400;   #Earth-Moon distance in [km]
}

#------------------------------------------------------------------------------------
# Data reading
#------------------------------------------------------------------------------------
#------------------------------------------------
# Filename to check
#------------------------------------------------
fileprefix = paste0(currentfolder, "Serv/Serv_tm_Energy_", Energy, "_order_", order, "_ofs_", ofs_order, METHOD)
filename = paste0(fileprefix, ".txt")

#------------------------------------------------
# Load csv source
#------------------------------------------------
# if (file.exists(filename))
# {
#   tmdf  = read.csv(filename, header = T, sep = ",")
# }else
# {
#   tmdf = data.frame()
# }


#------------------------------------------------
# Building the data.frame of results
#------------------------------------------------
tmdf = data.frame()
for (i in vorders)  #loop on the orders
{
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0(currentfolder, "Serv/Serv_tm_Energy_", Energy, "_order_", toString(i), "_ofs_", ofs_order, METHOD)
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
  
  #rbind in ttm_all
  tmdf = rbind(tmdf, ttm_c)
}


#------------------------------------------------------------------------------------
# Postprocessing
#------------------------------------------------------------------------------------

#--------------------
# new columns
#--------------------
# From NC to EM units
tmdf = NCtoC(tmdf, gamma)
# From EM to physical units
tmdf = CtoPH(tmdf, L)
# Radii from Li
tmdf$rNC = sqrt(tmdf$x^2+tmdf$y^2+tmdf$z^2)
tmdf$rPH = sqrt(tmdf$xPH^2+tmdf$yPH^2+tmdf$zPH^2)
#Parity of event
tmdf$parity = tmdf$number%%2
#log(ePm)
tmdf$log10ePm = log10(tmdf$ePm)

#--------------------
# Selecting particular values
#--------------------
#Starting points
tmdf0 = tmdf[which(tmdf$number == 0),]

#------------------------------------------------
# Plots
#------------------------------------------------
#--------------------------
# Computation standard deviation & mean
#--------------------------
tmdfe_sum = ddply(tmdf, .(label, order), summarize, dHw.sdpc = sd(dHw)/mean(dHw), dHw.sd = sd(dHw), dHw.mean = mean(dHw), maxT = max(t), maxN = max(number))

pSd = plotdf_point(tmdfe_sum , "dHw.mean", "dHw.sd", expression(paste("\n ", mu(dH))), expression(paste(sigma(dH), "\n")), "order", "Order", 1)
pSd
ggsave(pSd, file=paste0(currentfolder, "PlanarStrobMap_sd_vs_mean.eps")) #Save

pSd = plotdf_point(tmdfe_sum , "dHw.mean", "dHw.sdpc", 
                   expression(paste("\n ", mu(dH))), 
                   expression(paste(frac(sigma(dH), mu(dH)), "\n")), 
                   "order", "Order", 1)
pSd = pSd + scale_y_continuous(limits = c(0, 0.35))
pSd = pSd + theme(axis.title.y = element_text(angle = 0))
pSd
ggsave(pSd, file=paste0(currentfolder, "PlanarStrobMap_sdmean_vs_mean.eps")) #Save
