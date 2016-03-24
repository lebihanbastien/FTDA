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
source("source/source_folder.R")
source("source/source_plot.R")
source("source/source_routines.R")
source("source/multiplot.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li    = "L2"
MODEL = "QBCP"
FWRK  = "EM"
Type  = "s1s2s3s4" #selection or global
vorders = c(10, 15, 20, 25, 30)
ofs_order = 30
dHz = 0.05
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
isLegendOn = 0;
legendOnly = 0;

#------------------------------------------------
# Building the data.frame of results
#------------------------------------------------
#Function to read the data
read_imap_4d <- function(vorders)
{
  #New data frame
  ttm_all = data.frame()
  for (i in vorders)  #loop on the orders
  {
    # Filename to check
    fileprefix = paste0(currentfolder, "Serv/eIm_", Type, "_ofs_",ofs_order,"_order_", toString(i), "_hmax_", toString(dHz))
    filename = paste0(fileprefix, ".txt")
    
    # Load csv source
    if (file.exists(filename))
    {
      ttm_c  = read.csv(filename, header = T, sep = ",")
    }else
    {
      ttm_c = data.frame()
    }
    
    #rbind in ttm_all
    ttm_all = rbind(ttm_all, ttm_c)
  }
  
  #Return the result
  return(ttm_all)
  
}

ttm_all = read_imap_4d(vorders) #do it once !


#------------------------------------------------
# Postprocessing
#------------------------------------------------
# Compute -log10(precision)
ttm_all$log10eOm = log10(ttm_all$eOm)
#Get rid of the origin
isOrigin = ttm_all$s1 == 0 & ttm_all$s2 == 0 & ttm_all$s3 == 0 & ttm_all$s4 == 0
ttm_all = ttm_all[which(!isOrigin),]
# From NC to EM units
ttm_all = NCtoC(ttm_all, gamma)
# From EM to physical units
ttm_all = CtoPH(ttm_all, L)
# Radii from Li
ttm_all$rNC = sqrt(ttm_all$x^2+ttm_all$y^2+ttm_all$z^2)
ttm_all$rPH = sqrt(ttm_all$xPH^2+ttm_all$yPH^2+ttm_all$zPH^2)

#------------------------------------------------
# Function to perform cuts
#------------------------------------------------
imap_4d_cut <- function(si, sj, ind)
{
  if(si == "s1" && sj == "s3")
  {
    
    #Select all values of s2 == 1/2 svec[ind]
    ttm_cut = ttm_all[which(abs(ttm_all$s2-0.5*svec[ind]) < 1e-5),]
    #Select all values of s4 == 1/2 svec[ind]
    ttm_cut = ttm_cut[which(abs(ttm_cut$s4-0.5*svec[ind]) < 1e-5),]
    
  }else if(si == "s1" && sj == "s2")
  {
    #Select all values of s3 == svec[ind]
    ttm_cut = ttm_all[which(abs(ttm_all$s3-svec[ind]) < 1e-5),]
    #Select all values of s4 == 1/2 svec[ind]
    ttm_cut = ttm_cut[which(abs(ttm_cut$s4-0.5*svec[ind]) < 1e-5),]
    
  }else if(si == "s1" && sj == "s4")
  {
    #Select all values of s2 == 1/2 svec[ind]
    ttm_cut = ttm_all[which(abs(ttm_all$s2-0.5*svec[ind]) < 1e-5),]
    #Select all values of s3 == svec[ind]
    ttm_cut = ttm_cut[which(abs(ttm_cut$s3-svec[ind]) < 1e-5),]
    
  }else if(si == "s2" && sj == "s3")
  {
    #Select all values of s1 == svec[ind]
    ttm_cut = ttm_all[which(abs(ttm_all$s1-svec[ind]) < 1e-5),]
    #Select all values of s4 == 1/2 svec[ind]
    ttm_cut = ttm_cut[which(abs(ttm_cut$s4-0.5*svec[ind]) < 1e-5),]
    
  }else if(si == "s2" && sj == "s4")
  {
    #Select all values of s1 == svec[ind]
    ttm_cut = ttm_all[which(abs(ttm_all$s1-svec[ind]) < 1e-5),]
    #Select all values of s3 == svec[ind]
    ttm_cut = ttm_cut[which(abs(ttm_cut$s3-svec[ind]) < 1e-5),]
    
  }else if(si == "s3" && sj == "s4")
  {
    #Select all values of s1 == svec[ind]
    ttm_cut = ttm_all[which(abs(ttm_all$s1-svec[ind]) < 1e-5),]
    #Select all values of s2 == 1/2 svec[ind]
    ttm_cut = ttm_cut[which(abs(ttm_cut$s2-0.5*svec[ind]) < 1e-5),]
  }
  
  #Return the data frame
  return(ttm_cut)
}