############################################################
# R script to handle a precision map of the CRTBP
# Has to be used files such as e.g. PoincareMap_QBCP_EML2.R
#
# 01/2016
############################################################

#------------------------------------------------
# Init
#------------------------------------------------
source("source/init.R")

#------------------------------------------------
#Current working folder
#------------------------------------------------
currentfolder = paste0(printfolder(MODEL, FWRK, Li))

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK)