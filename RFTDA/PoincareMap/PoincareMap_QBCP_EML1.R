############################################################
# R script to handle a poincare map of the QBCP around EML1,2
############################################################

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L1"
MODEL     = "QBCP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
Energy    = "0.01"
nst       = "001"
add       = "t0_projtenth_" #t0_hardprec_ for 0.005 and 0.0075
order     = "20" #mainly 20
ofs_order = "30"
fileext    = ".bin"
fileprefix = "Serv_pm_Energy_0.01_order_20_ofs_30_SINGLE_INT"
fileprefixnodots  = "Serv_pm_Energy_001_order_20_ofs_30_SINGLE_INT"

#------------------------------------------------
# Limits for plots
#------------------------------------------------
pEM_limits   = FALSE
pEM_limits_x = c(-0.03, 0.04)
pEM_limits_y = c(-0.1, 0.1)



#------------------------------------------------
# Process data
#------------------------------------------------
source("PoincareMap/PoincareMap.R")