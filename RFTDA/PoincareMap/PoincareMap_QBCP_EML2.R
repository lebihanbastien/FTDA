############################################################
# R script to handle a poincare map of the QBCP around EML1,2
############################################################

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
MODEL     = "QBCP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
Energy    = "0.012" #available for 30: 0.012, 0.011, 0.01
nst       = "0012"
add       = "t0_" #"t0_projtenth_"
order     = "30" #order 20 for 0.005 and 0.0025
ofs_order = "30"
fileext   = ".txt"
fileprefix = "Serv_pm_t0_Energy_0.012_order_30_ofs_30_SINGLE_INT"
fileprefixnodots  = "Serv_pm_t0_Energy_0012_order_30_ofs_30_SINGLE_INT"

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