############################################################
# R script to handle a poincare map of the QBCP around EML1,2
############################################################

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
MODEL     = "QBCP"
FWRK      = "EM"
METHOD    = "_DUAL_INT"#_SINGLE_INT"
Energy    = "0.005" #available for 30: 0.012, 0.011, 0.01
nst       = "0005"
add       = "t0_hardprec_"
order     = "20" #order 20 for 0.005 and 0.0025
ofs_order = "30"
fileext    = ".txt"

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