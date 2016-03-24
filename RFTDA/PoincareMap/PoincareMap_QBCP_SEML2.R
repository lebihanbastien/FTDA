############################################################
# R script to handle a poincare map of the QBCP around SEML2
############################################################

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
MODEL     = "QBCP"
FWRK      = "SEM"
METHOD    = "_SINGLE_INT"
Energy    = "1e-05"
nst       = "1e-05"
add       = "t0_"
order     = "20"
ofs_order = "30"

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