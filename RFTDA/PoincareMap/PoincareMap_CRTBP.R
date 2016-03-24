############################################################
# R script to handle a poincare map of the QBCP around SEML2
############################################################

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
MODEL     = "RTBP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
dH0       = 0.015
nst       = "0015"
Energy    = toString(dH0)
order     = "30"
add       = ""
ofs_order = "0"

#------------------------------------------------
# Limits for plots
#------------------------------------------------
pEM_limits   = TRUE
pEM_limits_x = c(-0.03, 0.04)
pEM_limits_y = c(-0.1, 0.1)

#------------------------------------------------
# Process data
#------------------------------------------------
source("PoincareMap/PoincareMap.R")