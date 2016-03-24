############################################################
# R script to handle a precision map of the RTBP around EML2
############################################################

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
MODEL     = "RTBP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
Type      = "s1" #s1s2s3, planar0
order     = "20"
ofs_order = "0"
vorders   = c(10, 15, 20, 25, 30, 35, 40) #c(10, 15, 20, 25, 30) for planar0, c(10, 20, 30, 40) for s1s2s3
Energy    = 0

#------------------------------------------------
# Plotting parameters
#------------------------------------------------
pE_limits=FALSE
pS1_limits=FALSE
pX_limits=FALSE
pX_breaks_x=FALSE
pX_breaks_y=FALSE
pY_limits=FALSE
pY_breaks_x=FALSE
pY_breaks_y=FALSE

#------------------------------------------------
# Process data
#------------------------------------------------
source("Emap/Eplot_RTBP.R")

