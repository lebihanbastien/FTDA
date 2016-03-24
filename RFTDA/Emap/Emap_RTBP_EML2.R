############################################################
# R script to handle a poincare map of the CRTBP around EML1,2
############################################################

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
MODEL     = "RTBP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
Type      = "s1s2s3" #selection or global
si        = "s1"
sj        = "s3"
ofs_order = "0"
vorders   = c(40)

#------------------------------------------------
#Additionnal parameters
#------------------------------------------------
isLegendOn = 1;
legendOnly = 0;

#------------------------------------------------
# Process data
#------------------------------------------------
source("Emap/Emap.R")