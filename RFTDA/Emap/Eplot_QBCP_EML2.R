############################################################
# R script to handle a precision map of the QBCP around EML2
############################################################

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
Lit        = "$L_2$"
MODEL     = "QBCP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
Type      = "s1s2s3s4" #"s1s3"#"s1s2s3s4" #"s1eqms3" #s1, s1eqs3, s1eqms3, s3, s1eq2s2eqs3eq2s4
order     = "20"
ofs_order = "30"
Energy    = "0"#"-10" #0.01 #0.01 #0.0025, 0.005, 0.0075
vorders   = seq(6,26,4)

#------------------------------------------------
# Plotting parameters
#------------------------------------------------
pE_limits=TRUE
pE_limits_x=c(0.0, 0.02)
pE_limits_y=c(-12, -1)

pS1_limits=FALSE
pS1_limits_x=c(0, 15)
pS1_limits_y=c(-9, -3)

pX_limits=FALSE
pX_breaks_x =FALSE
pX_breaks_y =FALSE
pX_limits_x        = c(NaN, NaN)
pX_limits_y        = c(-12, -1)
pX_breaks_x_values = seq(-1.19,-1.10,0.01)
pX_breaks_y_values = seq(-10,0,2)

pY_limits=FALSE
pY_breaks_x =FALSE
pY_breaks_y =FALSE
pY_limits_x=c(-0.07, 0.07)
pY_limits_y=c(-11, -2)
pY_breaks_x_values = seq(-0.10,0.10,0.01)
pY_breaks_y_values = seq(-10,0,2)

pG_limits_y = c(-12, -1);
pG_breaks_y = seq(-16,0,2);

#------------------------------------------------
# Process data
#------------------------------------------------
#source("Emap/Eplot.R")
source("Emap/Eplot_TEX.R")