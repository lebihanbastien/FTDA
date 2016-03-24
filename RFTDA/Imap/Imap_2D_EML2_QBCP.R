############################################################
# R script to handle an invariance map of the QBCP around EML1
############################################################

#------------------------------------------------
# Select the plots
#------------------------------------------------
Li        = "L2"
MODEL     = "QBCP"
FWRK      = "EM"
si        = "s1"
sj        = "s3"
sis       = expression(s[1])
sjs       = expression(s[3])
Type      = "s1s3" #paste0(si, sj)
dHz       = 0.01 #0.0075 #0.05
ofs_order = 30


#------------------------------------------------
#Additionnal parameters
#------------------------------------------------
isLegendOn  = 0
legendOnly  = 0
legendOnTop = 0

#------------------------------------------------
# Plot the various orders
#------------------------------------------------
vorders   = c(10, 16, 20, 26, 30)
source("Imap/Imap_2D.R")

#------------------------------------------------
# Plot difference between orders
#------------------------------------------------
vorders   = c(10, 16, 20, 26, 30)
pph_withAxes = TRUE
pph_withNotes= FALSE
pph_withMoon = FALSE
pph_limits   = FALSE
pph_x_lim    = c(-0.9, -0.80)
pph_y_lim    = c(-0.15, 0.15)
source("Imap/Imap_2D_Build.R")
ordersToSub = c(16, 20)
source("Imap/Imap_2D_SUB.R")
pph
