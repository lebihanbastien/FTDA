############################################################
# R script to handle an invariance map of the QBCP around EML1
############################################################
#Saving 9.01 x 7.69 in image
#------------------------------------------------
# Select the plots
#------------------------------------------------
Li        = "L1"
MODEL     = "QBCP"
FWRK      = "EM"
si        = "s1"
sj        = "s3"
Type      = paste0(si, sj) #"s1fs2s3" #
dHz       = 0.06 #0.005 #0.05 or 0.06
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
vorders   = c(10, 15, 20, 25, 30)
source("Imap/Imap_2D.R")

#Select very high values of too high values (artifacts)
#ttm_bad = ttm_all[which(ttm_all$log10eOm > 0),] 

#------------------------------------------------
# Plot difference between orders
#------------------------------------------------
vorders   = c(10, 15, 20, 25, 30)
pph_withAxes = TRUE
pph_withNotes= FALSE
pph_withMoon = FALSE
pph_limits   = TRUE
pph_x_lim    = c(-0.9, -0.80)
pph_y_lim    = c(-0.15, 0.15)
source("Imap/Imap_2D_Build.R")
ordersToSub = c(25, 30)
source("Imap/Imap_2D_SUB.R")
pph
