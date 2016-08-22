# R script to handle a poincare map of the QBFBP around EML1,2
#---------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Init
#------------------------------------------------------------------------------------
# R options
#------------------------------------------------
options(digits = 15)

#------------------------------------------------
# Load libraries
#------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(tikzDevice)
library(latex2exp)


#------------------------------------------------
# Load Source files
#------------------------------------------------
source("source/folder.R")
source("source/plot.R")
source("source/routines.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L1"
MODEL     = "QBCP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
Energy    = "0.0025"
order     = "20" #available for EML2: 10,15,20,30
ofs_order = "30"

#Current working folder
currentfolder = paste0(printfolder(MODEL, FWRK, Li))

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
gamma = gamma(Li, FWRK);
c1 =  c1(Li, FWRK);
if(FWRK == "SEM")
{
  L = 149.60e6; #Sun-Earth distance in [km]
}else{
  L = 384400;   #Earth-Moon distance in [km]
}

#------------------------------------------------------------------------------------
# Data reading
#------------------------------------------------------------------------------------
#------------------------------------------------
# Filename to check
#------------------------------------------------
fileprefix = paste0(currentfolder, "Serv/Serv_tm_Energy_", Energy, "_order_", order, "_ofs_", ofs_order, METHOD)
filename = paste0(fileprefix, ".txt")

#------------------------------------------------
# Load csv source
#------------------------------------------------
if (file.exists(filename))
{
  tmdf  = read.csv(filename, header = T, sep = ",")
}else
{
  tmdf = data.frame()
}

#------------------------------------------------------------------------------------
# Postprocessing
#------------------------------------------------------------------------------------

#--------------------
# new columns
#--------------------
# From NC to EM units
tmdf = NCtoSYS(tmdf, gamma, c1)
# From NC to C units
tmdf = NCtoC(tmdf, gamma)
# From C to physical units
tmdf = CtoPH(tmdf, L)
# From SYS to physical units
tmdf = SYStoPH(tmdf, L)
# Radii from Li
tmdf$rNC = sqrt(tmdf$x^2+tmdf$y^2+tmdf$z^2)
tmdf$rC = sqrt(tmdf$xC^2+tmdf$yC^2+tmdf$zC^2)
tmdf$rPH = sqrt(tmdf$xPH^2+tmdf$yPH^2+tmdf$zPH^2)
#Parity of event
tmdf$parity = tmdf$number%%2
#log(ePm)
tmdf$log10ePm = log10(tmdf$ePm)

#--------------------
# Selecting particular values
#--------------------
#Starting points
tmdf0 = tmdf[which(tmdf$number == 0),]

#Only some labels
#condition = (tmdf$label%%2 == 0)
tmdf_lab = tmdf#[which(condition),] 

#Only a certain range of precision
precFlag = tmdf$log10ePm < -4 & tmdf$number > 0
tmdf_prec = tmdf[which(precFlag),]

#Select a given layer of energy
layer.value = 0.001
layer.width = 5e-4
tmdf_layer = tmdf[which(abs(tmdf$dHw - layer.value) < layer.width),]

#Only a certain range of energy
tmdf_lab =  tmdf_lab[which(tmdf_lab$dHw < 0.02),]

#Only certain precision
precFlag = tmdf_lab$log10ePm < -6 & tmdf_lab$number > 0
tmdf_lab$precFlag = precFlag
tmdf_lab_prec = tmdf_lab[which(precFlag),]

#------------------------------------------------
#Only the labels for which all points verifies precFlag
#------------------------------------------------
nonPrecFlag = tmdf_lab$log10ePm > Inf
badLabels  = tmdf_lab$label[which(nonPrecFlag)]
tmdf_lab_prec_strict = tmdf_lab[which(!(tmdf_lab$label %in% badLabels)),]

#Corresponding solution
pLab = plotdf_point(tmdf_lab_prec_strict, "xEM", "yEM", "X [-]", "Y [-]", pointSize = 1)
pLab = pLab+big_font_theme
pLab
ggsave(pLab, width = 18, height = 18, file=paste0(currentfolder, "PlanarStrobMap.eps")) #Save

#--------------------------
# Computation standard deviation & mean
#--------------------------
tmdfe_sum = ddply(tmdf_lab_prec_strict, .(label, order), summarize, rPH = mean(rPH), rC = mean(rC), rNC = mean(rNC), dHw.sdpc = sd(dHw)/mean(dHw), dHw.sd = sd(dHw), dHw.mean = mean(dHw), maxT = max(t), maxN = max(number))
tmdfe_sum = tmdfe_sum[which(tmdfe_sum$dHw.mean < 0.02),]

#Mean energy
pSd = plotdf_point(tmdfe_sum , "dHw.mean", "dHw.sd", TeX('$\\mu(\\delta H)$'), TeX('$\\sigma(\\delta H)$'), pointSize = 3)
pSd = pSd + scale_y_continuous(limits = c(0, 8e-4))
pSd = pSd + scale_x_continuous(limits = c(0, 0.011), breaks = seq(0, 0.125, 0.0025))
pSd
ggsave(pSd, width = xSize, height = ySize, file=paste0(currentfolder, "PlanarStrobMap_sd_vs_mean_order_", order, ".eps")) #Save

#Radius
pPH = plotdf_point(tmdfe_sum , "dHw.mean", "rNC", TeX('$\\mu(\\delta H)$'), TeX('$\\sigma(\\delta H)$'), pointSize = 3)
pPH = pPH + scale_x_continuous(limits = c(0, 0.011), breaks = seq(0, 0.125, 0.0025))
pPH

#Relative mean energy
pSdr = plotdf_point(tmdfe_sum , "dHw.mean", "dHw.sdpc", TeX('$\\mu(\\delta H)$'), TeX('$\\sigma(\\delta H)$'), pointSize = 3)
pSdr = pSdr + scale_x_continuous(limits = c(0, 0.011), breaks = seq(0, 0.125, 0.0025))
pSdr


stop();



#--------------------------
# Variations of the energy
#--------------------------
# Adding some new columns:
# dHw.sdpc = sd(dHw)/mean(dHw)
# dHw.sd = sd(dHw)
# dHw.mean = mean(dHw)
# maxT = max(t)
# maxN = max(number)
#--------------------------
tmdfe = ddply(tmdf, .(label), mutate, dHw.sdpc = sd(dHw)/mean(dHw), dHw.sd = sd(dHw), dHw.mean = mean(dHw), maxT = max(t), maxN = max(number))
tmdfe = tmdfe[which(tmdfe$dHw.mean > 0.005),]

#--------------------------
#Select the label for which the energy is maximum
#--------------------------
dHw.max.label = tmdfe$label[which(tmdfe$dHw == max(tmdfe$dHw))]
tmdf_en = tmdfe[which(tmdfe$label == dHw.max.label[1]),]
#Plot the corresponding energy variations
pEn = plotdf_path(tmdf_en, "t", "dHw", "\n time [adim]", "dH(t) [adim] \n", "label", "label", 0)
pEn = pEn + scale_color_continuous(guide = FALSE)
pEn


#--------------------------
#Select the label for which the sd is minimum
#--------------------------
dHw.sd.max.label = tmdfe$label[which(tmdfe$dHw.sd == min(tmdfe$dHw.sd))]
tmdf_en = tmdfe[which(tmdfe$label == dHw.sd.max.label[1]),]
#Plot the corresponding energy variations
pEn = plotdf_path(tmdf_en, "t", "dHw", "\n time [adim]", "dH(t) [adim] \n", "label", "label", 0)
pEn = pEn + scale_color_continuous(guide = FALSE)
pEn

#Corresponding solution
plotdf_point(tmdf_en, "xEM", "yEM", "xEM", "yEM", "label", "label", 1, pointSize = 1)


#--------------------------
#Select a given range of mu
#--------------------------
tmdf_en = tmdfe[which(abs(tmdfe$dHw.mean - 0.00625) < 0.00125),]
#Corresponding solution
plotdf_point(tmdf_en, "xEM", "yEM", "xEM", "yEM", "label", "label", 1, pointSize = 1)


stop();

#------------------------------------------------
# Plots
#------------------------------------------------
#--------------------------
#Lab: vs ePm
#--------------------------
pLab = plotdf_point(tmdf_lab[which(tmdf_lab$dHw < 0.018),], "xEM", "yEM", "x [-]", "y [-]", "precFlag", expression(e[P](t) < 10^{-6}),1, pointSize = 1)
#pLab = plotdf_point(tmdf_lab_prec, "xEM", "yEM", "xEM", "yEM", "log10ePm", expression(e[P](t)), 0, pointSize = 1)
pLab = pLab+big_font_theme
pLab
ggsave(pLab, width = 18, height = 18, file=paste0(currentfolder, "PlanarStrobMap_vs_ePm_", order, ".eps")) #Save

#--------------------------
#Lab: energy
# Font
# Computer Modern Unicode Serif
# Latin Modern Roman
#--------------------------
# Start the tikz device
pLab = plotdf_point(tmdf_lab[which(tmdf_lab$dHw < 0.02),], "xEM", "yEM", "x [-]", "x [-]", "dHw", "", 0, pointSize = 1)
pLab = pLab+big_font_theme
pLab = pLab+ scale_colour_gradient2(TeX('$\\delta H(t)$'), space="Lab", midpoint = 0.008, mid = "green", high = "blue")
#pLab = pLab + scale_x_continuous(limits = c(-0.0025, 0.005))
#pLab = pLab + scale_y_continuous(limits = c(-0.01, 0.01))
#pLab = pLab+theme(text=element_text(family="CM Sans"))
pLab
ggsave(plot=pLab, width = 18, height = 18, file=paste0(currentfolder, "PlanarStrobMap_vs_dH.eps")) #Save


#Annotations
pLab = pLab + annotate("text", x = -0.0245, y = 0.04, label = "22", size = 10) 
pLab = pLab + geom_segment(aes(x = -0.024, y = 0.037, xend = -0.0207, yend = 0.0), 
                               colour = "black", 
                               arrow = arrow(length = unit(0.3, "cm"), 
                               type = "closed"))

ggsave(pLab, file=paste0(currentfolder, "PlanarStrobMap_vs_dH.pdf")) #Save

#--------------------------
# EM coordinates
#--------------------------
#Global
pY = plotdf_point(tmdf, "xEM", "yEM", "xEM", "yEM", "log10ePm", "ePm", 0, pointSize = 1)
pY = pY + scale_x_continuous(limits = c(-0.005, 0.005))
pY = pY + scale_y_continuous(limits = c(-0.01, 0.01))
pY = pY + geom_point(data = tmdf_prec, aes(x = xEM, y = yEM), color = "red")
pY


#Energy
pE = plotdf_point(tmdf, "xEM", "yEM", "xEM", "yEM", "dHw", "dHw", 0, pointSize = 1)
pE = pE + geom_point(data = tmdf0, aes(x = xEM, y = yEM), color = "white")
pE = pE + geom_point(data = tmdf_prec, aes(x = xEM, y = yEM), color = "red")
pE = pE + geom_point(data = tmdf_layer, aes(x = xEM, y = yEM), color = "green")
pE

#Figure 4.4 of Andreu: invariant manifolds of periodic orbits 23/25, for t = 0[T]
#Label
pL = plotdf_point(tmdf, "xEM", "yEM", "xEM", "yEM", "label", "label", 1, pointSize = 1)
pL

#--------------------------
# NC coordinates
#--------------------------
plotdf_point(tmdf, "x", "y", "x", "y", "log10ePm", "log10ePm", 0, pointSize = 1)

#--------------------------
# RCM coordinates
#--------------------------
pS = plotdf_point(tmdf, "s1", "s3", "s1", "s3", pointSize = 1)
pS = pS+geom_point(data = tmdf0, aes(x = s1, y = s3), color = "black", size = 4)+geom_point(data = tmdf0, aes(x = s1, y = s3), color = "white", size = 3)
pS






