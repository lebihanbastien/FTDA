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

#------------------------------------------------
# Load Source files
#------------------------------------------------
source("source/source_folder.R")
source("source/source_plot.R")
source("source/source_routines.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li        = "L2"
MODEL     = "QBCP"
FWRK      = "EM"
METHOD    = "_SINGLE_INT"
Energy    = "0.0025"
order     = "20" #available for L2: 10,15,20
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
tmdf = data.frame()
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
  tmdf_c  = read.csv(filename, header = T, sep = ",")
}else
{
  tmdf_c = data.frame()
}
tmdf_c$steps = 5
tmdf = rbind(tmdf, tmdf_c)

#------------------------------------------------
# Filename to check
#------------------------------------------------
fileprefix = paste0(currentfolder, "Serv/Serv_tm_steps_2_Energy_", Energy, "_order_", order, "_ofs_", ofs_order, METHOD)
filename = paste0(fileprefix, ".txt")
#------------------------------------------------
# Load csv source
#------------------------------------------------
if (file.exists(filename))
{
  tmdf_c  = read.csv(filename, header = T, sep = ",")
}else
{
  tmdf_c = data.frame()
}
tmdf_c$steps = 2
tmdf = rbind(tmdf, tmdf_c)

#------------------------------------------------------------------------------------
# Postprocessing
#------------------------------------------------------------------------------------

#--------------------
# new columns
#--------------------
# From NC to EM units
tmdf = NCtoC(tmdf, gamma)
# From EM to physical units
tmdf = CtoPH(tmdf, L)
# Radii from Li
tmdf$rNC = sqrt(tmdf$x^2+tmdf$y^2+tmdf$z^2)
tmdf$rPH = sqrt(tmdf$xPH^2+tmdf$yPH^2+tmdf$zPH^2)
#Parity of event
tmdf$parity = tmdf$number%%2
#log(ePm)
tmdf$log10ePm = log10(tmdf$ePm)

#--------------------
# Selecting particular values
#--------------------
#Only full period
#tmdf = tmdf[which(tmdf$parity == 0),]

#Only a certain range of x
Xrange = tmdf$xEM < -0.01 | (tmdf$xEM > 0.01 & tmdf$xEM < 0.03 ) 
#tmdf = tmdf[which(Xrange),]

#Starting points
tmdf0 = tmdf[which(tmdf$number == 0),]

#Only a certain range of precision
precFlag = tmdf$log10ePm < -5 & tmdf$number > 0
tmdf_prec = tmdf[which(precFlag),]

#Select a given layer of energy
layer.value = 0.006
layer.width = 5e-4
tmdf_layer = tmdf[which(abs(tmdf$dHw - layer.value) < layer.width),]

#Only a certain range of s1 & s3
srange = tmdf0$s1 >= -30 & tmdf0$s1 <= 30 & tmdf0$s3 >= -30 & tmdf0$s3 <= 30
goodLabels = tmdf0$label[which(srange)]
tmdf_lab = tmdf[which(tmdf$label %in% goodLabels),]

#Only certain plots
#tmdf_lab = tmdf_lab[which(tmdf_lab$label %% 2 == 0),]

#Only certain precision
precFlag = tmdf_lab$log10ePm < -6 & tmdf_lab$number > 0
tmdf_lab$precFlag = precFlag
tmdf_lab_prec = tmdf_lab[which(precFlag),]

#Only the labels for which all points verifies precFlag
nonPrecFlag = tmdf_lab$log10ePm > -6 
badLabels  = tmdf_lab$label[which(nonPrecFlag)]
tmdf_lab_prec_strict = tmdf_lab[which(!(tmdf_lab$label %in% badLabels)),]
#Corresponding solution
plotdf_point(tmdf_lab_prec_strict, "xEM", "yEM", "xEM", "yEM")

#------------------------------------------------
# Plots
#------------------------------------------------
#--------------------------
#Lab: vs ePm
#--------------------------
pLab = plotdf_point(tmdf_lab, "xEM", "yEM", "xEM", "yEM", "precFlag", expression(e[P](t) < 10^{-6}),1)
pLab = pLab+big_font_theme
#pLab = pLab + geom_point(data = tmdf_lab_prec_strict, aes(x = xEM, y = yEM), color = "white") 
pLab
ggsave(pLab, width = 18, height = 18, file=paste0(currentfolder, "PlanarStrobMap_vs_ePm.eps")) #Save

#--------------------------
#Lab: labels
#--------------------------
pLab = plotdf_point(tmdf_lab, "xEM", "yEM", "xEM", "yEM", "label", "label", 1)
pLab = pLab+scale_colour_discrete(guide = FALSE)
pLab
#ggsave(pLab, width = 18, height = 18, file=paste0(currentfolder, "PlanarStrobMap_vs_label.eps")) #Save

#--------------------------
#Lab: energy
#--------------------------
pLab = plotdf_point(tmdf_lab, "xEM", "yEM", "xEM", "yEM", "dHw", "dH(t)", 0)
pLab = pLab+big_font_theme
#Annotations
pLab = pLab + annotate("text", x = -0.0245, y = 0.04, label = "22", size = 10) 
pLab = pLab + geom_segment(aes(x = -0.024, y = 0.037, xend = -0.0207, yend = 0.0), 
                               colour = "black", 
                               arrow = arrow(length = unit(0.3, "cm"), 
                               type = "closed"))
pLab
ggsave(pLab, width = 18, height = 18, file=paste0(currentfolder, "PlanarStrobMap_vs_dH.eps")) #Save

#ggsave(pLab, file=paste0(currentfolder, "PlanarStrobMap_vs_dH.eps")) #Save

#--------------------------
# EM coordinates
#--------------------------
#Global
pY = plotdf_point(tmdf, "xEM", "yEM", "xEM", "yEM", "log10ePm", "ePm", 0)
#pY = pY + scale_x_continuous(limits = c(-0.05, 0.05))
#pY = pY + scale_y_continuous(limits = c(-0.15, 0.15))
pY = pY + geom_point(data = tmdf_prec, aes(x = xEM, y = yEM), color = "red")
pY


#Energy
pE = plotdf_point(tmdf, "xEM", "yEM", "xEM", "yEM", "dHw", "dHw", 0)
pE = pE + geom_point(data = tmdf0, aes(x = xEM, y = yEM), color = "white")
#pE = pE + geom_point(data = pmdf_prec, aes(x = xEM, y = yEM), color = "red")
pE = pE + geom_point(data = tmdf_layer, aes(x = xEM, y = yEM), color = "green")
pE

#Figure 4.4 of Andreu: invariant manifolds of periodic orbits 23/25, for t = 0[T]
#Label
pL = plotdf_point(tmdf, "xEM", "yEM", "xEM", "yEM", "label", "label", 1)
pL

#--------------------------
# NC coordinates
#--------------------------
plotdf_point(tmdf, "x", "y", "x", "y", "rNC", "rNC", 0)

#--------------------------
# RCM coordinates
#--------------------------
pS = plotdf_point(tmdf, "s1", "s3", "s1", "s3")
pS = pS+geom_point(data = tmdf0, aes(x = s1, y = s3), color = "black", size = 4)+geom_point(data = tmdf0, aes(x = s1, y = s3), color = "white", size = 3)
pS



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
#Select the label for which the sd is maximum
#--------------------------
dHw.sd.max.label = tmdfe$label[which(tmdfe$dHw.sd == max(tmdfe$dHw.sd))]
tmdf_en = tmdfe[which(tmdfe$label == dHw.sd.max.label[1]),]
#Plot the corresponding energy variations
pEn = plotdf_path(tmdf_en, "t", "dHw", "\n time [adim]", "dH(t) [adim] \n", "label", "label", 0)
pEn = pEn + scale_color_continuous(guide = FALSE)
pEn

#--------------------------
#Select a random label
#--------------------------
tmdf_en = tmdfe[which(tmdfe$label == 3),]
#Plot the corresponding energy variations
pEn = plotdf_path(tmdf_en, "t", "dHw", "\n time [adim]", "dH(t) [adim] \n", "label", "label", 0)
pEn = pEn + scale_color_continuous(guide = FALSE)
pEn

#Corresponding solution
plotdf_point(tmdf_en, "xEM", "yEM", "xEM", "yEM", "label", "label", 1)

#--------------------------
# Computation standard deviation & mean
#--------------------------
tmdfe_sum = ddply(tmdf, .(label, order, steps), summarize, dHw.sdpc = sd(dHw)/mean(dHw), dHw.sd = sd(dHw), dHw.mean = mean(dHw), maxT = max(t), maxN = max(number))

pSd = plotdf_point(tmdfe_sum , "dHw.mean", "dHw.sd", expression(paste("\n ", mu(dH))), expression(paste(sigma(dH), "\n")), "steps", "proj/period", 1)
pSd
ggsave(pSd, file=paste0(currentfolder, "PlanarStrobMap_sd_vs_mean_order_Steps.eps")) #Save

pSd = plotdf_point(tmdfe_sum , "dHw.mean", "dHw.sdpc", 
                   expression(paste("\n ", mu(dH))), 
                   expression(paste(frac(sigma(dH), mu(dH)), "\n")),
                   "steps", "proj/period", 1)
pSd = pSd + scale_y_continuous(limits = c(0, 0.35))
pSd = pSd + theme(axis.title.y = element_text(angle = 0))
pSd
ggsave(pSd, file=paste0(currentfolder, "PlanarStrobMap_sdmean_vs_mean_Steps.eps")) #Save
