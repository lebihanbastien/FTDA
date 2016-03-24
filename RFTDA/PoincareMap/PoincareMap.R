############################################################
# R script to handle a poincare map of the QBCP
# Has to be used files such as e.g. PoincareMap_QBCP_EML2.R
#
# 01/2016
############################################################

#------------------------------------------------
# Init
#------------------------------------------------
source("source/init.R")

#------------------------------------------------
#Current working folder
#------------------------------------------------
currentfolder = paste0(printfolder(MODEL, FWRK, Li))

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
muR = muR(FWRK);
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK)

#------------------------------------------------
#Size of plots
#------------------------------------------------
plotW = 18
plotH = 18

#------------------------------------------------------------------------------------
# Data reading
#------------------------------------------------------------------------------------
#------------------------------------------------
# Filename to check
#------------------------------------------------
#Old version
#filesuffix = paste0("Serv_pm_", add, "Energy_", Energy, "_order_", order, "_ofs_", ofs_order, METHOD);
#filesuffixnodots = paste0("Serv_pm_", add, "Energy_", nst, "_order_", order, "_ofs_", ofs_order, METHOD);
#New version
#filesuffix = paste0("Serv_pm_", add, "Energy_", Energy, "_order_", order, "_ofs_", ofs_order, "_proj_", projFreq, "_max_events_", maxEvents, METHOD);
#filesuffixnodots = paste0("Serv_pm_", add, "Energy_", nst, "_order_", order, "_ofs_", ofs_order, "_proj_", projFreq, "_max_events_", maxEvents, METHOD);
#Direct name
filename   = paste0(currentfolder, "Serv/", fileprefix, "/", fileprefix, fileext)

#------------------------------------------------
# Load csv source
#------------------------------------------------
if (file.exists(filename))
{
  if (fileext == ".txt")
  {
    pmdf  = read.csv(filename, header = T, sep = ",")
  }else if (fileext == ".bin")
  {
    #Read the data in binary form
    newdata = file(filename, "rb")
    datavals = readBin(newdata, double(), n = 18 * 100000)
    close(newdata)
    #Columns names
    dimnames <- list(
      row = seq(1,length(datavals) / 18),
      name = c(
        "label", "x", "y", "z", "px", "py", "pz",
        "s1", "s2", "s3", "s4", "t", "dHz", "dHw",
        "dHz.dHw", "pmap.dHv", "ePm", "number"
      )
    )
    #Create matrix from array
    matavals = matrix(datavals, ncol = 18, byrow = TRUE, dimnames = dimnames)
    #Create data.frame from matrix
    pmdf = data.frame(matavals)
  }else{
    #Send an error if the file extension is unknown
    stop('The file extension is unknown (not .txt nor .bin)')
  }
}else
{
  pmdf = data.frame()
}

#------------------------------------------------
#New fileprefix to avoid "." (dot) in outputs
#------------------------------------------------
filename = paste0(currentfolder, "Serv/", fileprefix, "/", fileprefixnodots)

#------------------------------------------------------------------------------------
# Postprocessing
#------------------------------------------------------------------------------------
#--------------------
# new columns
#--------------------
# From NC to EM units
pmdf = NCtoC(pmdf, gamma, c1)
# From EM to physical units
pmdf = CtoPH(pmdf, L)
# Radii from Li
pmdf$rNC = sqrt(pmdf$x ^ 2 + pmdf$y ^ 2 + pmdf$z ^ 2)
pmdf$rPH = sqrt(pmdf$xPH ^ 2 + pmdf$yPH ^ 2 + pmdf$zPH ^ 2)
#Parity of event
pmdf$parity = pmdf$number %% 2
#log(ePm)
pmdf$log10ePm = log10(pmdf$ePm)

#--------------------
#Starting points
#--------------------
pmdf0 = pmdf[which(pmdf$number == 0),]

#--------------------
# Selecting only values
# for which pz > 0! (true Poincaré map)
#--------------------
pmdf = pmdf[which(pmdf$pz >= 0),]

#--------------------
# Get rid of too big values of energy
# !! NEEDS IMPROVEMENT IN ORIGINAL C++ CODE !!
#--------------------
pmdf = pmdf[which(pmdf$dHw < 1.0),]

#--------------------
# Selecting particular values
#--------------------
#Maximum/Minimum time
maxnumber = max(pmdf$number)
minnumber = min(pmdf$number)

#--------------------
#Endpoints
#--------------------
pmdf1 = pmdf[which(pmdf$number == maxnumber),]

#--------------------
#Select some labels to "lighten" the plot
#--------------------
pmdfback = pmdf

#--------------------
#Find the halo orbits
#--------------------
if(Li == "L1")
{
  pmdfhaloindices = pmdfback[which(abs(pmdfback$yEM) > 0.057),]
}else{
  pmdfhaloindices = pmdfback[which(abs(pmdfback$yEM) > 0.07),]
}

#pmdfhalo   = pmdfback[which(pmdfback$label %in% pmdfhaloindices$label),]
#pmdfhalo_0 = pmdfhalo[which(pmdfhalo$number == 0),]

#--------------------
#Get rid of unwanted divergences
#--------------------
if(Li == "L2")
{
  pmdfdiv = pmdfback[which(pmdfback$xEM > -1.1),]
}

#--------------------
#Given condition
#--------------------
#condition = (pmdfback$label%%2 == 1 | pmdfback$label %in% pmdfhaloindices$label)
if(Li == "L2")
{
  condition = (pmdfback$label %% 2 == 1  | pmdfback$label %% 2 == 0 | pmdfback$label %in% pmdfhaloindices$label) & !(pmdfback$label %in% pmdfdiv$label)
  #condition = pmdfback$label %in% pmdfhaloindices$label
}else{
  condition = (pmdfback$label %% 2 == 1 | pmdfback$label %% 2 == 0 | pmdfback$label %in% pmdfhaloindices$label)
}
pmdf = pmdfback[which(condition),]

#--------------------
# Only a certain range of precision:
#--------------------
desiredPrec     = -6      #careful, if changed, pEM is obsolete!
desiredHardPrec = -3      #all kept solutions will have this precions

#--------------------
# pmdf_lab contains all complete solutions which the precision
# is better than desiredPrec  (for all points of the solution)
#--------------------
precFlag = pmdf$log10ePm < desiredHardPrec
wrongLabels = pmdf[which(!precFlag),]$label

#--------------------
# Make a difference between desiredPrec &  desiredHardPrec in pmdf_lab
#--------------------)
pmdf_lab = pmdf[which(!(pmdf$label %in% wrongLabels)),]
precFlag = pmdf_lab$log10ePm < desiredPrec
pmdf_lab$precFlag = precFlag


#------------------------------------------------------------------------------------
# Plots
#------------------------------------------------------------------------------------
#--------------------------
#Lab: vs ePm
#--------------------------
pEM = plotdf_point(pmdf_lab, "xEM", "yEM", "x [-]", "y [-]", "precFlag", expression(e[P](t) < 10 ^{-6}),1, pointSize = 2)
#Colors
pEM = pEM + scale_colour_manual(values = c("#000000", "#279B61"), guide = FALSE);#, "#CC6666", "#9999CC"), guide = FALSE) #to get the proper order!
#Theme
pEM = pEM + big_font_theme
#Limits if necessary
if (pEM_limits)
{
  pEM = pEM + scale_x_continuous(limits = pEM_limits_x)
  pEM = pEM + scale_y_continuous(limits = pEM_limits_y)
}
#Display
pEM
#Save
ggsave(
  pEM, width = plotW, height = plotH, file = paste0(filename, "_EM.eps")
) #Save



#--------------------------
# EM with prescribed precision
#--------------------------
pEMprec = plotdf_point(pmdf_lab, "xEM", "yEM", "x [-]", "y [-]", pointSize = 2)
#Theme
pEMprec = pEMprec + big_font_theme
#Display
pEMprec

#--------------------------
#Lab: labels
#--------------------------
pLab = plotdf_point(pmdf_lab, "xEM", "yEM", "x [-]", "y [-]", "label", "label", 1, pointSize = 2)
#Colors
pLab = pLab + scale_colour_discrete(guide = FALSE)
#Theme
pLab = pLab + big_font_theme
#Display
pLab
#Save
ggsave(
  pLab, width = plotW, height = plotH, file = paste0(filename, "_Labels.eps")
) #Save


#--------------------------
#Time between each projection
#--------------------------
pmdf_proj = ddply(pmdf, .(label, order), summarize, projTime = mean(diff(t, lag = 1)), x = mean(abs(x)), y = mean(abs(y)))
pProj = plotdf_point(pmdf_proj, "x", "y", "x", "y", "projTime", "projTime", 0, pointSize = 3)
if(Li=="L2")
{
  pProj = pProj + scale_colour_gradient2("projTime", space="Lab", midpoint = 1.74, mid = "white", high = muted("blue"))
}else{
  pProj = pProj + scale_colour_gradient2("projTime", space="Lab", midpoint = 2.76, mid = "white", high = muted("blue"))
}
pProj

pProjX = plotdf_point(pmdf_proj, "x", "projTime", "x", "projTime", pointSize = 1)
pProjX = pProjX + scale_y_continuous(limits=c(1.3,2))
pProjX

pProjY = plotdf_point(pmdf_proj, "y", "projTime", "y", "projTime", pointSize = 1)
pProjY = pProjY + scale_y_continuous(limits=c(1.3,2))
pProjY
stop()

#--------------------------
# PH coordinates
#--------------------------
#Global
pPH = plotdf_point(pmdf_lab, "xPH", "yPH", "x [km]", "y [km]", "log10ePm", "ePm", 0, pointSize = 1)
#Add the moon
primaryR    =  1737.10  #Moon's radius
primaryPos  =  (muR - 1) * L  #Moon position wrt to Li
pPH = addMoon(
  pPH, x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2
)
pPH = pPH + annotate(
  "text", x =  primaryPos - 500, y = 4000, label = "Moon", size = 10
)
#Theme
pPH = pPH + big_font_theme + coord_fixed(ratio = 1) + scale_color_continuous(guide = FALSE)
#Display
pPH
#Save
ggsave(
  pPH, width = plotW, height = plotH, file = paste0(filename, "_Moon.eps")
) #Save


#--------------------------
# RCM coordinates
#--------------------------
pRCM = plotdf_point(pmdf0, "s1", "s3", "s1", "s3", pointSize = 2)
#Theme
pRCM = pRCM + big_font_theme
#Display
pRCM + scale_y_continuous(breaks = seq(-40,40,10))
pRCM + scale_y_continuous(breaks = seq(-2,2,0.1)) + scale_x_continuous(breaks = seq(-2,2,0.1))

#--------------------------
# A particular solution
#--------------------------
pmdf_sol = pmdf[which]
