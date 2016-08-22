############################################################
# R script to handle a precision map of the QBCP/CRTBP
# Has to be used files such as e.g. Eplot_QBCP_EML1.R
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
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK)

#------------------------------------------------------------------------------------
# Building the data.frame of results: GS
#------------------------------------------------------------------------------------
ttm_all = data.frame()

# Filename to check
#------------------------------------------------
filename = paste0("~/BackUpBox/PhD/OOFTDA/fprint/QBCP/EM/L2/Serv/eOm_s1s2s3s4_ofs_30_order_6", ".txt")

# Load csv source
#------------------------------------------------
ttm_c  = read.csv(filename, header = T, sep = ",")
ttm_c$type = 'GS'

#rbind in ttm_all
ttm_all = rbind(ttm_all, ttm_c)

#------------------------------------------------------------------------------------
# Building the data.frame of results: NF
#------------------------------------------------------------------------------------
# Filename to check
#------------------------------------------------
filename = paste0("~/BackUpBox/PhD/OOFTDA/fprint/QBCP/EM/L2/Serv/eOm_NF_s1s2s3s4_ofs_30_order_6", ".txt")

# Load csv source
#------------------------------------------------
ttm_c  = read.csv(filename, header = T, sep = ",")
ttm_c$type = 'NF'

#rbind in ttm_all
ttm_all = rbind(ttm_all, ttm_c)

#------------------------------------------------------------------------------------
# Postprocessing
#------------------------------------------------------------------------------------
#--------------------
# New columns
#--------------------
# Compute -log10(precision)
ttm_all$log10eOm = log10(ttm_all$eOm)
# Get rid of the origin
isOrigin = ttm_all$s1 == 0 &
  ttm_all$s2 == 0 & ttm_all$s3 == 0 & ttm_all$s4 == 0
ttm_all = ttm_all[which(!isOrigin),]
# From NC to EM units
ttm_all = NCtoSYS(ttm_all, gamma, c1)
# From NC to C units
ttm_all = NCtoC(ttm_all, gamma)
# From C to physical units
ttm_all = CtoPH(ttm_all, L)
# Radii from Li
ttm_all$rNC = sqrt(ttm_all$x^2+ttm_all$y^2+ttm_all$z^2)
ttm_all$rEM = sqrt(ttm_all$xEM^2+ttm_all$yEM^2+ttm_all$zEM^2)
ttm_all$rCPH = sqrt(ttm_all$xCPH^2+ttm_all$yCPH^2+ttm_all$zCPH^2)

#------------------------------------------------
# Y label
#------------------------------------------------
ylabel = TeX("$log_{10}(e_O^1)$")

#------------------------------------------------
#Plot the precision as a function of the energy gap
#------------------------------------------------
xlabel  = TeX("$\\delta H_0$");
pE = plotdf_line(ttm_all, "dHz", "log10eOm", xlabel, ylabel, "type", "Type", 1)+legend_inside_theme

#------------------------------------------------
#Plot the precision as a function of the energy gap
#------------------------------------------------
pErCPH = plotdf_line(ttm_all, "dHz", "rCPH","dHz", "rCPH", "type", "Type", 1)+legend_inside_theme

#------------------------------------------------
# pXYCPH
#------------------------------------------------
pXYCPH = plotdf_point(ttm_all, "xCPH", "yCPH", "X", "Y", "type", "Type", 1)

#------------------------------------------------
# pXZCPH
#------------------------------------------------
pXZCPH = plotdf_point(ttm_all, "xCPH", "zCPH", "X", "Z", "type", "Type", 1)

#------------------------------------------------
# pXY
#------------------------------------------------
pXZ = plotdf_line(ttm_all, "xEM", "zEM", "xEM", "zEM", "type", "Type", 1)
#------------------------------------------------
# pXY
#------------------------------------------------
pXY = plotdf_point(ttm_all, "xEM", "yEM", "xEM", "yEM", "type", "Type", 1)
if(Li == "L2")
{
  #------------------------------------------------
  # Orbit 24 (rough param)
  #------------------------------------------------
  xc = c(-1.129, -1.13, -1.136, -1.14, -1.15, -1.16, -1.17, -1.175, -1.177)
  yc = c(0, 0.031, 0.060, 0.068, 0.069, 0.060, 0.040, 0.020, 0)
  df24 = data.frame(x = c(xc, xc),
                    y = c(yc, -yc))
  pXY = pXY + geom_point(data=df24, aes(x=x, y=y), size=3)
  pXY
  
  #------------------------------------------------
  # Divergence point in NF
  #------------------------------------------------
  xdiv = c(-1.133, -1.1668, -1.148, -1.148, -1.13, -1.1756)
  ydiv = c(0.049, -0.046, 0.066, -0.066, 0, 0)
  dfdiv = data.frame(x=xdiv, y = ydiv)
  pXY = pXY + geom_point(data=dfdiv, aes(x=x, y=y), size=3, color="red")
  pXY
}

#------------------------------------------------
#Save in eps file
#------------------------------------------------
ggsave(pS1, width = xSize, height = ySize,  file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_RCM.eps"))
ggsave(pX, width = 12.3, height = 8.92,   file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_EM.eps"))
ggsave(pPH, width = xSize, height = ySize,  file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_NC.eps"))
ggsave(pE, width = xSize, height = ySize,   file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order,  "_", Type, "_Energy.eps"))