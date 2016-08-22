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
# Building the data.frame of results
#------------------------------------------------------------------------------------
emap = data.frame()
for (i in vorders)  #loop on the orders
{
  # Filename to check
  #------------------------------------------------
  if(missing(Energy) || Energy == 0)
  {
    fileprefix = paste0(currentfolder, "Serv/eOm_", Type, "_ofs_",ofs_order,"_order_", toString(i))
  }else{
    fileprefix = paste0(currentfolder, "Serv/eOm_", Type, "_ofs_",ofs_order,"_order_", toString(i), "_hmax_", Energy)
  }
  filename = paste0(fileprefix, ".txt")
  
  
  # Load csv source
  #------------------------------------------------
  if (file.exists(filename))
  {
    ttm_c  = read.csv(filename, header = T, sep = ",")
  }else
  {
    ttm_c = data.frame()
  }
  
  #rbind in emap
  emap = rbind(emap, ttm_c)
}

#------------------------------------------------------------------------------------
# Postprocessing
#------------------------------------------------------------------------------------
#--------------------
# New columns
#--------------------
# Compute -log10(precision)
emap$log10eOm = log10(emap$eOm)
# Get rid of the origin
isOrigin = emap$s1 == 0 &
  emap$s2 == 0 & emap$s3 == 0 & emap$s4 == 0
emap = emap[which(!isOrigin),]
# From NC to EM units
emap = NCtoSYS(emap, gamma, c1)
# From NC to C units
emap = NCtoC(emap, gamma)
# From C to physical units
emap = CtoPH(emap, L)
# Radii from Li
emap$rNC = sqrt(emap$x^2+emap$y^2+emap$z^2)
emap$rEM = sqrt(emap$xEM^2+emap$yEM^2+emap$zEM^2)
emap$rCPH = sqrt(emap$xCPH^2+emap$yCPH^2+emap$zCPH^2)


#--------------------
# Remove values above 0
#--------------------
#emap = emap[which(emap$log10eOm < 0),]

#--------------------
# Reverse x-axis if necessary
#--------------------
if(emap$xEM[1] - emap$xEM[2] > 0)
{
  emap$xEM = -emap$xEM
}

#------------------------------------------------
# Labels
#------------------------------------------------
ylabel = "$\\log_{10}(E_O)$"
llabel = "Order"


#------------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------------
#Plot the precision as a function of rCPH
#------------------------------------------------
xlabel  = paste0("Distance from ", Lit, " [km]");
pPH = plotdf_line(emap, "rCPH", "log10eOm", xlabel, ylabel, "order", llabel, 1)
pPH = pPH + legend_inside_theme
pPH = pPH + scale_y_continuous(limits = pG_limits_y, breaks = pG_breaks_y)


#------------------------------------------------
#Plot the precision as a function of the energy gap
#------------------------------------------------
xlabel  = "$\\delta H_0$";
pE = plotdf_line(emap, "dHz", "log10eOm", xlabel, ylabel, "order", llabel, 1)
pE = pE + legend_inside_theme
pE = pE + scale_y_continuous(limits = pG_limits_y, breaks = pG_breaks_y)


#------------------------------------------------
#Plot the precision as a function of xEM
#------------------------------------------------
xlabel  = "$|X|$ $[$-$]$";
pX = plotdf_line(emap, "xEM", "log10eOm", xlabel, ylabel, "order", llabel, 1)
pX = pX + scale_y_continuous(limits = pG_limits_y, breaks = pG_breaks_y)
pX = pX + scale_x_continuous(limits = pX_limits_x)
if(emap$xEM[1] - emap$xEM[2] < 0)
{
  pX = pX + legend_left_theme
}else{
  pX = pX + legend_inside_theme
}


#--------------------------------------------------------------------------------------------------------------------------
#Save in tikz
#--------------------------------------------------------------------------------------------------------------------------
# ggplot2tikz(pX, xSize, ySize, file = paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", Type, "_EM.tex"))
# ggplot2tikz(pPH, width = xSize, height = ySize,  file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", Type, "_PH.tex"))
ggplot2tikz(pE, width = xSize, height = ySize,   file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order,  "_", Type, "_Energy.tex"))
