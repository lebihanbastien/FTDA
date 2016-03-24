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
#------------------------------------------------
#Absolute energy associated with the Li
#------------------------------------------------
if(Li == "L2")
{
  HLi = -1.586080214753386e+00
}else
{
  HLi = -1.594170540495413e+00
}

#------------------------------------------------------------------------------------
# Building the data.frame of results
#------------------------------------------------------------------------------------
ttm_all = data.frame()
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
  
  #rbind in ttm_all
  ttm_all = rbind(ttm_all, ttm_c)
}

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
ttm_all = NCtoC(ttm_all, gamma, c1)
# From EM to physical units
ttm_all = CtoPH(ttm_all, L)
# Radii from Li
ttm_all$rNC = sqrt(ttm_all$x^2+ttm_all$y^2+ttm_all$z^2)
ttm_all$rEM = sqrt(ttm_all$xEM^2+ttm_all$yEM^2+ttm_all$zEM^2)
ttm_all$rPH = sqrt(ttm_all$xPH^2+ttm_all$yPH^2+ttm_all$zPH^2)
#Add the absolute energy
ttm_all$dHz = ttm_all$dHz +HLi

#--------------------
# Select only positive s1/s3
#--------------------
if(Type == "s3"){
  ttm_mid = ttm_all[which(ttm_all$s3 <=0),]
}else{
  ttm_mid = ttm_all[which(ttm_all$s1 >=0),]
}

#ttm_all = ttm_mid
#--------------------
# Remove values above 0
#--------------------
#ttm_all = ttm_all[which(ttm_all$log10eOm < 0),]

#------------------------------------------------
# Y label
#------------------------------------------------
ylabel = TeX("$log_{10}(e_O^1)$")

#------------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------------
#Plot the precision as a function of rPH
#------------------------------------------------
xlabel  = paste0("Distance from ", Li, " [km]");
pPH = plotdf_line(ttm_mid, "rPH", "log10eOm", xlabel, ylabel, "order", "Order", 1)+legend_inside_theme
pPH


#Plot the precision as a function of rPH
#------------------------------------------------
xlabel  = paste0("Distance from ", Li, " [-]");
pEM = plotdf_line(ttm_mid, "rEM", "log10eOm", xlabel, ylabel, "order", "Order", 1)+legend_inside_theme
pEM

#------------------------------------------------
#Plot the precision as a function of s1
#------------------------------------------------
xlabel  = expression(s[1]);
pS1 = plotdf_line(ttm_all, "s1", "log10eOm", xlabel, ylabel, "order", "Order", 1)+legend_inside_theme
if(pS1_limits)
{
  pS1 = pS1 + scale_x_continuous(limits = pS1_limits_x)
  pS1 = pS1 + scale_y_continuous(limits = pS1_limits_y)
}

#------------------------------------------------
#Plot the precision as a function of s2
#------------------------------------------------
xlabel  = expression(s[2]);
pS2 = plotdf_line(ttm_all, "s2", "log10eOm", xlabel, ylabel, "order", "Order", 1)


#------------------------------------------------
#Plot the precision as a function of s3
#------------------------------------------------
xlabel  = expression(s[3]);
pS3 = plotdf_line(ttm_all, "s3", "log10eOm", xlabel, ylabel, "order", "Order", 1)



#------------------------------------------------
#Plot the precision as a function of the energy gap
#------------------------------------------------
xlabel  = TeX("$\\delta H_0$");
pE = plotdf_line(ttm_mid, "dHz", "log10eOm", xlabel, ylabel, "order", "Order", 1)+legend_inside_theme
#Postprocessing
pE = pE + big_font_theme
pE = pE + theme(axis.title.x=element_text(vjust=0.01))
#Breaks & limits
if(Li == "L1") pE = pE + scale_x_continuous(breaks = seq(-1.60, -1.4, 0.01), limits = c(NA,-1.54))
if(Li == "L2") pE = pE + scale_x_continuous(breaks = seq(-1.60, -1.4, 0.002))



#------------------------------------------------
#Plot the precision as a function of xEM
#------------------------------------------------
xlabel  = "\n x [-]";
pX = plotdf_line(ttm_all, "xEM", "log10eOm", xlabel, ylabel, "order", "Order", 1)
if(pX_limits)
{
  pX = pX + scale_x_continuous(limits = pX_limits_x)
  pX = pX + scale_y_continuous(limits = pX_limits_y)
}

if(pX_breaks_x)
{
  pX = pX + scale_x_continuous(breaks = pX_breaks_x_values)
}

if(pX_breaks_y)
{
  pX = pX + scale_y_continuous(breaks = pX_breaks_y_values)
}

if(pX_limits && pX_breaks_x && pX_breaks_y)
{
  pX = pX + scale_x_continuous(limits = pX_limits_x, breaks = pX_breaks_x_values)
  pX = pX + scale_y_continuous(limits = pX_limits_y, breaks = pX_breaks_y_values)
}
pX = pX + legend_inside_theme

#------------------------------------------------
#Plot the precision as a function of yEM
#------------------------------------------------
xlabel  = "\n y [-]";
pY = plotdf_line(ttm_all, "yEM", "log10eOm", xlabel, ylabel, "order", "Order", 1)
if(pY_limits && pY_breaks_x && pY_breaks_y)
{
  pY = pY + scale_x_continuous(limits = pY_limits_x, breaks = pY_breaks_x_values)
  pY = pY + scale_y_continuous(limits = pY_limits_y, breaks = pY_breaks_y_values)
}

if(pY_breaks_x)
{
  pY = pY + scale_x_continuous(breaks = pY_breaks_x_values)
}

if(pY_breaks_y)
{
  pY = pY + scale_y_continuous(breaks = pY_breaks_y_values)
}


pY

#------------------------------------------------
#Plot the precision as a function of yEM
#------------------------------------------------
xlabel  = expression("zEM");
pZ = plotdf_line(ttm_all, "zEM", "log10eOm", xlabel, ylabel, "order", "Order", 1)


#------------------------------------------------
#Plot the energy gap as a function of s1
#------------------------------------------------
plotdf_line(ttm_all, "s1", "dHz", expression(s[1]), "Energy gap from the origin", "order", "Order", 1)
plotdf_line(ttm_all, "s3", "dHz", expression(s[3]), "Energy gap from the origin", "order", "Order", 1)
plotdf_line(ttm_all, "xEM", "dHz", expression(x[EM]), "Energy gap from the origin", "order", "Order", 1)

#------------------------------------------------
#yPH as a function of xPH
#------------------------------------------------
plotdf_line(ttm_all, "xPH", "yPH", "xPH", "yPH", "order", "Order", 1)
pXZ = plotdf_line(ttm_all, "xEM", "zEM", "xEM", "zEM", "order", "Order", 1)
pS1S2 = plotdf_line(ttm_all, "s1", "s2", "s1", "s2", "order", "Order", 1)

pXY = plotdf_point(ttm_all, "xEM", "yEM", "xEM", "yEM", "order", "Order", 1)

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


#--------------------
# Error on the order
#--------------------
errorOrder = data.frame();

if(Type == "s3")
{
  #Loop on all orders
  for(h in vorders)
  {
    #Select order h
    ttm_o = ttm_all[which(ttm_all$order == h),]
    ttm_o = ttm_o[order(ttm_o$s3),]
    
    #Build the n array  
    n = NULL;
    for (i in seq(1, length(ttm_o$s3)-1))
    {
      n[i] = log(ttm_o$eOm[i]/ttm_o$eOm[i+1])/log(ttm_o$s3[i]/ttm_o$s3[i+1])
    }
    #Build the corresponding dataframe
    errorOrder_c = data.frame(n = n, s3 = ttm_o$s3[1:length(ttm_o$s3)-1], order = h);
    #Bind
    #rbind in ttm_all
    errorOrder = rbind(errorOrder, errorOrder_c)
  }
  pOrder = plotdf_line(errorOrder, "s3", "n", expression(s[3]), "n", "order", "Order", 1)+legend_inside_theme
  pOrder = pOrder + scale_y_continuous(breaks = c(0, 1, 10, 11, 15, 16, 20, 21, 25, 26, 30, 31))
}else{
  
  #Loop on all orders
  for(h in vorders)
  {
    #Select order h
    ttm_o = ttm_all[which(ttm_all$order == h),]
    ttm_o = ttm_o[order(ttm_o$s1),]
    
    #Build the n array  
    n = NULL;
    for (i in seq(1, length(ttm_o$s1)-1))
    {
      n[i] = log(ttm_o$eOm[i]/ttm_o$eOm[i+1])/log(ttm_o$s1[i]/ttm_o$s1[i+1])
    }
    #Build the corresponding dataframe
    errorOrder_c = data.frame(n = n, s1 = ttm_o$s1[1:length(ttm_o$s1)-1], order = h);
    #Bind
    #rbind in ttm_all
    errorOrder = rbind(errorOrder, errorOrder_c)
  }
  pOrder = plotdf_line(errorOrder, "s1", "n", expression(s[1]), "n", "order", "Order", 1)+legend_inside_theme
  #pOrder = pOrder + scale_y_continuous(breaks = seq(1,40))
  pOrder = pOrder + scale_y_continuous(breaks = c(0, 1, 10, 11, 15, 16, 20, 21, 25, 26, 30, 31))
  
}