# R script to handle a precision map of the QBFBP around EML1,2
#---------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------
# Init
#------------------------------------------------
source("source/init.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li    = "L2"
Lit    = "$L_2$"
MODEL = "QBCP"
FWRK  = "EM"
Type  = "rand" #selection or global
Energy = 0
vorders = c(5, 10, 15, 20, 25, 30);
vorders = c(5, 7, 10, 12, 15, 17, 20, 22,  25, 27, 30);
currentfolder = paste0(printfolder(MODEL, FWRK, Li))

#------------------------------------------------
#Normalized units (gamma, c1)
#------------------------------------------------
muR = muR(FWRK);
gamma = gamma(Li, FWRK);
c1    =  c1(Li, FWRK);
L     = Ldist(FWRK);
if(FWRK == "EM")
{
  primaryR    =  1737.10      #m2
}else{
  primaryR    =  6378.10      #m2 
}

#------------------------------------------------
#Additionnal parameters
#------------------------------------------------
isLegendOn = 1;
legendOnly = 0;

#------------------------------------------------
# Building the data.frame of results
#------------------------------------------------
imap = data.frame()
for (i in vorders)  #loop on the orders
{
  if(Energy == 0)
  {
    fileprefix = paste0(currentfolder, "Serv/eIm_", Type, "_ofs_30_order_",toString(i));
    #fileprefix = paste0(currentfolder, "Serv/eIm_", Type, "_ofs_30_order_",toString(i), "_t0_3.3956");
  }else
  {
    fileprefix = paste0(currentfolder, "Serv/eIm_", Type, "_ofs_30_order_",toString(i), "_energy_", toString(Energy));
  }
  filename = paste0(fileprefix, ".bin")
  
  
  # Load csv source
  #------------------------------------------------
  if (file.exists(filename))
  {
    names = c("label", "x", "y", "z", "px", "py", "pz",
              "xEM", "yEM", "zEM", "pxEM", "pyEM", "pzEM",
              "s1", "s2", "s3", "s4", "t", "dHz", "eIm");
    imapc = dffbinary(filename, 20, names);
    imapc$order = i;
    
  }else
  {
    imapc = data.frame()
  }
  
  #rbind in ttm_all
  imap = rbind(imap, imapc);
}

#------------------------------------------------
# Postprocessing
#------------------------------------------------
# Using EIm (EM units) instead of eIm (in NC units) 
imap$EIm = gamma*imap$eIm 
# To centered units
imap = NCtoC(imap, gamma)
# To physical units
imap = CtoPH(imap, L)
# To physical units
imap = EMtoPH(imap, L)
# Compute -log10(precision)
imap$log10eOm = log10(imap$EIm)
# Compute -log10(precision)
imap$flog10eOm = floor(log10(imap$EIm))
#Distance to the center
imap$rC = sqrt(imap$xC^2+imap$yC^2+imap$zC^2)
imap$rCPH = sqrt(imap$xCPH^2+imap$yCPH^2+imap$zCPH^2)
imap$frCPH = floor(imap$rCPH*1e-4)
#Abs of zC
imap$zCn = abs(imap$zC)

#------------------------------------------------
# Get rid of bad precision
#------------------------------------------------
imap = imap[which(imap$log10eOm < -1),]

#------------------------------------------------
# Select only positive some given value in the results
#------------------------------------------------
imapsys = data.frame();
# veps = c(seq(1,8,2)*1e-5, seq(1,8,2)*1e-6, seq(1,8,2)*1e-7, seq(1,8,2)*1e-8, seq(1,8,2)*1e-9);
if(Li == "L2")
{
  veps = 10^(-seq(5,9));
}else{
  veps = 10^(-seq(5,9));
}
for (eps in veps)  #loop on the orders
{
  isCloseToEPS = abs(imap$EIm - eps) < 0.2*eps;
  imapr = imap[which(isCloseToEPS),]
  
  #Ordered version
  imapr = imapr[order(imapr$order),]
  
  # Compute the mean position
  imapm = ddply(imapr, .(order), summarize, mdHz = mean(dHz), mrCPH = mean(rCPH), mrC = mean(rC), mxC = mean(abs(xC)), myC = mean(abs(yC)),  mzC = mean(abs(zC)))
  imapm$eps = log10(eps)
  
  
  #rbind in ttm_all
  imapsys = rbind(imapsys, imapm);
}
#imapsys = imapsys[order(imapsys$order),]

#------------------------------------------------
# Compute the mean position
#------------------------------------------------
ppmd = plotdf_line(imapsys, "mrCPH", "order", paste0("mean distance to ", Lit), "Order", colorCol = "eps", colorLabel = "Precision", isColorFac = TRUE)
ppmd = ppmd + scale_y_continuous(breaks = seq(5,30,5))
ppmd

# To put in the tikz file for legend title : \\raisebox{0.3ex}{\\scriptsize{$\\log_{10}(E_I)$}}
ppme = plotdf_path(imapsys, "mdHz", "order", "$\\delta \\bar{H}_0$", "Order $N$", colorCol = "eps", colorLabel = "$\\raisebox{0.3ex}{\\scriptsize{$\\log_{10}(E_I)$}}$", isColorFac = TRUE)
ppme = ppme + scale_y_continuous(breaks = seq(5,30,5))
ppme = ppme + legend_inside_theme# + theme(legend.title = element_text(margin = margin(10,10,10,50)))
ppme

ppmer = plotdf_line(imapsys, "mdHz", "eps", "mean($\\delta H_0$)", "$log_{10}(e_I)$", colorCol = "order", colorLabel = "Order", isColorFac = TRUE)
ppmer = ppmer + scale_y_continuous(breaks = seq(-10,-1,1))
ppmer

# plotdf_smooth(imapm, "mxC", "order", paste0("mean $X$ to ", Lit), "Order")
# plotdf_smooth(imapm, "myC", "order", paste0("mean $Y$ to ", Lit), "Order")
# plotdf_smooth(imapm, "mzC", "order", paste0("mean $Z$ to ", Lit), "Order")

#---------------------------------------------------------------------------------------------------------------------------
# Plot 
#---------------------------------------------------------------------------------------------------------------------------
#---------------------
# Only a given precision
#---------------------
pplc = plotdf_point(imapr, "xEM", "yEM", "X", "Y", pointSize = 3, colorCol = "order", colorLabel = "order", isColorFac = TRUE)
pplc
#Add m2
primaryPos  =  -(1-muR)  #Add m2
pplc = addMoon(pplc, x = primaryPos, y = 0, primaryR/L, surfSize = 0.4, cratSize = 0.2
)
pplc

pplz = plotdf_point(imapr, "xEM", "zEM", "X", "Z", pointSize = 3, colorCol = "order", colorLabel = "order", isColorFac = TRUE)
pplz

#---------------------
# Adding the resonnances
#---------------------
if(Li == "L2")
{
  #------------------------------------------------
  # Orbit 24 (rough param)
  #------------------------------------------------
  xc  = c(-1.129, -1.13, -1.136, -1.14, -1.15, -1.16, -1.17, -1.175, -1.177)
  xc2 = c(-1.177, -1.175, -1.17, -1.16, -1.15, -1.14, -1.136, -1.13, -1.129)
  yc  = c(0, 0.031, 0.060, 0.068, 0.069, 0.060, 0.040, 0.020, 0)
  yc2 = c(0, 0.020, 0.040, 0.060, 0.069, 0.068, 0.060, 0.031, 0)
  df24 = data.frame(x = c(xc,  xc2),
                    y = c(yc, -yc2))
  pplc = pplc + geom_path(data=df24, aes(x=x, y=y), size=3)
}

#---------------------
# Only a given precision, in physical units
#---------------------
pplcph = plotdf_point(imapr, "xCPH", "yCPH", "X [km]", "Y [km]", pointSize = 3, colorCol = "order", colorLabel = "order", isColorFac = TRUE)
pplcph

pplph = plotdf_point(imapr, "xPH", "yPH", "X [km]", "Y [km]", pointSize = 3, colorCol = "order", colorLabel = "order", isColorFac = TRUE)
#Add m2
primaryPos  =  -(1-muR) * L #Add m2
pplph = addMoon(pplph, x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2
)
pplph

#---------------------
# Only a given order, at a given z
#---------------------
# y vs x, plot the precision
ppl = plotdf_point(imap[which(imap$order == 20),], "xEM", "yEM", "X", "Y", pointSize = 3, colorCol = "flog10eOm", colorLabel = "log10eOm", isColorFac = TRUE)
ppl

# z vs precision, plot the radius
ppz = plotdf_point(imap[which(imap$order == 30),], "zCn", "log10eOm", "|zC|", "log10eOm", pointSize = 3, colorCol = "frCPH", colorLabel = "rC", isColorFac = TRUE)
ppz

# z vs precision, plot the energy
pph = plotdf_point(imap[which(imap$order == 20),], "zCn", "log10eOm", "|zC|", "log10eOm", pointSize = 3, colorCol = "dHz", colorLabel = "dHz", isColorFac = FALSE)
pph + scale_color_gradient2(space="Lab", midpoint = 0.02, mid = "white", high = muted("red"), low = muted("blue"))

# z vs x, plot the precision
ppr = plotdf_point(imap[which(imap$order == 20),], "xC", "zC", "rC", "zC", pointSize = 3, colorCol = "log10eOm", colorLabel = "log10eOm", isColorFac = FALSE)
ppr + scale_color_gradient2(space="Lab", midpoint = -6, mid = "white", high = muted("blue"))

#---------------------
# With respect to the energy
#---------------------
#Select only some given value in the results
eps = 0.015;
isCloseToEPS = abs(imap$dHz - eps) < 1e-4;
imape = imap[which(isCloseToEPS),]
ppe = plotdf_point(imap[which(imap$order == 20),], "xEM", "yEM", "X", "Y", pointSize = 3, colorCol = "flog10eOm", colorLabel = "log10eOm", isColorFac = TRUE)
ppe

#---------------------------------------------------------------------------------------------------------------------------
# Adding the CR3BP zero-velocity curves
#---------------------------------------------------------------------------------------------------------------------------
if(Energy > 0)
{
  #---------------------
  #Build vectors
  #---------------------
  x0   = seq(-2,-1,0.001)
  y0   = seq(-0.15,0.15,0.001)
  
  z0   = expand.grid(x0, y0)
  potn = apply(z0, 1, cr3bpot_EM)
  
  
  #Build data frame
  df = data.frame(x1 = z0[1], y1 = z0[2], pot = potn)
  names(df) = c('x', 'y', 'pot')
  
  #---------------------
  #Energy associated with the Libration points
  #---------------------
  l1 = c(-0.836915145386502, 0, 0)
  l2 = c(-1.155682150113637, 0, 0)
  l1pot = cr3bpot(l1[1], l1[2], muR)
  l2pot = cr3bpot(l2[1], l2[2], muR)
  
  #---------------------
  #Select a given energy
  #---------------------
  eps = l2pot+Energy;
  isCloseToEPS = abs(df$pot - eps) < 1e-4;
  dfr = df[which(isCloseToEPS),]
  
  #---------------------
  #Plot
  #---------------------
  pplcp = pplc + geom_point(data=dfr, aes(x=x, y=y), size=3)
  pplcp+ coord_fixed(ratio=1)
}


#---------------------------------------------------------------------------------------------------------------------------
# 3D plot
#---------------------------------------------------------------------------------------------------------------------------
# scatter3D(imapr$xEM, imapr$yEM, imapr$zEM, colvar = imapr$order, pch = 16, cex = 1.5)
# scatter3D(imap$xEM, imap$yEM, imap$zEM, colvar = imap$flog10eOm, pch = 16, cex = 1.5)

#Save in tikz
#--------------------------------------------------------------------------------------------------------------------------
# ggplot2tikz(ppme, xSize, ySize, file = paste0(currentfolder, "EIm_", Type, "_ofs_30_order_",toString(i), "_3D.tex"))
