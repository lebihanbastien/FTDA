############################################################
# R script to handle a precision map of the CRTBP
# Has to be used files such as e.g. Emap_RTBP_EML2.R
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
ttm_all = data.frame()
for (i in vorders)  #loop on the orders
{
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0(currentfolder, "Serv/eOm_", Type, "_ofs_",ofs_order,"_order_", toString(i))
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
#Get rid of the origin
isOrigin = ttm_all$s1 == 0 &
  ttm_all$s2 == 0 & ttm_all$s3 == 0 & ttm_all$s4 == 0
ttm_all = ttm_all[which(!isOrigin),]
# From NC to EM units
ttm_all = NCtoC(ttm_all, gamma)
# From EM to physical units
ttm_all = CtoPH(ttm_all, L)
# Radii from Li
ttm_all$rNC = sqrt(ttm_all$x^2+ttm_all$y^2+ttm_all$z^2)
ttm_all$rPH = sqrt(ttm_all$xPH^2+ttm_all$yPH^2+ttm_all$zPH^2)
#Select only some values
ttm_all = ttm_all[which(ttm_all$log10eOm <= 0),]

#------------------------------------------------------------------------------------
# Plot (Tile)
#------------------------------------------------------------------------------------
#------------------------------------------------
# Complet Plot (eOm = f(s1, s2))
#------------------------------------------------
Ri     = 15000;
dRi    = 400;
Rcolor = muted("green")
Rsize  = 0.2
dHzi   = 0.0025
ddHzi  = 1e-4

#Multiplot for all orders
#---------------------
ttm_l   = list();
ttm_lr  = list();
ttm_le  = list();
ttm_str = list();
ppl     = list();
for(i in 1:length(vorders))
{
  order = vorders[i];
  ttm_l[[i]]  = ttm_all[which(ttm_all$order == order),]      #select order
  ttm_lr[[i]] = ttm_l[[i]][which(abs(ttm_l[[i]]$rPH-Ri) < dRi),]  #select a certain radius within order
  ttm_le[[i]] = ttm_l[[i]][which(abs(ttm_l[[i]]$dHz-dHzi) < ddHzi),]  #select a certain energy within order
  strangeLine = abs(ttm_l[[i]]$log10eOm+6) < 3 & ttm_l[[i]]$s2 > 10
  ttm_str[[i]]= ttm_l[[i]][which(strangeLine),]
  ppl[[i]]    = plotdf_tile(ttm_l[[i]], si, sj, isLegendOn, -6)
  #ppl[[i]]    = ppl[[i]] + geom_tile(data = ttm_lr[[i]],  aes(s1, s2), color = Rcolor, alpha = 0, size = Rsize) #plot of one radius
  #ppl[[i]]    = ppl[[i]] + geom_tile(data = ttm_le[[i]],  aes(s1, s2), color = "red", alpha = 0, size = Rsize)  #plot of one energy
}
#Actual plot in multiplot format
pMult = multiplot(plotlist = ppl, cols = 2)

#Legend only (don't forget to display the legend in plotfd_tile!)
#--------------------
if(legendOnly == 1)
{
  dev.off()                    #clear all plots
  legend <- g_legend(ppl[[1]]) #select the legend
  pLegend = grid.draw(legend)  #displays it
}
#Save in eps file
#------------------------------------------------
if(legendOnly == 0)
{
  #plot the error map
  ggsave(ppl[[1]], file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, ".eps"))
}else{
  #plot the legend only
  ggsave(ppl[[1]], file=paste0(currentfolder, "legend_", Type, ".eps"))
}

#######################################################################################################
stop()
#######################################################################################################

#------------------------------------------------------------------------------------
# Plot (XY)
#------------------------------------------------------------------------------------
#------------------------------------------------
# Complet Plot (eOm = f(x, y))
#------------------------------------------------
xi = "xPH"
yi = "yPH"
xist = "x [km]"
yist = "y [km]"

primaryFactor = 1.0#-gamma*Le
legTitle = expression(log[10]~bgroup("(",e[O]~bgroup("(",frac(T,2),")"),")") ~~~~ "")


if(FWRK == "EM"){
  #Multiplot for all orders
  #---------------------
  ttm_l   = list();
  ppl     = list();
  for(i in 1:length(vorders))
  {
    #Plot
    order = vorders[i];
    ttm_l[[i]]  = ttm_all[which(ttm_all$order == order),]      #select order
    ppl[[i]]    = plotdf_point(ttm_l[[i]], xi, yi, xist, yist, "log10eOm", "log10eOm", 0, pointSize = 3)
    ppl[[i]]    = ppl[[i]]+scale_colour_gradient2(legTitle, space="Lab", midpoint = -6, mid = "white", high = muted("blue"), guide = FALSE)
    
    if(xi == "xPH")
    {
      #Earth direction
      ppl[[i]]    = ppl[[i]]+ annotate("text", x = sign(primaryFactor)*60000, y = -28000, label = "to Earth \ ", size = 6)
      ppl[[i]]    = ppl[[i]]+ geom_segment(aes(x = 55000, y = -25000, xend = 65000, yend = -25000), 
                                           colour = "black", 
                                           arrow = arrow(length = unit(0.3, "cm"), 
                                                         type = "closed"))
      #Add the moon
      primaryR    = 1737.10/abs(primaryFactor)  #Moon's radius
      primaryPos  = gamma*L/primaryFactor       #Moon position wrt to Li
      ppl[[i]]    = addMoon(ppl[[i]], x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2)
      ppl[[i]]    = ppl[[i]]+ annotate("text", x = sign(primaryFactor)*gamma*L/primaryFactor, y = 5000, label = "Moon", size = 6)
    }
  }
  #Actual plot in multiplot format
  pMult_ph = multiplot(plotlist = ppl, cols = 2)
}else
{
  #Multiplot for all orders
  #---------------------
  ttm_l   = list();
  ppl     = list();
  for(i in 1:length(vorders))
  {
    #Plot
    order = vorders[i];
    ttm_l[[i]]  = ttm_all[which(ttm_all$order == order),]      #select order
    ppl[[i]]    = plotdf_point(ttm_l[[i]], xi, yi, xist, yist, "log10eOm", "log10eOm", 0)
    ppl[[i]]    = ppl[[i]]+scale_colour_gradient2(legTitle, space="Lab", midpoint = -6, mid = "white", high = muted("blue"), guide = FALSE)
    
    
    #Add the Earth
    if(xi == "xPH")
    {
      primaryR    = 2*6478/abs(primaryFactor)  #Moon's radius
      primaryPos  = gamma*L/primaryFactor       #Moon position wrt to Li
      ppl[[i]]    = addMoon(ppl[[i]], x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2)
      ppl[[i]]    = ppl[[i]]+ annotate("text", x = sign(primaryFactor)*gamma*L/primaryFactor, y = 50000, label = "Earth", size = 6)
    }
  }
  #Actual plot in multiplot format
  pMult_ph = multiplot(plotlist = ppl, cols = 2)
}

#Save to eps file
ggsave(ppl[[1]], width= 9.86, height = 7.78, file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_withMoon_", xi, yi,".pdf"))
