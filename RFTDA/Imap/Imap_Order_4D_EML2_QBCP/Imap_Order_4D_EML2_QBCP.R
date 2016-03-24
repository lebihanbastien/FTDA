#------------------------------------------------
# Import data (just once)
#------------------------------------------------
#source('Imap/Imap_Order_4D_EML2_QBCP/Imap_Order_4D_EML2_QBCP_HEADER.R')

#------------------------------------------------
# Select the plots
#------------------------------------------------
si = "s2"
sj = "s4"
vorders = (10)

#------------------------------------------------
# Find all positive values of s1
#------------------------------------------------
svec = sort(unique(ttm_all$s1))
svec = svec[which(svec >= 0)]
lvec = length(svec)

#------------------------------------------------
# Perform cuts
#------------------------------------------------
#Current indix
ind = 9
#Cut
ttm_cut = imap_4d_cut(si, sj, ind)

#------------------------------------------------
# Plot (Tile)
#------------------------------------------------
#------------------------------------------------
# Complet Plot (eOm = f(s1, s2))
#------------------------------------------------
legTitle = expression(log[10]~bgroup("(",e[I]~bgroup("(",0,")"),")") ~~~~ "")

#Multiplot for all orders
#---------------------
ttm_l   = list();
ppl     = list();
for(i in 1:length(vorders))
{
  order = vorders[i];
  ttm_l[[i]]  = ttm_cut[which(ttm_cut$order == order),]      #select order
  ppl[[i]]    = plotdf_tile(ttm_l[[i]], si, sj, isLegendOn, -8, title = legTitle, legendOnTop = legendOnly)
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
  ggsave(ppl[[1]], file=paste0(currentfolder, "R_eIm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_", ind, ".eps"))
}else{
  #plot the legend only
  ggsave(pLegend, file=paste0(currentfolder, "R_eIm_legend_", Type, ".eps"))
}

stop()

#------------------------------------------------
# Plot (XY)
#------------------------------------------------
#------------------------------------------------
# Complet Plot (eOm = f(x, y))
#------------------------------------------------
xi = "xPH"
yi = "yPH"
xist = "x [km]"
yist = "y [km]"

primaryFactor = 1.0#-gamma*Le



if(FWRK == "EM"){
  #Multiplot for all orders
  #---------------------
  ttm_l   = list();
  ppl     = list();
  for(i in 1:length(vorders))
  {
    #Plot
    order = vorders[i];
    ttm_l[[i]]  = ttm_cut[which(ttm_cut$order == order),]      #select order
    ppl[[i]]    = plotdf_point(ttm_l[[i]], xi, yi, xist, yist, "log10eOm", "log10eOm", 0)
    ppl[[i]]    = ppl[[i]]+scale_colour_gradient2(legTitle, space="Lab", midpoint = -6, mid = "white", high = muted("blue"), guide = FALSE)
    
    if(xi == "xPH")
    {
      #Earth direction
      ppl[[i]]    = ppl[[i]]+ annotate("text", x = -sign(primaryFactor)*4000, y = -28000, label = "to Earth \ ", size = 6)
      ppl[[i]]    = ppl[[i]]+ geom_segment(aes(x = -1000, y = -25000, xend = -9000, yend = -25000), 
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
    ttm_l[[i]]  = ttm_cut[which(ttm_cut$order == order),]      #select order
    ppl[[i]]    = plotdf_point(ttm_l[[i]], xi, yi, xist, yist, "log10eOm", "log10eOm", 0)
    ppl[[i]]    = ppl[[i]]+scale_colour_gradient2(legTitle, space="Lab", midpoint = -6, mid = "white", high = muted("blue"), guide = FALSE)
    
    #Earth direction
    #ppl[[i]]    = ppl[[i]]+ annotate("text", x = -sign(primaryFactor)*4000, y = -28000, label = "to Earth \ ", size = 6)
    #ppl[[i]]    = ppl[[i]]+ geom_segment(aes(x = -1000, y = -25000, xend = -9000, yend = -25000), 
    #                                     colour = "black", 
    #                                     arrow = arrow(length = unit(0.3, "cm"), 
    #                                                   type = "closed"))
    
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
ggsave(ppl[[1]], file=paste0(currentfolder, "R_eOm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_withMoon_", xi, yi,".eps"))


stop()


#------------------------------------------------
# Plot (dH vs XY)
#------------------------------------------------
#------------------------------------------------
# Complet Plot (eOm = f(x, y))
#------------------------------------------------
xi = "xPH"
yi = "yPH"
xist = "x [km]"
yist = "y [km]"

primaryFactor = 1.0#-gamma*Le

#Multiplot for all orders
#---------------------
ttm_l   = list();
ppl     = list();
for(i in 1:length(vorders))
{
  #Plot
  order = vorders[i];
  ttm_l[[i]]  = ttm_cut[which(ttm_cut$order == order),]      #select order
  ppl[[i]]    = plotdf_point(ttm_l[[i]], xi, yi, xist, yist, "dHz", "dHz", 0)
  ppl[[i]]    = ppl[[i]]+scale_colour_gradient2(legTitle, space="Lab", midpoint = 0.01, mid = "white", high = muted("blue"))
  
  if(xi == "xPH")
  {
    #Earth direction
    ppl[[i]]    = ppl[[i]]+ annotate("text", x = -sign(primaryFactor)*4000, y = -28000, label = "to Earth \ ", size = 6)
    ppl[[i]]    = ppl[[i]]+ geom_segment(aes(x = -1000, y = -25000, xend = -9000, yend = -25000), 
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


