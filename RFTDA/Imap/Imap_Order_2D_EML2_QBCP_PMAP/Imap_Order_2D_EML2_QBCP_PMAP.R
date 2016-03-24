#------------------------------------------------
# Import data (just once)
#------------------------------------------------
source('Imap/Imap_Order_2D_EML2_QBCP_PMAP/Imap_Order_2D_EML2_QBCP_PMAP_HEADER.R')

#------------------------------------------------
#Additionnal parameters
#------------------------------------------------
isLegendOn  = 0;
legendOnly  = 0;
#vorders = (30) #redefine the vector of orders

#------------------------------------------------
#Color gradient
#------------------------------------------------
colorCol      = "log10eOm";
colorMidPoint = -8;
colorLabel    = expression(log[10]~bgroup("(",e[I]~bgroup("(",0,")"),")") ~~~~ "")

#------------------------------------------------
# Postprocessing
#------------------------------------------------
#Select very high values of too high values (artifacts)
ttm_bad = ttm_all[which(ttm_all$log10eOm > 0),] 
#Get rid of these values in the current df
ttm_all = ttm_all[which(ttm_all$log10eOm < 0),] 
#Select some given values
#value = ttm_all$s1 == 14 & ttm_all$s3 == -2
#ttm_all = ttm_all[which(value),]

#------------------------------------------------
# Plot (Tile)
#------------------------------------------------
#------------------------------------------------
# Complet Plot (eOm = f(s1, s2))
#------------------------------------------------


#Multiplot for all orders
#---------------------
ttm_l   = list();
ppl     = list();
for(i in 1:length(vorders))
{
  order = vorders[i];
  ttm_l[[i]]  = ttm_all[which(ttm_all$order == order),]      #select order
  ppl[[i]]    = plotdf_tile_n(ttm_l[[i]], si, sj, sil, sjl, colorCol, colorLabel, colorMidPoint, isLegendOn, "right")
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
  ggsave(ppl[[1]], file=paste0(currentfolder, "R_eIm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, ".eps"))
}else{
  #plot the legend only
  ggsave(pLegend, file=paste0(currentfolder, "R_eIm_legend_", Type, ".eps"))
}


########################################################################################################################################
stop()
########################################################################################################################################

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

#Multiplot for all orders
#---------------------
ttm_l   = list();
ppl     = list();
for(i in 1:length(vorders))
{
  #Plot
  order = vorders[i];
  ttm_l[[i]]  = ttm_all[which(ttm_all$order == order),]      #select order
  
  #------------------------------------------------
  # Selecting specific lines
  #------------------------------------------------
  ttm_s1= ttm_l[[i]][which(ttm_l[[i]]$s3 == 0),]
  ttm_s3= ttm_l[[i]][which(ttm_l[[i]]$s1 == 0),]
  ttm_s1s3= ttm_l[[i]][which(ttm_l[[i]]$s1 == ttm_l[[i]]$s3),]
  ttm_s1ms3= ttm_l[[i]][which(ttm_l[[i]]$s1 == -ttm_l[[i]]$s3),]
  
  #------------------------------------------------
  # Plot
  #------------------------------------------------
  ppl[[i]]    = plotdf_point(ttm_l[[i]], xi, yi, xist, yist, "log10eOm", "log10eOm", 0)
  ppl[[i]]    = ppl[[i]]+scale_colour_gradient2(colorLabel, space="Lab", midpoint = -8, mid = "white", high = muted("blue"), guide = FALSE)
  
  #Add s1/s3 axes
  ppl[[i]] = ppl[[i]] + geom_point(data = ttm_s1,     aes_string(x = xi, y = yi), color = "black", shape = 16, size = 2)
  ppl[[i]] = ppl[[i]] + geom_point(data = ttm_s3,     aes_string(x = xi, y = yi), color = "black", shape = 16, size = 2)
  ppl[[i]] = ppl[[i]] + geom_point(data = ttm_s1s3,   aes_string(x = xi, y = yi), color = "brown", shape = 16, size = 2)
  ppl[[i]] = ppl[[i]] + geom_point(data = ttm_s1ms3,  aes_string(x = xi, y = yi), color = "brown", shape = 16, size = 2)
  
  #Add the moon
  if(xi == "xPH")
  {
    primaryR    = 1737.10/abs(primaryFactor)  #Moon's radius
    primaryPos  = gamma*L/primaryFactor       #Moon position wrt to Li
    ppl[[i]]    = addMoon(ppl[[i]], x = primaryPos+5000, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2)
    ppl[[i]]    = ppl[[i]]+ annotate("text", x = sign(primaryFactor)*gamma*L/primaryFactor, y = 5000, label = "Moon", size = 6)
    ppl[[i]]    = ppl[[i]]+ coord_fixed(ratio = 1)
  }
  
  #Limits
  #ppl[[i]]    = ppl[[i]]+ scale_y_continuous(limits = c(-55000, 55000), breaks=c(-30000, 0, 30000))
}
#Actual plot in multiplot format
pMult_ph = multiplot(plotlist = ppl, cols = 2)


#Save to eps file
ggsave(ppl[[1]], file=paste0(currentfolder, "R_eIm_", "ofs_" , ofs_order, "_", vorders[1], "_", Type, "_withMoon_", xi, yi,".eps"))
