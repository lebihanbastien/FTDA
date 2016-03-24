############################################################
# R script to plot data for a precision map of the QBCP:
# Difference between two orders given in the array ordersToSub
#
# 01/2016
############################################################

#------------------------------------------------------------------------------------
# Building the subtraction
#------------------------------------------------------------------------------------
ttm_sub = data.frame()
ttm_l = list();

#------------------------------------------------
# Loop
#------------------------------------------------
for (i in c(1,2))  #loop on the orders
{
  #Select proper order
  ttm_l[[i]] = ttm_all[which(ttm_all$order == ordersToSub[i]),]
  #Sort
  ttm_l[[i]] = ttm_l[[i]][order(ttm_l[[i]]$s1, ttm_l[[i]]$s2, ttm_l[[i]]$s3, ttm_l[[i]]$s4),]
}


#------------------------------------------------
# Perform subtraction
#------------------------------------------------
ttm_sub = data.frame(s1=ttm_l[[1]]$s1, 
                     s2=ttm_l[[1]]$s2, 
                     s3=ttm_l[[1]]$s3, 
                     s4=ttm_l[[1]]$s4,
                     xPH = ttm_l[[1]]$xPH, 
                     yPH = ttm_l[[1]]$yPH, 
                     zPH = ttm_l[[1]]$zPH, 
                     xEM = ttm_l[[1]]$xEM, 
                     yEM = ttm_l[[1]]$yEM, 
                     zEM = ttm_l[[1]]$zEM,
                     log10eOm=ttm_l[[1]]$log10eOm-ttm_l[[2]]$log10eOm)

ttm_sub$sign = ttm_sub$log10eOm > -0

#------------------------------------------------
# Selecting one line
#------------------------------------------------
#ttm_sub_line = ttm_sub[which(ttm_sub$s3 == ttm_sub$s1),]
#ttm_sub_line = ttm_sub[which(ttm_sub$s4 == 0),]

#------------------------------------------------
# Plot (Tile)
#------------------------------------------------
legTitle = expression(Delta~log[10]~bgroup("(",e[I]~bgroup("(",0,")"),")") ~~~~ "")

#--------------------
#Plot in s space
#--------------------
ppl = plotdf_tile(ttm_sub , si, sj, isLegendOn, 0)

#--------------------
#Plot in ohysical space
#--------------------
xi = "xEM"
yi = "yEM"
xist = "x [-]"
yist = "y [-]"
primaryFactor = 1.0#-gamma*Le

#Plot the difference in precision
if(isLegendOn){
  pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "log10eOm", "log10eOm", 0)+scale_colour_gradient2(legTitle, space="Lab", midpoint = 0, mid = "white", high = muted("blue"))
  
}else{
  pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "log10eOm", "log10eOm", 0)+scale_colour_gradient2(legTitle, space="Lab", midpoint = 0, mid = "white", high = muted("blue"), guide = FALSE)
  
}
#Plot the sign of the difference
#pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "sign", "sign", 1)

#------------------------------------------------
# Selecting specific lines
#------------------------------------------------
ttm_s1= ttm_sub[which(ttm_sub$s3 == 0),]
ttm_s3= ttm_sub[which(ttm_sub$s1 == 0),]
ttm_s1s3= ttm_sub[which(ttm_sub$s1 == ttm_sub$s3),]
ttm_s1ms3= ttm_sub[which(ttm_sub$s1 == -ttm_sub$s3),]

#------------------------------------------------
# Plot
#------------------------------------------------
# #Add s1/s3 axes
if(pph_withAxes)
{
    pph = pph + geom_point(data = ttm_s1,  aes_string(x = xi, y = yi), color = "red", shape = 16, size = 2)
    pph = pph + geom_point(data = ttm_s3,   aes_string(x = xi, y = yi), color = "blue", shape = 16, size = 2)
    pph = pph + geom_point(data = ttm_s1s3,  aes_string(x = xi, y = yi), color = "green", shape = 16, size = 2)
    pph = pph + geom_point(data = ttm_s1ms3,  aes_string(x = xi, y = yi), color = "black", shape = 16, size = 2)
    if(!pph_withNotes) pph    = pph+ scale_y_continuous(breaks = c(-0.1, -0.05, 0.00, 0.05, 0.1))
}
# #Annotate
lsize = 8

if(pph_withNotes)
{
  if(Li == "L2")
  {
    #EML2
    pph    = pph + annotate("text", x = -1.11,  y =  0.14, label = "s[1] == s[3]", size = lsize, parse=TRUE)
    pph    = pph + annotate("text", x = -1.10,  y =  0.015, label = "s[1] == -s[3]", size = lsize, parse=TRUE)
    pph    = pph + annotate("text", x = -1.17,  y = +0.08, label = "s[1] == 0", size = lsize, parse=TRUE)
    pph    = pph + annotate("text", x = -1.17,  y = -0.08, label = "s[3] == 0", size = lsize, parse=TRUE)
  }else
  {
    #EML1
    pph    = pph + annotate("text", x = -0.888,  y =  0.14, label = "s[1] == s[3]", size = lsize, parse=TRUE)
    pph    = pph + annotate("text", x = -0.895,  y =  0.015, label = "s[1] == -s[3]", size = lsize, parse=TRUE)
    pph    = pph + annotate("text", x = -0.825,  y = -0.09, label = "s[1] == 0", size = lsize, parse=TRUE)
    pph    = pph + annotate("text", x = -0.825,  y = +0.09, label = "s[3] == 0", size = lsize, parse=TRUE)
  }
}

#Add the moon
if(pph_withMoon)
{
  if(xi == "xPH")
  {
    primaryR    = 1737.10/abs(primaryFactor)  #Moon's radius
    primaryPos  = gamma*L/primaryFactor       #Moon position wrt to Li
    pph    = addMoon(pph, x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2)
    pph    = pph+ annotate("text", x = sign(primaryFactor)*gamma*L/primaryFactor, y = 5000, label = "Moon", size = 6)
    pph    = pph+ coord_fixed(ratio = 1)
  }
}


#Limits
if(pph_limits)
{
  pph    = pph+ scale_x_continuous(limits = pph_x_lim)
  pph    = pph+ scale_y_continuous(limits = pph_y_lim)
}

#Legend only (don't forget to display the legend in plotfd_tile!)
#--------------------
if(legendOnly == 1)
{
  if(legendOnTop == 1)
  {
    pph = pph+theme(legend.position="top", legend.title.align= 0.5)
  }
  dev.off()                    #clear all plots
  legend <- g_legend(pph) #select the legend
  pLegend = grid.draw(legend)  #displays it
}

#Save to eps file
if(pph_withAxes && pph_withNotes){
  ggsave(pph, file=paste0(currentfolder, "R_eIm_sub_", ordersToSub[1], "_", ordersToSub[2], "_ofs_" , ofs_order, "_", vorders[1], "_", Type, "_", xi, yi,"_withlinesannoted.eps"))

}else if(pph_withAxes){
  ggsave(pph, file=paste0(currentfolder, "R_eIm_sub_", ordersToSub[1], "_", ordersToSub[2], "_ofs_" , ofs_order, "_", vorders[1], "_", Type, "_", xi, yi,"_withlines.eps"))
}else
{
  ggsave(pph, file=paste0(currentfolder, "R_eIm_sub_", ordersToSub[1], "_", ordersToSub[2], "_ofs_" , ofs_order, "_", vorders[1], "_", Type, "_", xi, yi,".eps"))
}