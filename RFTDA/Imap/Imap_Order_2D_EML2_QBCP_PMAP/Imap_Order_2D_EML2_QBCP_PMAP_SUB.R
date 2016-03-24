#------------------------------------------------
# Import data (just once)
#------------------------------------------------
source('Imap/Imap_Order_2D_EML2_QBCP/Imap_Order_2D_EML2_QBCP_HEADER.R')

#------------------------------------------------
# Orders to Subs
#------------------------------------------------
ordersToSub = c(25, 30)

#------------------------------------------------
#Additionnal parameters
#------------------------------------------------
isLegendOn = 0;
legendOnly = 0;
legendOnTop = 0;

#------------------------------------------------
# Get the data in proper form for subtraction
#------------------------------------------------
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
#Plot in z space
#--------------------
xi = "xEM"
yi = "yEM"
xist = "x [km]"
yist = "y [km]"
primaryFactor = 1.0#-gamma*Le

#Plot the difference in precision
if(isLegendOn){
  pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "log10eOm", "log10eOm", 0)+scale_colour_gradient2(legTitle, space="Lab", midpoint = 0, mid = "white", high = muted("blue"))
  
}else{
  pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "log10eOm", "log10eOm", 0)+scale_colour_gradient2(legTitle, space="Lab", midpoint = 0, mid = "white", high = muted("blue"), guide = FALSE)
  
}

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
# pph = pph + geom_point(data = ttm_s1,  aes_string(x = xi, y = yi), color = "black", shape = 16, size = 2)
# pph = pph + geom_point(data = ttm_s3,   aes_string(x = xi, y = yi), color = "black", shape = 16, size = 2)
# pph = pph + geom_point(data = ttm_s1s3,  aes_string(x = xi, y = yi), color = "brown", shape = 16, size = 2)
# pph = pph + geom_point(data = ttm_s1ms3,  aes_string(x = xi, y = yi), color = "brown", shape = 16, size = 2)
# #Annotate
# lsize = 8
# pph    = pph+ annotate("text", x = 22000, y = 52000, label = "s[1] == s[3]", size = lsize, parse=TRUE)
# pph    = pph+ annotate("text", x = 28000, y = 12000, label = "s[1] == -s[3]", size = lsize, parse=TRUE)
# pph    = pph+ annotate("text", x = 22000, y = -24000, label = "s[1] == 0", size = lsize, parse=TRUE)
# pph    = pph+ annotate("text", x = -5000, y = -34000, label = "s[3] == 0", size = lsize, parse=TRUE)

#Plot the sign of the difference
#pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "sign", "sign", 1)

#Add the moon
# if(xi == "xPH")
# {
#   primaryR    = 1737.10/abs(primaryFactor)  #Moon's radius
#   primaryPos  = gamma*L/primaryFactor       #Moon position wrt to Li
#   pph    = addMoon(pph, x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2)
#   pph    = pph+ annotate("text", x = sign(primaryFactor)*gamma*L/primaryFactor, y = 5000, label = "Moon", size = 6)
#   pph    = pph+ coord_fixed(ratio = 1)
# }

#Add the line
#pph = pph + geom_point(data = ttm_sub_line,  aes_string(x = xi, y = yi), color = "black")

#Limits
#pph    = pph+ scale_y_continuous(limits = c(-55000, 55000), breaks=c(-30000, 0, 30000))
pph    = pph+ scale_x_continuous(breaks=c(-0.02, 0, 0.03))

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
ggsave(pph, file=paste0(currentfolder, "R_eIm_sub_", ordersToSub[1], "_", ordersToSub[2], "_ofs_" , ofs_order, "_", vorders[1], "_", Type, "_", xi, yi,"_withlinesannoted.eps"))
