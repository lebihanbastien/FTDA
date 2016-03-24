#------------------------------------------------
# Import data (just once)
#------------------------------------------------
#source('Imap/Imap_Order_4D_EML2_QBCP/Imap_Order_4D_EML2_QBCP_HEADER.R')

#------------------------------------------------
# Select the plots
#------------------------------------------------
si = "s1"
sj = "s2"

#------------------------------------------------
# Orders to Subs
#------------------------------------------------
ordersToSub = c(25, 30)

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
ind = 8
#Cut
ttm_cut = imap_4d_cut(si, sj, ind)


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
  ttm_l[[i]] = ttm_cut[which(ttm_cut$order == ordersToSub[i]),]
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
                     log10eOm=ttm_l[[1]]$log10eOm-ttm_l[[2]]$log10eOm)

ttm_sub$sign = ttm_sub$log10eOm > -0

#------------------------------------------------
# Selecting one line
#------------------------------------------------
ttm_sub_line = ttm_sub[which(ttm_sub$s1 == -ttm_sub$s3),]

#------------------------------------------------
# Plot (Tile)
#------------------------------------------------

#--------------------
#Plot in s space
#--------------------
ppl = plotdf_tile(ttm_sub , si, sj, isLegendOn, 0)

#--------------------
#Plot in z space
#--------------------
xi = "xPH"
yi = "yPH"
xist = "x [km]"
yist = "y [km]"
primaryFactor = 1.0#-gamma*Le

#Plot the difference in precision
pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "log10eOm", "log10eOm", 0, pointSize = 4)+scale_colour_gradient2("log10eIm", space="Lab", midpoint = 0, mid = "white", high = muted("blue"))

#Plot the sign of the difference
#pph = plotdf_point(ttm_sub, xi, yi, xist, yist, "sign", "sign", 1)

if(xi == "xPH")
{
  #Earth direction
  pph    = pph+ annotate("text", x = 57000, y = -25000, label = "to Earth \ ", size = 6)
  pph    = pph+ geom_segment(aes(x = 50000, y = -30000, xend = 60000, yend = -30000), 
                             colour = "black", 
                             arrow = arrow(length = unit(0.3, "cm"), 
                                           type = "closed"))
  #Add the moon
  primaryR    = 1737.10/abs(primaryFactor)  #Moon's radius
  primaryPos  = gamma*L/primaryFactor       #Moon position wrt to Li
  pph    = addMoon(pph, x = primaryPos, y = 0, primaryR, surfSize = 0.4, cratSize = 0.2)
  pph    = pph+ annotate("text", x = sign(primaryFactor)*gamma*L/primaryFactor, y = 5000, label = "Moon", size = 6)
}

#Add the line
pph = pph + geom_point(data = ttm_sub_line,  aes(x = xPH, y = yPH), color = "black")
