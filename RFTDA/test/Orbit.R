# Script to plot an orbit along with its associated precisions
# in the parameterization method of the QBFBP/RTBP around L1/L2 of the Earth-Moon system
#----------------------------------------------------------------------------------------

# Load libraries
#------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
#Libration point
Li = "L2"
MODEL = "QBFBP"
size = "50_planar"
ltype = ifelse(MODEL=="QBFBP", "solid", "dashed")
Period = ifelse(MODEL=="QBFBP", 6.79119387190792, 2*pi)
#Constants
tsize = 20
lsize = 20
psize = 2
#Theme
custom_theme = theme(
  axis.text.x  = element_text(colour="grey20", size=lsize,angle=0,hjust=.5,vjust=.5),
  axis.text.y  = element_text(colour="grey20", size=lsize),
  axis.title.x = element_text(colour="grey20", size=tsize, vjust = 0.5),
  axis.title.y = element_text(colour="grey20", size=tsize, vjust = 0.5),
  legend.text  = element_text(colour="grey20", size=lsize),
  legend.title = element_text(colour="grey20", size=tsize),
  legend.title.align = 0.5,
  legend.justification=c(1,0), legend.position=c(1,0),
  legend.key.size  = unit(1, 'cm'),
  legend.key.width = unit(2, 'cm')
)

blank_theme = theme(
  axis.text.x  = element_blank(),
  axis.text.y  = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.text  = element_blank(),
  legend.title = element_blank(),
  plot.title = element_text(colour="grey20", size=tsize),
  plot.background=element_blank(),
  axis.ticks=element_blank()
)

#Index order
#if(MODEL == "QBFBP") dfindex = c(2,5,10,18) else dfindex = c(2,5,10,20)
#dfnumber = 1:4
if(MODEL == "QBFBP") dfindex = c(5,10,18) else dfindex = c(5,10,20)
dfnumber = 1:3

#--------------------------------------------------------------------------------------------------------------------------
#                                             eH
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "eH_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0("~/BackUpBox/PhD/OOFTDA/plot/",MODEL,"/",Li,"/orbits/",dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("x", "y")
  dflist[[i]]$order = dfindex[i]
}
# Concatenate the results
eH = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){eH = rbind(eH, dflist[[i]])}

#Plot
#------------------------------------------------
pH = ggplot()
#Data
pH = pH + geom_line(data = eH, aes(x, y, colour = factor(order)), size = psize, linetype = ltype)
#Scaling
if(MODEL == "QBFBP") pH = pH + scale_y_log10(limits=c(1e-11, 1e3), breaks = 10^(-8:1*2)) else 
  pH = pH + scale_y_log10(limits=c(1e-8, 1e3), breaks = 10^(-8:1*2))
pH = pH + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), 
                             labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"))  #cont scale on the x axis 
#Theme
#pH = pH + theme_bw()
pH = pH + custom_theme
#Labels
# pH = pH + labs(x = expression( ~ '\n' ~ abs(W(s,0))), y = "Error on the potential [-]\n")
pH = pH + labs(x = "\n time (%T)", y = "eH [-] \n")
pH = pH + scale_colour_discrete(name="Order")
#Display the plot
pH


spH = ggplot()
#Data
spH = spH + geom_line(data = eH, aes(x, y, colour = factor(order)), size = psize, linetype = ltype)
spH = spH + scale_colour_discrete(guide = F)
#Scaling
if(MODEL == "QBFBP") spH = spH + scale_y_log10(limits=c(5e-4, 5e-2), breaks = 10^(-8:1)) else 
  spH = spH + scale_y_log10(limits=c(1e-8, 1e3), breaks = 10^(-8:1*2))
spH = spH + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), 
                             labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"), limits = c(Period-0.03, Period+0.01))  #cont scale on the x axis 
#Theme
spH = spH + theme_bw()
spH = spH + custom_theme
spH = spH + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.text  = element_blank(),
  legend.title = element_blank(),
  plot.title = element_text(colour="grey20", size=tsize),
  plot.background=element_blank(),
  axis.ticks=element_blank()
)

#Labels
# spH = spH + labs(x = expression( ~ '\n' ~ abs(W(s,0))), y = "Error on the potential [-]\n")
spH = spH + labs(x = "\n time (%T)", y = "eH [-] \n", title = "Zoom at t = T")
spH = spH + scale_colour_discrete(name="Order")
#Display the plot
spH

#Add a viewport and print
vp <- viewport(width = 0.35, height = 0.35, x = 0.7, y = 0.7, just = c("right","top"))
print(pH)
print(spH, vp = vp)



#--------------------------------------------------------------------------------------------------------------------------
#                                             eI
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "eI_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0("~/BackUpBox/PhD/OOFTDA/plot/",MODEL,"/",Li,"/orbits/",dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("x", "y")
  dflist[[i]]$order = dfindex[i]
}
# Concatenate the results
eI = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){eI = rbind(eI, dflist[[i]])}

#Plot
#------------------------------------------------
pI = ggplot()
#Data
pI = pI + geom_line(data = eI, aes(x, y, colour = factor(order)), size = psize, linetype = ltype)
#Scaling
# if(MODEL == "QBFBP") pI = pI + scale_y_log10(limits=c(1e-6, 1e0), breaks = 10^(-8:1*2)) else 
#   pI = pI + scale_y_log10(limits=c(1e-6, 1e0), breaks = 10^(-8:1*2))
#---
# if(MODEL == "QBFBP") pI = pI + scale_y_log10(limits=c(1e-6, 1e-2), breaks = 10^(-8:1*2)) else 
#   pI = pI + scale_y_log10(limits=c(1e-6, 1e0), breaks = 10^(-8:1*2))
#---
if(MODEL == "QBFBP") pI = pI + scale_y_log10(limits=c(1e-8, 1e-6), breaks = 10^(-8:1)) else 
  pI = pI + scale_y_log10(limits=c(1e-6, 1e0), breaks = 10^(-8:1*2))
#---
pI = pI + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), 
                             labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"))  #cont scale on the x axis 

#Theme
#pI = pI + theme_bw()
pI = pI + custom_theme
#Labels
# pI = pI + labs(x = expression( ~ '\n' ~ abs(W(s,0))), y = "Error on the potential [-]\n")
pI = pI + labs(x = "\n time (%T)", y = "eI [-] \n")
pI = pI + scale_colour_discrete(name="Order")
#Display the plot
pI


#--------------------------------------------------------------------------------------------------------------------------
#                                             eO
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "eO_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0("~/BackUpBox/PhD/OOFTDA/plot/",MODEL,"/",Li,"/orbits/",dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("x", "y")
  dflist[[i]]$order = dfindex[i]
}
# Concatenate the results
eO = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){eO = rbind(eO, dflist[[i]])}

#Plot
#------------------------------------------------
pO = ggplot()
#Data
pO = pO + geom_line(data = eO, aes(x, y, colour = factor(order)), size = psize, linetype = ltype)
#Scaling
# if(MODEL == "QBFBP") pO = pO + scale_y_log10(limits=c(1e-7, 1e2), breaks = 10^(-8:1*2)) else 
#   pO = pO + scale_y_log10(limits=c(1e-7, 1e2), breaks = 10^(-8:1*2))

if(MODEL == "QBFBP") pO = pO + scale_y_log10(limits=c(1e-9, 1e-1), breaks = 10^(-8:1*2)) else 
  pO = pO + scale_y_log10(limits=c(1e-7, 1e2), breaks = 10^(-8:1*2))

pO = pO + scale_x_continuous(breaks = seq(0.0, Period, 0.25*Period), 
                             labels = c("0.0", "0.25 T", "0.5 T", "0.75 T", "T"))  #cont scale on the x axis 
#Theme
#pO = pO + theme_bw()
pO = pO + custom_theme
#Labels
# pO = pO + labs(x = expression( ~ '\n' ~ abs(W(s,0))), y = "Error on the potential [-]\n")
pO = pO + labs(x = "\n time (%T)", y = "eO [-] \n")
pO = pO + scale_colour_discrete(name="Order")
#Display the plot
pO


#--------------------------------------------------------------------------------------------------------------------------
#                                             XY
#--------------------------------------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------
dflist = list();
dftype = "XYZ_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0("~/BackUpBox/PhD/OOFTDA/plot/",MODEL,"/",Li,"/orbits/",dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("t", "x", "y", "z")
  dflist[[i]]$order = dfindex[i]
}

# Concatenate the results
orbit = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){orbit = rbind(orbit, dflist[[i]])}

# Load data (PM)
#------------------------------------------------
dflist = list();
dftype = "XYZ_PM_Order_"
for(i in dfnumber)
{
  dflist[[i]] = read.table(paste0("~/BackUpBox/PhD/OOFTDA/plot/",MODEL,"/",Li,"/orbits/",dftype,dfindex[i],"_Size_", size, ".txt"), header = F)
  colnames(dflist[[i]]) = c("t", "x", "y", "z")
  dflist[[i]]$order = dfindex[i]
}

# Concatenate the results
orbit_pm = rbind(dflist[[1]], dflist[[2]])
for(i in 3:length(dfnumber)){orbit_pm = rbind(orbit_pm, dflist[[i]])}

#Select half period time
orbit_half = ddply(orbit, ~order, function(x){x[which.min(abs(x$t-0.25*Period)),]})
#Select when orbital precision is > maxprec
maxprec = 1e-6
eO_prec = ddply(eO, ~order, function(x){x[which.min(abs(x$y-maxprec)),]})  #time for which eO ~ maxprec
orbit_prec = ddply(orbit, ~order, function(orb){orb[which.min(abs(orb$t-eO_prec$x[which(eO_prec$order == 18)])),]}) #select order 18

# Plot
#------------------------------------------------
porb = ggplot()
#Data
# Different orders
porb = porb + geom_path(data = orbit, aes(x, y, colour = factor(order)), size = psize, linetype = ltype)
porb = porb + scale_colour_grey(start = 0.9, end = 0.0, name="Order", guide = F)
#Center manifold
porb = porb + geom_path(data  = orbit_pm[which(orbit_pm$order==18),], aes(x, y, color = "PM"), size = psize, linetype = "dashed", colour = "limegreen")
# Starting point
orbit18 = orbit[which(orbit$order == 18),]
porb = porb + geom_point(data = orbit18[which(orbit18$t == 0.0),], aes(x, y), size = 4, colour = "black", fill = "white", pch = 21)
# Point at a given precision
porb = porb + geom_point(data = orbit_prec[which(orbit_prec$order==18),], aes(x, y), size = 4, pch = 22, colour = "black", fill = "white", solid = T)
#Scaling
#porb = porb + scale_x_continuous(limits = c(-0.35, 0.3))  #cont scale on the x axis 
#porb = porb + scale_y_continuous(limits = c(-0.3, 0.4))  #cont scale on the y axis
#---
#porb = porb + scale_x_continuous(limits = c(-0.05, 0.1))  #cont scale on the x axis 
#porb = porb + scale_y_continuous(limits = c(-0.15, 0.15))  #cont scale on the y axis 
#---
# porb = porb + scale_x_continuous(limits = c(-0.5, 1))  #cont scale on the x axis 
# porb = porb + scale_y_continuous(limits = c(-1, 1))  #cont scale on the y axis 
#---
porb = porb + scale_x_continuous(limits = c(-1e-3, 0.002))  #cont scale on the x axis 
porb = porb + scale_y_continuous(limits = c(-0.0025, 0.0025))  #cont scale on the y axis 

#Theme
porb = porb + theme_bw()
porb = porb + custom_theme
#Labels
porb = porb + labs(x = "\nx (NC)", y = "y (NC)\n")
#Annotation
#porb = porb + annotate("text", x = -0.3, y = -0.21, label = "CM \n Order 18", parse = F, size = 8, colour = "limegreen", fontface = "bold")
#porb = porb + geom_segment(aes(x = -0.25, y = -0.2, xend = -0.19, yend = -0.2), arrow = arrow(length = unit(0.5, "cm")), size = psize)
#Display the plot
porb


# Plot
#------------------------------------------------
sporb = ggplot()
#Data
# Different orders
sporb = sporb + geom_path(data = orbit, aes(x, y, colour = factor(order)), size = 1, linetype = ltype)
sporb = sporb + scale_colour_grey(start = 0.9, end = 0.0, name="Order", guide = F)
#Scaling
#Theme
sporb = sporb + theme_bw()
sporb = sporb + blank_theme
#Labels
sporb = sporb + labs(x = "\nx (NC)", y = "y (NC)\n")
#Title
#sporb = sporb + annotate("text", x = -10, y = -6, label = "Full data", parse = F, size = 7, colour="grey20")
#sporb = sporb + annotate("text", x = 0.15, y = -0.5, label = "Full data", parse = F, size = 7, colour="grey20")
#sporb = sporb + annotate("text", x = -5, y = -4, label = "Full data", parse = F, size = 7, colour="grey20")
sporb = sporb + annotate("text", x = 0.006, y = -0.001, label = "Full data", parse = F, size = 7, colour="grey20")
#Add a viewport and print
vp <- viewport(width = 0.35, height = 0.35, x = 0.96, y = 0.96, just = c("right","top"))
print(porb)
print(sporb, vp = vp)


# Plot
#------------------------------------------------
pyz = ggplot()
#Data
# Different orders
pyz = pyz + geom_path(data = orbit, aes(y, z, colour = factor(order)), size = psize, linetype = ltype)
pyz = pyz + scale_colour_grey(start = 0.9, end = 0.0, name="Order")
#Center manifold
pyz = pyz + geom_path(data  = orbit_pm[which(orbit_pm$order==18),], aes(y, z, color = "PM"), size = psize, linetype = "dashed", colour = "limegreen")
# Starting point
orbit18 = orbit[which(orbit$order == 18),]
pyz = pyz + geom_point(data = orbit18[which(orbit18$t == 0.0),], aes(y, z), size = 4, colour = "black", fill = "white", pch = 21)
# Point at a given precision
pyz = pyz + geom_point(data = orbit_prec[which(orbit_prec$order==18),], aes(y, z), size = 4, pch = 22, colour = "black", fill = "white", solid = T)
#Scaling
pyz = pyz + scale_x_continuous(limits = c(-0.5, 0.3))  #cont scale on the y axis 
pyz = pyz + scale_y_continuous(limits = c(-0.7, 0.7))  #cont scale on the z axis 
#Theme
pyz = pyz + theme_bw()
pyz = pyz + custom_theme
pyz = pyz + theme(legend.justification=c(1,1), legend.position=c(1,1))
#Labels
pyz = pyz + labs(x = "\ny (NC)", y = "z (NC)\n")
#Annotation
#pyz = pyz + annotate("text", x = -0.3, y = 0.22, label = "CM \n Order 18", parse = F, size = 8, colour = "limegreen", fontface = "bold")
#pyz = pyz + geom_segment(aes(x = -0.25, y = 0.22, xend = -0.2, yend = 0.5), arrow = arrow(length = unit(0.5, "cm")), size = psize)
#Display the plot
pyz


# Plot
#------------------------------------------------
spyz = ggplot()
#Data
# Different orders
spyz = spyz + geom_path(data = orbit, aes(y, z, colour = factor(order)), size = 1, linetype = ltype)
spyz = spyz + scale_colour_grey(start = 0.9, end = 0.0, name="Order", guide = F)
#Scaling
#Theme
spyz = spyz + theme_bw()
spyz = spyz + blank_theme
#Labels
spyz = spyz + labs(x = "\ny (NC)", y = "z (NC)\n")
#Title
spyz = spyz + annotate("text", x = -10, y = -0.5, label = "Full data", parse = F, size = 7, colour="grey20")
#spyz = spyz + annotate("text", x = 0.006, y = -0.001, label = "Full data", parse = F, size = 7, colour="grey20")
#Add a viewport and print
vp <- viewport(width = 0.35, height = 0.35, x = 0.12, y = 0.12, just = c("left","bottom"))
print(pyz)
print(spyz, vp = vp)

#--------------------------------------------------------------------------------------------------------------------------
#Save in eps file
#--------------------------------------------------------------------------------------------------------------------------
ggsave(pH, file = paste0("~/BackUpBox/PhD/OOFTDA/plot/",MODEL,"/",Li,"/orbits/R_eH_Size_", size, ".eps"))
ggsave(pI, file = paste0("~/BackUpBox/PhD/OOFTDA/plot/",MODEL,"/",Li,"/orbits/R_eI_Size_", size, ".eps"))
ggsave(pO, file = paste0("~/BackUpBox/PhD/OOFTDA/plot/",MODEL,"/",Li,"/orbits/R_eO_Size_", size, ".eps"))