# R script to handle a precision map of the QBFBP around EML1,2
#---------------------------------------------------------------------------------------------------------------------------
# Constants to set
#------------------------------------------------
Li = "L3"      #Libration point
Type = "planar0" #selection or global


# Load libraries
#------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)

# Constants
#------------------------------------------------
psize = 1
lsize = 1
bsize = 7
ssize = 0.1

# Plot function
#------------------------------------------------
#Tile plot
plotdf_tile<- function(ttm_c, st1, st2) {
  
  #Build the -log10(eOm) if it does not exist
  if(is.null(ttm_c$log10eOm))
  {
    ttm_c$log10eOm = -log10(ttm_c$eOm)
  }
  #Ggplot init
  pt = ggplot()
  # Plot tiles
  pt = pt + geom_tile(data = ttm_c, aes_string(st1, st2, fill = "log10eOm"))
  #Get the proper color gradient
  pt = pt + scale_fill_gradient2(expression(-log[10](e[O])), space="Lab", midpoint = min(ttm_c$log10eOm)-1, mid = "white", high = muted("blue")) 
  #Plot all the tiles for which the precision is < 1e-6
  #pt = pt + geom_tile(data = ttm6_all, aes(s1, s3),  color = "red", alpha = 0, size = ssize)
  #Labels
  pt = pt + labs(x = st1, y = st2)
  #Theme
  #pt = pt + theme_bw()
  pt = pt + custom_theme
  
  return(pt)
}

#Point plot
plotdf_point<- function(ttm_c, st1, st2) {
  
  #Build the -log10(eOm) if it does not exist
  if(is.null(ttm_c$log10eOm))
  {
    ttm_c$log10eOm = -log10(ttm_c$eOm)
  }
  #Ggplot init
  px = ggplot()
  #Full dataset
  px = px + geom_point(data = ttm_c, aes_string(st1, st2, colour= "log10eOm"), size = bsize)
  #Precision < 10^{-6}
  #px = px + geom_point(data = ttm6_all, aes(x, y, colour = -log10(eOm)), size = bsize)
  #Get the proper color gradient
  px = px + scale_color_gradient2(expression(-log[10](e[O])), space="Lab", midpoint = min(ttm_c$log10eOm)-1, mid = "white", high = muted("blue")) 
  #Labels
  px = px + labs(x = st1, y = st2)
  #Theme
  px= px + custom_theme
  #Display
  px
  
  return(px)
}

plotdf_line<- function(ttm, st1, st2) {
  
  #Ggplot init
  px = ggplot()
  #Changing temporarily the names of the desired columns
  i1 = which(colnames(ttm) == st1)
  i2 = which(colnames(ttm) == st2)
  colnames(ttm)[i1] = "temp1"
  colnames(ttm)[i2] = "temp2"
  
  #Plot
  px = px + geom_line(data = ttm, aes(x = temp1, y = temp2, color = factor(type), linetype = factor(type)), size = 1.2)
  #px = px + geom_line(data = ttm, aes(x = temp1, y = temp2, color = factor(type)), size = 1.2)
  
  #Colour scheme
  px = px + scale_color_discrete(name="Order")
  #Labels
  #px = px + labs(x = st1, y = st2)
  px = px + labs(x = "Distance from L2 [km]", y = "-log(orbital error) after T/2")
  #Theme
  px= px + custom_theme
  
  #Display
  px
  
  return(px)
}

plotdf_pointf<- function(ttm, st1, st2) {
  
  #Ggplot init
  px = ggplot()
  #Changing temporarily the names of the desired columns
  i1 = which(colnames(ttm) == st1)
  i2 = which(colnames(ttm) == st2)
  colnames(ttm)[i1] = "temp1"
  colnames(ttm)[i2] = "temp2"
  
  #Plot
  px = px + geom_point(data = ttm, aes(x = temp1, y = temp2, color = factor(order)))
  
  
  #Colour scheme
  px = px + scale_color_discrete(name="Order")
  #Labels
  px = px + labs(x = st1, y = st2)
  #Theme
  px= px + custom_theme
  
  #Display
  px
  
  return(px)
}



#Libration point
L1 = "L1"
L2 = "L2"
L3 = "L3"

#Normalized units (gamma, c1)
if(Li == L1)
{
  gamma = 0.150934272990064
  c1    = -5.5448979798087
}
if(Li == L2)
{
  gamma = 0.16783273173707
  c1    = -6.8859163415402
}
if(Li == L3)
{
  gamma = +9.929120655179929e-01
  c1    = +1.012237318938304e+00
}


#Theme
tsize = 20  #title size
lsize = 20  #legend size

custom_theme = theme(
  axis.text.x  = element_text(colour="grey20", size=lsize,angle=0,hjust=.5,vjust=.5),
  axis.text.y  = element_text(colour="grey20", size=lsize),
  axis.title.x = element_text(colour="grey20", size=tsize, vjust = 0.5),
  axis.title.y = element_text(colour="grey20", size=tsize, vjust = 0.5),
  legend.text  = element_text(colour="grey20", size=lsize),
  legend.title = element_text(colour="grey20", size=tsize),
  legend.title.align = 0.5,
  legend.justification=c(1,0), #legend.position=c(1,0),
  legend.key.size  = unit(1, 'cm'),
  legend.key.width = unit(2, 'cm')
)

# R options
#------------------------------------------------
options(digits = 15)


ttm_all = data.frame()
i = 5

  # Filename to check
  #------------------------------------------------
  fileprefix = paste0("~/BackUpBox/PhD/OOFTDA/fprint/QBCP/", Li, "/eOm_", Type, "_ofs_30_order_", toString(i))
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
  
  # -log10(precision)
  ttm_c$log10eOm = -log10(ttm_c$eOm)
  
  #Type of COC computation
  ttm_c$type = "Old"
  
  #rbind in ttm_all
  ttm_all = rbind(ttm_all, ttm_c)
  #------------------------------------------------
  
  
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0("~/BackUpBox/PhD/OOFTDA/fprint/QBCP/", Li, "/eOm_", Type, "_ofs_30_order_", toString(i), "_test")
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
  
  # -log10(precision)
  ttm_c$log10eOm = -log10(ttm_c$eOm)
  
  #Type of COC computation
  ttm_c$type = "New"
  
  #rbind in ttm_all
  ttm_all = rbind(ttm_all, ttm_c)
  #------------------------------------------------


# From NC to EM units
ttm_all$xEM = -gamma*(ttm_all$x)
ttm_all$yEM = -gamma*(ttm_all$y)
ttm_all$zEM = +gamma*(ttm_all$z)

# From EM to physical units
L = 384400
ttm_all$xPH = L*ttm_all$xEM
ttm_all$yPH = L*ttm_all$yEM
ttm_all$zPH = L*ttm_all$zEM

#radius
ttm_all$rNC = sqrt(ttm_all$x^2+ttm_all$y^2)
ttm_all$rPH = sqrt(ttm_all$xPH^2+ttm_all$yPH^2)

#Select only positive s1
ttm_all = ttm_all[which(ttm_all$s1 >=0),]

#Get rid of the origin
isOrigin = ttm_all$s1 == 0 &
  ttm_all$s2 == 0 & ttm_all$s3 == 0 & ttm_all$s4 == 0
#ttm_all = ttm_all[which(!isOrigin),]


#Plot the precision as a function of (s1=s2=s3=s4)
plotdf_line(ttm_all, "s1", "log10eOm")#+scale_x_continuous(limits=c(0, 10))+scale_y_continuous(limits=c(3.5, 4.5))

plotdf_line(ttm_all, "rPH", "log10eOm")

plotdf_line(ttm_all, "rNC", "log10eOm")
plotdf_line(ttm_all, "rPH", "log10eOm")+scale_x_continuous(limits=c(200, 5000))+scale_y_continuous(limits=c(3, 4.5))

ttm_c[which(ttm_c$log10eOm > 3.9),]
