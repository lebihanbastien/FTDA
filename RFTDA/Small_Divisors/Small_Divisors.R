# R script to handle a precision map of the QBFBP around EML1,2
#---------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------
# R options
#------------------------------------------------
options(digits = 15)

#------------------------------------------------
# Load libraries
#------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)

#------------------------------------------------
# Load Source files
#------------------------------------------------
source("source/source_folder.R")
source("source/source_plot.R")
source("source/source_routines.R")
source("source/multiplot.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li    = "L2"
MODEL = "QBCP"
FWRK  = "EM"
eps = "0.5"
currentfolder = paste0(printfolder(MODEL, FWRK, Li))


#----------------------------------------------------------------------------------
# Non-Autonomous divisors
#----------------------------------------------------------------------------------
#------------------------------------------------
# Filename to check
#------------------------------------------------
fileprefix = paste0(currentfolder, "Serv/nonaut_small_divisors_",eps)
filename = paste0(fileprefix, ".txt")

#------------------------------------------------
# Load csv source
#------------------------------------------------
if (file.exists(filename))
{
  nasd  = read.csv(filename, header = T, sep = ",")
}else
{
  nasd = data.frame()
}

#------------------------------------------------
# Postprocess
#------------------------------------------------
nasd = nasd[which(nasd$Order < 16),]
nasd$Type = "Non-autonomous"

#----------------------------------------------------------------------------------
# Autonomous divisors
#----------------------------------------------------------------------------------
#------------------------------------------------
# Filename to check
#------------------------------------------------
fileprefix = paste0(currentfolder, "Serv/aut_small_divisors_",eps)
filename = paste0(fileprefix, ".txt")

#------------------------------------------------
# Load csv source
#------------------------------------------------
if (file.exists(filename))
{
  asd  = read.csv(filename, header = T, sep = ",")
}else
{
  asd = data.frame()
}

#------------------------------------------------
# Postprocess
#------------------------------------------------
asd = asd[which(asd$Order < 16),]
asd$Type = "Autonomous"

#------------------------------------------------
#rbind in sd
#------------------------------------------------
sd = rbind(asd, nasd)

#------------------------------------------------
# Histogram
#------------------------------------------------
his = ggplot()
his = his + geom_bar(data = sd, aes(x=Order, y = y, fill = factor(Type)), stat = "identity")
#Labels
his = his+labs(x = "Order of the homogeneous polynomial", y = "Number of small divisors")
#Legend
his = his+scale_fill_discrete(name="Type")
#Theme
his+custom_theme

#----------------------------------------------------------------------------------
# All cases
#----------------------------------------------------------------------------------
#------------------------------------------------
# Building the data.frame of results
#------------------------------------------------
veps = c(0.01, 0.05, 0.1, 0.5)
sdall = data.frame()
for (eps in veps)  #loop on the orders
{
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0(currentfolder, "Serv/aut_small_divisors_", toString(eps))
  filename = paste0(fileprefix, ".txt")
  # Load csv source
  #------------------------------------------------
  ac  = read.csv(filename, header = T, sep = ",")
  ac$type = " Autonomous  " 
  
  # Filename to check
  #------------------------------------------------
  fileprefix = paste0(currentfolder, "Serv/nonaut_small_divisors_", toString(eps))
  filename = paste0(fileprefix, ".txt")
  # Load csv source
  #------------------------------------------------
  nac  = read.csv(filename, header = T, sep = ",")
  nac$type = " Non-autonomous  " 
  
  #rbind in ttm_c
  ttm_c = rbind(ac, nac)
  
  #Epsilon
  ttm_c$eps = paste0("eps = ",toString(eps))
  
  #rbind in sdall
  sdall = rbind(sdall, ttm_c)
}

#------------------------------------------------
# Postprocess
#------------------------------------------------
sdall = sdall[which(sdall$Order < 16),]

#------------------------------------------------
# Histogram
#------------------------------------------------
hisall = ggplot()
#hisall = hisall + geom_bar(data = sdall, aes(x=Order, y = y, fill = factor(eps)), color = "black", stat = "identity", position="dodge")
hisall = hisall + geom_bar(data = sdall, aes(x=Order, y = y, fill = factor(type)), stat = "identity")+facet_grid(.~eps)
hisall = hisall+labs(x = "\n Order of the homogeneous polynomial", y = "Number of small divisors \n")
#Legend
hisall = hisall+scale_fill_discrete(name="Type:  ")
#Theme
hisall = hisall+custom_theme+theme(legend.position="top", legend.title.align= 0.5)
#Theme for facets
hisall = hisall+theme(strip.text.x = element_text(size=10,  face="bold"),
      strip.text.y = element_text(size=10,  face="bold"))

#------------------------------------------------
# Save in eps file
#------------------------------------------------
ggsave(hisall, file=paste0(currentfolder, "small_divisors.eps"))
