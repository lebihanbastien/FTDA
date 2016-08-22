# General test on the QBTBP vector field
#---------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------
# Init
#------------------------------------------------
source("source/init.R")

#------------------------------------------------
# Select Models & libration point
#------------------------------------------------
Li    = "L2"
MODEL = "QBCP"
FWRK  = "EM"

#------------------------------------------------
#Normalized units (gamma_li, c1_li)
#------------------------------------------------
muR      = muR(FWRK);
gamma_li = gamma(Li, FWRK);
c1_li    =  c1(Li, FWRK);
L        = Ldist(FWRK);
Tsem     = SEMperiod(FWRK);


#------------------------------------------------
# Get data from file
#------------------------------------------------
prefix = "QBTBP_EM_vs_IN"

dfe = data.frame();

#Earth
suffix = "E.txt"
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", prefix, "_", suffix);
res = read.table(filename)
colnames(res) = c("t", "Error")
res$Primary = "Earth"
dfe = rbind(dfe, res)

#Moon
suffix = "M.txt"
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", prefix, "_", suffix);
res = read.table(filename)
colnames(res) = c("t", "Error")
res$Primary = "Moon"
dfe = rbind(dfe, res)


#Sun
suffix = "S.txt"
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", prefix, "_", suffix);
res = read.table(filename)
colnames(res) = c("t", "Error")
res$Primary = "Sun"
dfe = rbind(dfe, res)

#------------------------------------------------
# Plots
#------------------------------------------------
#Initial EML2 orbit
ppe = plotdf_line(dfe, "t", "Error", "time [$\\times T$]", "State vector error [EM units]", "Primary", "Primary", isColorFac = TRUE)
ppe = ppe + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))
ppe = ppe + legend_inside_theme
ppe

#------------------------------------------------
# Save
#------------------------------------------------
ggplot2tikz(ppe, xSize, ySize, file = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", prefix, ".tex"))