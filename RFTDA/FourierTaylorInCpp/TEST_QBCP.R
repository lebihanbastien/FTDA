# General test on the QBCP (four-body) vector field
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
# Get data from file (EM)
#------------------------------------------------
dfem = data.frame();

#NCEM
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", "QBCP_NCEM_vs_IN.txt");
res = read.table(filename)
colnames(res) = c("t", "Error")
res$type = "EMNC"
res$fwrk = "EM"
dfem = rbind(dfem, res)

#VNCEM
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", "QBCP_VNCEM_vs_IN.txt");
res = read.table(filename)
colnames(res) = c("t", "Error")
res$type = "EMNCV"
res$fwrk = "EM"
dfem = rbind(dfem, res)

#EM
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", "QBCP_EM_vs_IN.txt");
res = read.table(filename)
colnames(res) = c("t", "Error")
res$type = "EM"
res$fwrk = "EM"
dfem = rbind(dfem, res)


#------------------------------------------------
# Get data from file (SEM)
#------------------------------------------------
#NCSEM
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", "QBCP_NCSEM_vs_IN.txt");
res = read.table(filename)
colnames(res) = c("t", "Error")
res$type = "SEMNC"
res$fwrk = "SEM"
dfem = rbind(dfem, res)

#VNCSEM
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", "QBCP_VNCSEM_vs_IN.txt");
res = read.table(filename)
colnames(res) = c("t", "Error")
res$type = "SEMNCV"
res$fwrk = "SEM"
dfem = rbind(dfem, res)

#SEM
filename = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", "QBCP_SEM_vs_IN.txt");
res = read.table(filename)
colnames(res) = c("t", "Error")
res$type = "SEM"
res$fwrk = "SEM"
dfem = rbind(dfem, res)




#------------------------------------------------
# Plots
#------------------------------------------------
ppe = ggplot(dfem, aes(x=t, y=Error, group = type))
ppe = ppe + geom_line(aes(linetype=type, color = type), size = 1.5)
ppe = ppe + scale_linetype_manual(values=c("solid", "solid", "solid",  "dashed",  "dashed",  "dashed")) # Change linetypes
ppe = ppe + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))
ppe = ppe + legend_inside_theme
ppe

#------------------------------------------------
# Save
#------------------------------------------------
ggplot2tikz(ppe, xSize, ySize, file = paste0("/home/b.le-bihan/BackUpBox/PhD/FourierTaylorInCpp/plot/QBCP/EM/L2/", "QBCP_ALL_vs_IN.tex"))