#--------------------------------------------------------------------------------------
# Potential of the CR3BP
#--------------------------------------------------------------------------------------
# Init
#---------------------
source("source/init.R")
FWRK  = "EM"
muR = muR(FWRK);

#---------------------
#Build vectors
#---------------------
x0   = seq(-2,-1,0.001)
y0   = seq(-0.15,0.15,0.001)

z0   = expand.grid(x0, y0)
potn = apply(z0, 1, cr3bpot_EM)


#Build data frame
df = data.frame(x1 = z0[1], y1 = z0[2], pot = potn)
names(df) = c('x', 'y', 'pot')

#---------------------
#Energy associated with the Libration points
#---------------------
l1 = c(-0.836915145386502, 0, 0)
l2 = c(-1.155682150113637, 0, 0)
l1pot = cr3bpot(l1[1], l1[2], muR)
l2pot = cr3bpot(l2[1], l2[2], muR)

#---------------------
#Select a given energy
#---------------------
eps = l2pot+0.01;
isCloseToEPS = abs(df$pot - eps) < 5e-4;
dfr = df[which(isCloseToEPS),]

#---------------------
# Plot
#---------------------
ppl = plotdf_point(dfr, "x", "y", "X", "Y", pointSize = 3)
ppl
