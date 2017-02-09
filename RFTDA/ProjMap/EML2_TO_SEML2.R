################################################################################
# R script to handle a projection map (connections between EML2 and SEML1,2)
################################################################################

#===============================================================================
# Initialization
#===============================================================================
source("source/init.R")

#-------------------------------------------------------------------------------
# Select parameters
#-------------------------------------------------------------------------------
LIB_POINT_EM  = "L2"
LIB_POINT_SEM = "L2"
DYN_MODEL     = "QBCP"
FWRK          = "EM"
DATA_SOURCE   = "FROM_SERVER" #FROM_SERVER or LOCAL
ORDER         =  switch(DATA_SOURCE, "LOCAL" = "16", "FROM_SERVER" = "20");
LOCAL_ORDER   = "20";

#-------------------------------------------------------------------------------
# Constraints on the center/center-unstable manifolds
#-------------------------------------------------------------------------------
SNMAX    = 0.6;
YNMAX    = 0.6;
GSIZEMAX = 35;

#-------------------------------------------------------------------------------
# Which type of data?
#-------------------------------------------------------------------------------
# Subfolder
#-------------------------
FILE_SUBFOLDER = switch(DATA_SOURCE, "LOCAL" = "plot/QBCP/EM/L2/", "FROM_SERVER" = "plot/QBCP/EM/L2/Serv/")

# Prefix
#-------------------------
FILE_PREFIX     = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "projcu_order_", ORDER, "_dest_", LIB_POINT_SEM),
                     "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "projcu_order_", ORDER) )
FILE_PREFIX_SOL = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "sortprojintcu_order_", ORDER, "_dest_", LIB_POINT_SEM),
                     "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "sortprojintcu_order_", ORDER) )


# Which data file?
#-------------------------
###### NEW ###### 
#FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "_tspan_075T_T_FINAL" )
FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "_tspan_05T_075T_FINAL" )
#FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "_tspan_025T_05T_FINAL" )
#FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "_tspan_0T_025T_FINAL" )
#FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "" )

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
# EM framework
CST_MASS_RATIO_EM  = muR("EM");
CST_GAMMA_LIB_EM   = gamma(LIB_POINT_EM, "EM");
CST_C1_LIB_EM      = c1(LIB_POINT_EM, "EM");
CST_DIST_PRIM_EM   = Ldist("EM");
CST_SEM_PERIOD_EM  = SEMperiod("EM");

# SEM framework
CST_DIST_PRIM_SEM  = Ldist("SEM");
CST_GAMMA_LIB_SEM  = gamma(LIB_POINT_SEM, "SEM");
CST_C1_LIB_SEM     = c1(LIB_POINT_SEM, "SEM");
CST_SEM_PERIOD_SEM = SEMperiod("SEM");
CST_NC_TO_SI_FACTOR_SEM = CST_DIST_PRIM_SEM*CST_GAMMA_LIB_SEM;

#===============================================================================
# Get data from file
#===============================================================================
source("ProjMap/EML2_TO_SEML2_PROJMAP.R")

#===============================================================================
# Continuation
#===============================================================================
source("ProjMap/EML2_TO_SEML2_CONTINUATION.R")
#source("ProjMap/EML2_TO_SEML2_CONTINUATION_MULTI.R")

#===============================================================================
# Postprocess
#===============================================================================
source("ProjMap/EML2_TO_SEML2_POSTPROCESS.R")

#===============================================================================
# PLOTS
#===============================================================================
source("ProjMap/EML2_TO_SEML2_PLOTS.R")
#source("ProjMap/EML2_TO_SEML2_PLOTS_LATEX.R")

#===============================================================================
# MULTIPLOTS
#===============================================================================
#source("ProjMap/EML2_TO_SEML2_MULTIPLOTS.R")
