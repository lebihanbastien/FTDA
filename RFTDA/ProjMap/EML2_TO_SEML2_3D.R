######################################################################################
# R script to handle a projection map (connections between EML2 and SEML1,2)
######################################################################################

#===============================================================================
# Initialization
#===============================================================================
source("source/init.R")

#------------------------------------------------
# Select paramaters
#------------------------------------------------
LIB_POINT_EM = "L2"
DYN_MODEL    = "QBCP"
FWRK         = "EM"
DATA_SOURCE  = "LOCAL" #FROM_SERVER or LOCAL
ORDER        = switch(DATA_SOURCE, "LOCAL" = "16", "FROM_SERVER" = "20");


#------------------------------------------------
# Constraints on the center/center-unstable manifolds
#------------------------------------------------
SNMAX    = 0.6;
YNMAX    = 0.6;
GSIZEMAX = 35;

#------------------------------------------------
# Which type of data?
#------------------------------------------------
# Subfolder
#-------------------------
FILE_SUBFOLDER = switch(DATA_SOURCE, "LOCAL" = "plot/QBCP/EM/L2/", "FROM_SERVER" = "plot/QBCP/EM/L2/Serv/")

# From server or not?
#-------------------------
FILE_PREFIX     = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "projcu_3d_order_", ORDER),
                     "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "projcu_3d_order_", ORDER) )

# Which data file?
#-------------------------
###### NEW ###### 
FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = ".bin", "FROM_SERVER" = ".bin" )

#------------------------------------------------
# Constants
#------------------------------------------------
# EM framework
CST_MASS_RATIO_EM  = muR("EM");
CST_GAMMA_LIB_EM   = gamma(LIB_POINT_EM, "EM");
CST_C1_LIB_EM      = c1(LIB_POINT_EM, "EM");
CST_DIST_PRIM_EM   = Ldist("EM");
CST_SEM_PERIOD_EM  = SEMperiod("EM");

# SEM framework
CST_DIST_PRIM_SEM  = Ldist("SEM");
CST_GAMMA_LIB_SEM  = gamma("L2", "SEM");
CST_C1_LIB_SEM     = c1("L2", "SEM");
CST_SEM_PERIOD_SEM = SEMperiod("SEM");
CST_NC_TO_SI_FACTOR_SEM = CST_DIST_PRIM_SEM*CST_GAMMA_LIB_SEM;

#===============================================================================
# Get data from file
#===============================================================================
#------------------------------------------------
# Projection map
#------------------------------------------------
filename = paste0(FILE_PREFIX, FILE_SUFFIX);
if (file.exists(filename))
{
  names = c("t0_CMU_EM",  "x0_CMU_EM", "y0_CMU_EM", "z0_CMU_EM", "px0_CMU_EM", "py0_CMU_EM", "pz0_CMU_EM",
            "x0_CM_SEM", "y0_CM_SEM", "z0_CM_SEM", "px0_CM_SEM", "py0_CM_SEM", "pz0_CM_SEM",
            "s1_CMU_EM", "s2_CMU_EM", "s3_CMU_EM", "s4_CMU_EM",  "s5_CMU_EM",  "pmin_dist_SEM", "dv_at_projection_SEM", "tf_man_SEM",
            "xf_CM_SEM", "yf_CM_SEM", "zf_CM_SEM", "pxf_CM_SEM", "pyf_CM_SEM", "pzf_CM_SEM",
            "xp_CM_SEM", "yp_CM_SEM", "zp_CM_SEM", "pxp_CM_SEM", "pyp_CM_SEM", "pzp_CM_SEM",
            "s1_CM_SEM", "s2_CM_SEM", "s3_CM_SEM", "s4_CM_SEM");
  proj_map_source = dffbinary(filename, 37, names)
}else
{
  proj_map_source = data.frame()
}


#------------------------------------------------
# Best solution
#------------------------------------------------
proj_map_sol = data.frame()


#------------------------------------------------
# Continuation
#------------------------------------------------
proj_map_cont = data.frame()

#------------------------------------------------
# Continuation trajectories
#------------------------------------------------
proj_map_cont_traj = data.frame()


#===============================================================================
# Postprocess
#===============================================================================
source("ProjMap/EML2_TO_SEML2_POSTPROCESS_3D.R")

#===============================================================================
# PLOTS
#===============================================================================
source("ProjMap/EML2_TO_SEML2_PLOTS_3D.R")

#===============================================================================
# MULTIPLOTS
#===============================================================================
#source("ProjMap/EML2_TO_SEML2_MULTIPLOTS.R")
