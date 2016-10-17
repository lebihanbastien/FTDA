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
FILE_PREFIX     = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "projcu_order_", ORDER),
                     "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "projcu_order_", ORDER) )
FILE_PREFIX_SOL = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "sortprojintcu_order_", ORDER),
                     "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "sortprojintcu_order_", ORDER) )

FILE_PREFIX_CONT      = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "cont_atf_order_", ORDER),
                         "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "cont_atf_order_", ORDER) )
FILE_PREFIX_CONT_TRAJ = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "cont_atf_traj_order_", ORDER),
                          "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "cont_atf_traj_order_", ORDER) )

# Which data file?
#-------------------------
###### NEW ###### 
#FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "_tspan_0T_025T_FINAL" )
#FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "_tspan_075T_T_FINAL" )
#FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "_tspan_05T_075T_FINAL" )
#FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "_tspan_025T_05T_FINAL" )
FILE_SUFFIX = switch(DATA_SOURCE, "LOCAL" = "", "FROM_SERVER" = "" )

###### OLD ###### 
#FILE_SUFFIX = "_tspan_0_T"
#FILE_SUFFIX = "_tspan_083T_097T"
#FILE_SUFFIX = "_tspan_053T_061T"
#FILE_SUFFIX = "_tspan_0942T"

#### LOCAL ####
#FILE_SUFFIX = "_tspan_0T"

# Suffix for continuation file
#-------------------------
FILE_SUFFIX_CONT      = "_t0_0.99"
FILE_SUFFIX_CONT_TRAJ = "_t0_0.99"

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
# Best solution
#------------------------------------------------
filepre  = paste0(FILE_PREFIX_SOL, FILE_SUFFIX);
filename = paste0(filepre, ".bin");
if (file.exists(filename))
{
  names = c("label", "s1_CMU_EM", "s3_CMU_EM", "s1_CM_SEM", "s3_CM_SEM",  "pmin_dist_SEM",
            "t0_CMU_EM", "x_man_SEM", "y_man_SEM", "z_man_SEM", 
            "px_man_SEM", "py_man_SEM", "pz_man_SEM", "h_man_SEM", 
            "t_orb_eml_SEM", "x_orb_eml_SEM", "y_orb_eml_SEM", "z_orb_eml_SEM", 
            "px_orb_eml_SEM", "py_orb_eml_SEM", "pz_orb_eml_SEM", "h_orb_eml_SEM",
            "t_orb_seml_SEM", "x_orb_seml_SEM", "y_orb_seml_SEM", "z_orb_seml_SEM", 
            "px_orb_seml_SEM", "py_orb_seml_SEM", "pz_orb_seml_SEM", "h_orb_seml_SEM");
  
  # OLD VERSION (SERVER)
  # names = c("label", "t0_CMU_EM", "x_man_SEM", "y_man_SEM", "z_man_SEM", "px_man_SEM", "py_man_SEM", "pz_man_SEM",
  #           "s1_CMU_EM", "s3_CMU_EM", "s1_CM_SEM", "s3_CM_SEM", "pmin_dist_SEM", "h_man_SEM", 
  #           "t_orb_eml_SEM", "x_orb_eml_SEM", "y_orb_eml_SEM", "z_orb_eml_SEM", 
  #           "px_orb_eml_SEM", "py_orb_eml_SEM", "pz_orb_eml_SEM", "h_orb_eml_SEM",
  #           "t_orb_seml_SEM", "x_orb_seml_SEM", "y_orb_seml_SEM", "z_orb_seml_SEM", "px_orb_seml_SEM", "py_orb_seml_SEM", "pzO2", "h_orb_seml_SEM");
  
  proj_map_sol = dffbinary(filename, 30, names)
}else
{
  proj_map_sol = data.frame()
}

#------------------------------------------------
# Projection map
#------------------------------------------------
filename = paste0(FILE_PREFIX, FILE_SUFFIX, ".bin");
if (file.exists(filename))
{
  names = c("t0_CMU_EM", "x0_CMU_NCEM", "y0_CMU_NCEM", "z0_CMU_NCEM", "px0_CMU_NCEM", "py0_CMU_NCEM", "pz0_CMU_NCEM",
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
# Continuation
#------------------------------------------------
filename_cont = paste0(FILE_PREFIX_CONT, FILE_SUFFIX_CONT, ".txt");
if (file.exists(filename_cont))
{
  proj_map_cont = read.table(file = filename_cont, header = TRUE)
  
  #Small postprocess
  proj_map_cont$x0_CMS_SEM = -CST_GAMMA_LIB_SEM*(proj_map_cont$x0_CMS_NCSEM - CST_C1_LIB_SEM)
  proj_map_cont$y0_CMS_SEM = -CST_GAMMA_LIB_SEM*(proj_map_cont$y0_CMS_NCSEM - 0)
  proj_map_cont$z0_CMS_SEM = +CST_GAMMA_LIB_SEM*(proj_map_cont$z0_CMS_NCSEM - 0)
}else
{
  proj_map_cont = data.frame()
}

#------------------------------------------------
# Continuation trajectories
#------------------------------------------------
filepre_cont_tra = paste0(FILE_PREFIX_CONT_TRAJ, FILE_SUFFIX_CONT_TRAJ);
filename_cont_traj = paste0(filepre_cont_tra, ".bin");
if (file.exists(filename_cont_traj))
{
  #proj_map_cont_traj = read.table(file = filename_cont_traj, header = TRUE)
  
  #label t_CMU_SEM x_CMS_NCSEM y_CMS_NCSEM z_CMS_NCSEM px_CMS_NCSEM py_CMS_NCSEM pz_CMS_NCSEM 
  names = c("label",  "t_CMU_SEM", "x_CMS_NCSEM", "y_CMS_NCSEM", "z_CMS_NCSEM", "px_CMS_NCSEM", "py_CMS_NCSEM", "pz_CMS_NCSEM");
  proj_map_cont_traj = dffbinary(filename_cont_traj, 8, names)
  
  #Small postprocess
  proj_map_cont_traj$x_CMS_SEM = -CST_GAMMA_LIB_SEM*(proj_map_cont_traj$x_CMS_NCSEM - CST_C1_LIB_SEM)
  proj_map_cont_traj$y_CMS_SEM = -CST_GAMMA_LIB_SEM*(proj_map_cont_traj$y_CMS_NCSEM - 0)
  proj_map_cont_traj$z_CMS_SEM = +CST_GAMMA_LIB_SEM*(proj_map_cont_traj$z_CMS_NCSEM - 0)
  
  #Label only
  proj_map_cont_label = proj_map_cont_traj[which(proj_map_cont_traj$t_CMU_SEM == min(proj_map_cont_traj$t_CMU_SEM)),]
  
  #Copy label into 
  proj_map_cont$label = proj_map_cont_label$label
  
}else
{
  proj_map_cont_traj = data.frame()
}

#===============================================================================
# Postprocess
#===============================================================================
source("ProjMap/EML2_TO_SEML2_POSTPROCESS.R")

#===============================================================================
# PLOTS
#===============================================================================
source("ProjMap/EML2_TO_SEML2_PLOTS_LATEX.R")
#source("ProjMap/EML2_TO_SEML2_PLOTS.R")

#===============================================================================
# MULTIPLOTS
#===============================================================================
#source("ProjMap/EML2_TO_SEML2_MULTIPLOTS.R")
