######################################################################################
# R script used in EML2_TO_SEML.R for the postprocess 
# of a projection map (connections between EML2 and SEML1,2)
######################################################################################
#------------------------------------------------
# Define limit of projection
#------------------------------------------------
#hard
projection_lim_max = Inf;
projection_lim_mid = 5e-4;

#===============================================================================
# Define limit of projection for colors
#===============================================================================
projection_color_lim = c(0, 5e-3);

#easy
# projection_lim_max = 5e-3;
# projection_lim_mid = 1e-3;

# Projection in SI
projection_lim_max_SI = projection_lim_max*CST_DIST_PRIM_SEM
projection_lim_mid_SI = projection_lim_mid*CST_DIST_PRIM_SEM 
projection_lim_mid_SI


#------------------------------------------------
# Time index
#------------------------------------------------
time_index = -1;


#===============================================================================
#Keep only values with a projection distance below projection_lim_max.
#===============================================================================
proj_map = proj_map_source[which(proj_map_source$pmin_dist_SEM <= projection_lim_max),] 

#===============================================================================
#Errors and norms
#===============================================================================
#Get the error in position in km
proj_map$pmin_dist_SI = proj_map$pmin_dist_SEM*CST_DIST_PRIM_SEM
#Norm of yv
proj_map$nf_CM_NCSEM = sqrt((-1/CST_GAMMA_LIB_SEM*proj_map$xf_CM_SEM + CST_C1_LIB_SEM)^2 + (1/CST_GAMMA_LIB_SEM*proj_map$yf_CM_SEM)^2)
#Get the distance wrt to the center of s1_CM_SEM/s3_CM_SEM
proj_map$ns_CM_SEM = sqrt(proj_map$s1_CM_SEM^2 + proj_map$s3_CM_SEM^2)
#Get the distance wrt to the center of s1EM/s3EM
proj_map$ns_CMU_EM = abs(proj_map$s1_CMU_EM)
#Position gap on axes
proj_map$dx_CM_SEM  = proj_map$xf_CM_SEM  - proj_map$xp_CM_SEM
proj_map$dy_CM_SEM  = proj_map$yf_CM_SEM  - proj_map$yp_CM_SEM
proj_map$dpx_CM_SEM = proj_map$pxf_CM_SEM - proj_map$pxp_CM_SEM
proj_map$dpy_CM_SEM = proj_map$pyf_CM_SEM - proj_map$pyp_CM_SEM

#===============================================================================
#DV
#===============================================================================
#Get the dv_at_projection_SEM in position in km/s
proj_map$dv_at_projection_SI = proj_map$dv_at_projection_SEM*CST_DIST_PRIM_SEM*2*pi/Tcrtbp("SEM")

#===============================================================================
#Minimum in position
#===============================================================================
proj_map_min = proj_map[which(proj_map$pmin_dist_SEM == min(proj_map$pmin_dist_SEM)),] 
#Associated position gap in km
proj_map_min$pmin_dist_SI
#Associated dv_at_projection_SEM in m/s
1e3*proj_map_min$dv_at_projection_SI

#===============================================================================
#Time
#===============================================================================
#Get the time as a percentage of T
proj_map$t0_CMU_EMT = proj_map$t0_CMU_EM/CST_SEM_PERIOD_EM
#Get the time as a percentage of T
proj_map$tf_man_SEMT = proj_map$tf_man_SEM/CST_SEM_PERIOD_SEM
#tof_EM (in T)
proj_map$tof_EM = (proj_map$tf_man_SEMT-proj_map$t0_CMU_EMT)
#tof_EM (days)
proj_map$tof_SI = (proj_map$tf_man_SEMT-proj_map$t0_CMU_EMT)*SEMperiod("EM")*Tcrtbp("EM")/(2*pi*3600*24)
#Keep only a certain time
#proj_map = proj_map[which(proj_map$t0_CMU_EMT == proj_map$t0_CMU_EMT[8]),]

#===============================================================================
# Unique
#===============================================================================
s1_CMU_EM_vec = unique(proj_map$s1_CMU_EM)
s2_CMU_EM_vec = unique(proj_map$s2_CMU_EM)
s3_CMU_EM_vec = unique(proj_map$s3_CMU_EM)
s4_CMU_EM_vec = unique(proj_map$s4_CMU_EM)

t0_CMU_EM_vec = unique(proj_map$t0_CMU_EM)

#===============================================================================
# Select a given value of time
#===============================================================================
if(time_index > 0)
{
  proj_map_tem = proj_map[which(proj_map$t0_CMU_EM == t0_CMU_EM_vec[time_index]),] 
}else{
  proj_map_tem = proj_map
}

