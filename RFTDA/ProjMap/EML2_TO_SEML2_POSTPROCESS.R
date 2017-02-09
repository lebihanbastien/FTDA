################################################################################
# R script used in EML2_TO_SEML.R for the postprocess 
# of a projection map (connections between EML2 and SEML1,2)
#
# WARNING: EML2_TO_SEML2.R or EML2_TO_SEML2_PROJMAP.R must be loaded once before.

# This scripts basically uses the dataframe proj_map_source created in 
# EML2_TO_SEML2_PROJMAP.R and updates the following dataframes:
#
# - proj_map: solutions of a multi projection procedure, with the projection 
#   distance under a certain value (projection_lim_max)
# - proj_map_tem: subset of solutions only with a given initial time
# - proj_map_s1: subset of solutions only with a given initial s1
# - proj_map_s3: subset of solutions only with a given initial s3
#
################################################################################

#===============================================================================
# Define limit of projection
#===============================================================================
#hard
projection_lim_max = 1e-2; #1e-1
projection_lim_mid = 5e-4;

# Projection in SI
projection_lim_max_SI = projection_lim_max*CST_DIST_PRIM_SEM
projection_lim_mid_SI = projection_lim_mid*CST_DIST_PRIM_SEM 
projection_lim_mid_SI

#===============================================================================
# Define limit of projection for colors
#===============================================================================
projection_color_lim = c(0, 1e-2);

#===============================================================================
# Indexes for cuts in data
#===============================================================================
s1_index = 50;
s3_index = 50;

#===============================================================================
#===============================================================================
# Postprocess on proj_map_source (new variables, etc)
#===============================================================================
#===============================================================================

#===============================================================================
#Errors and norms
#===============================================================================
#Get the error in position in km
proj_map_source$pmin_dist_SI = proj_map_source$pmin_dist_SEM*CST_DIST_PRIM_SEM
#Norm of yv
proj_map_source$nf_CM_NCSEM = sqrt((-1/CST_GAMMA_LIB_SEM*proj_map_source$xf_CM_SEM + CST_C1_LIB_SEM)^2 + (1/CST_GAMMA_LIB_SEM*proj_map_source$yf_CM_SEM)^2)
#Get the distance wrt to the center of s1_CM_SEM/s3_CM_SEM
proj_map_source$ns_CM_SEM = sqrt(proj_map_source$s1_CM_SEM^2 + proj_map_source$s3_CM_SEM^2)
#Get the distance wrt to the center of s1EM/s3EM
proj_map_source$ns_CMU_EM = abs(proj_map_source$s1_CMU_EM)
#Position gap on axes
proj_map_source$dx_CM_SEM  = proj_map_source$xf_CM_SEM  - proj_map_source$xp_CM_SEM
proj_map_source$dy_CM_SEM  = proj_map_source$yf_CM_SEM  - proj_map_source$yp_CM_SEM
proj_map_source$dpx_CM_SEM = proj_map_source$pxf_CM_SEM - proj_map_source$pxp_CM_SEM
proj_map_source$dpy_CM_SEM = proj_map_source$pyf_CM_SEM - proj_map_source$pyp_CM_SEM

#===============================================================================
#DV
#===============================================================================
#Get the dv_at_projection_SEM in position in km/s
proj_map_source$dv_at_projection_SI = proj_map_source$dv_at_projection_SEM*CST_DIST_PRIM_SEM*2*pi/Tcrtbp("SEM")

#===============================================================================
#Minimum in position
#===============================================================================
proj_map_source_min = proj_map_source[which(proj_map_source$pmin_dist_SEM == min(proj_map_source$pmin_dist_SEM)),] 
#Associated position gap in km
proj_map_source_min$pmin_dist_SI
#Associated dv_at_projection_SEM in m/s
1e3*proj_map_source_min$dv_at_projection_SI

#===============================================================================
#Radius
#===============================================================================
proj_map_source$rf_CM_SEM = sqrt(proj_map_source$xf_CM_SEM^2 + proj_map_source$yf_CM_SEM^2 + proj_map_source$zf_CM_SEM^2)
rf_CM_SEM_mid = mean(proj_map_source$rf_CM_SEM)

proj_map_source$sf_CM_SEM = sqrt(proj_map_source$s1_CM_SEM^2 + proj_map_source$s3_CM_SEM^2)
sf_CM_SEM_mid = mean(proj_map_source$sf_CM_SEM)

#===============================================================================
#Time
#===============================================================================
#Get the time as a percentage of T
proj_map_source$t0_CMU_EMT = proj_map_source$t0_CMU_EM/CST_SEM_PERIOD_EM
#Get the time as a percentage of T
proj_map_source$tf_man_SEMT = proj_map_source$tf_man_SEM/CST_SEM_PERIOD_SEM
#tof_EM (in T)
proj_map_source$tof_EM = (proj_map_source$tf_man_SEMT-proj_map_source$t0_CMU_EMT)
#tof_EM (days)
proj_map_source$tof_SI = (proj_map_source$tf_man_SEMT-proj_map_source$t0_CMU_EMT)*SEMperiod("EM")*Tcrtbp("EM")/(2*pi*3600*24)

#===============================================================================
#NCEM to EM
#===============================================================================
proj_map_source$x0_CMU_EM <- -CST_GAMMA_LIB_EM * (proj_map_source$x0_CMU_NCEM-CST_C1_LIB_EM)
proj_map_source$y0_CMU_EM <- -CST_GAMMA_LIB_EM * (proj_map_source$y0_CMU_NCEM)
proj_map_source$z0_CMU_EM <- +CST_GAMMA_LIB_EM * (proj_map_source$z0_CMU_NCEM)


#===============================================================================
#===============================================================================
# Select a subgroup (only certain distance of projection, time, etc)
#===============================================================================
#===============================================================================

#===============================================================================
#Keep only values with a projection distance below projection_lim_max.
#===============================================================================
proj_map = proj_map_source[which(proj_map_source$pmin_dist_SEM <= projection_lim_max),] 

#===============================================================================
# Select a given value of time
#===============================================================================
t0_CMU_EM_vec = sort(unique(proj_map$t0_CMU_EM))

# If time_desired > 0, we use it to define time_index
if(time_desired >= 0)
{
  time_index =  which.min(abs(t0_CMU_EM_vec-time_desired))
}

if(time_index > 0)
{
  proj_map_tem = proj_map[which(proj_map$t0_CMU_EM == t0_CMU_EM_vec[time_index]),] 
}else{
  proj_map_tem = proj_map
}

#===============================================================================
# Just the collisions within tem
#===============================================================================
proj_map_tem_coll = proj_map_tem[which(proj_map_tem$dv_at_projection_SEM %in% c(301, 399)),]


#===============================================================================
# Select a given value of s1
#===============================================================================
s1_CMU_EM_vec = sort(unique(proj_map$s1_CMU_EM))

if(s1_index > 0)
{
  #proj_map_s1 = proj_map[which(proj_map$s1_CMU_EM == s1_CMU_EM_vec[s1_index]),] 
  proj_map_s1 = proj_map[which(proj_map$s1_CMU_EM == 0.0),] 
}else{
  proj_map_s1 = proj_map
}

#===============================================================================
# Select a given value of s3
#===============================================================================
s3_CMU_EM_vec = sort(unique(proj_map$s3_CMU_EM))

if(s3_index > 0)
{
  #proj_map_s3 = proj_map[which(proj_map$s3_CMU_EM == s3_CMU_EM_vec[s3_index]),] 
  proj_map_s3 = proj_map[which(proj_map$s3_CMU_EM == 0.0),] 
}else{
  proj_map_s3 = proj_map
}



