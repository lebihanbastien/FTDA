######################################################################################
# R script used in EML2_TO_SEML.R for the video
# of a projection map (connections between EML2 and SEML1,2)
######################################################################################
#===============================================================================
# Define limit of projection
#===============================================================================
#hard
projection_lim_max = 1e0;
projection_lim_mid = 8e-4;

#===============================================================================
# Define limit of projection for colors
#===============================================================================
projection_color_lim = c(0, 5e-3);

#===============================================================================
#Keep only values with a projection distance below projection_lim_max.
#===============================================================================
proj_map = proj_map_source[which(proj_map_source$pmin_dist_SEM <= projection_lim_max),] 

#===============================================================================
# Select all unique values of time
#===============================================================================
t0_CMU_EM_vec = unique(proj_map$t0_CMU_EM)

#===============================================================================
# Loop
#===============================================================================
for(t0 in t0_CMU_EM_vec)
{
  #=================================
  # Select the data
  #=================================
  proj_map_tem = proj_map[which(proj_map$t0_CMU_EM == t0),] 

  #=================================
  # Plot
  #=================================
  pp_tiles_s1EM_s3EM_eP = plotdf_tile_n(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", expression("s"[1]), expression("s"[3]), "pmin_dist_SEM", "Projection \ndistance", projection_lim_mid, FALSE, colorLimits = projection_color_lim)
  pp_tiles_s1EM_s3EM_eP = pp_tiles_s1EM_s3EM_eP + scale_x_continuous(limits = c(-35, 35)) 
  pp_tiles_s1EM_s3EM_eP = pp_tiles_s1EM_s3EM_eP + scale_y_continuous(limits = c(-35, 35)) 
  
  #=================================
  # Title
  #=================================
  pp_tiles_s1EM_s3EM_eP   = pp_tiles_s1EM_s3EM_eP + ggtitle(paste0("t = ", toString(t0/CST_SEM_PERIOD_EM), "T"))
  
  #=================================
  # Save
  #=================================
  filename = paste0(getwd(), "/ProjMap/VIDEO/", "pp_tiles_s1EM_s3EM_eP_t0_", sprintf("%2.4f", t0), '.pdf')
  ggsave(pp_tiles_s1EM_s3EM_eP, width = xSize, height = ySize, file = filename)
}