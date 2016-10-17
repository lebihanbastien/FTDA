######################################################################################
# R script used in EML2_TO_SEML.R for the video
# of a projection map (connections between EML2 and SEML1,2)
######################################################################################
#===============================================================================
# Define limit of projection
#===============================================================================
#hard
projection_lim_max = 1e0;
projection_lim_mid = 0.0;

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
index = 0;
for(t0 in t0_CMU_EM_vec)
{
  #=================================
  # Select the data
  #=================================
  proj_map_tem = proj_map[which(proj_map$t0_CMU_EM == t0),] 

  #=================================
  # Plot
  #=================================
  pp_tiles_x0EM_y0EM_eP = plotdf_point(proj_map_tem, "x0_CMU_EM", "y0_CMU_EM", expression("x"[0]^"em"), expression("y"[0]^"em"), "pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 3)
  pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + scale_colour_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)
  pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + scale_x_continuous(limits = c(-0.3, 0.3)) 
  pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + scale_y_continuous(limits = c(-0.5, 0.5)) 
  
  #=================================
  # Title
  #=================================
  pp_tiles_x0EM_y0EM_eP   = pp_tiles_x0EM_y0EM_eP + ggtitle(paste0("t = ", toString(t0/CST_SEM_PERIOD_EM), "T"))
  
  #=================================
  # Save
  #=================================
  #filename = paste0(getwd(), "/ProjMap/VIDEO/", "pp_tiles_x0EM_y0EM_eP_t0_", sprintf("%2.4f", t0), '.pdf')
  filename = paste0(getwd(), "/ProjMap/VIDEO/", "pp_tiles_x0EM_y0EM_eP_t0_", sprintf("%i", index), '.png')
  ggsave(pp_tiles_x0EM_y0EM_eP, width = xSize, height = ySize, file = filename)
  
  index = index +1;
}