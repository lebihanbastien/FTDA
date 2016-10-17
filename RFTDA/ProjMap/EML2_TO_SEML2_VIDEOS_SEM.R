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
for(t0 in t0_CMU_EM_vec)
{
  #=================================
  # Select the data
  #=================================
  proj_map_tem = proj_map[which(proj_map$t0_CMU_EM == t0),] 

  #=================================
  # Plot
  #=================================
  pp_pts_x0SEM_y0SEM_eP = plotdf_point(proj_map_tem, "xf_CM_SEM", "yf_CM_SEM", expression("X"[f]^"sem"), expression("Y"[f]^"sem"), "pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 3)
  pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + scale_colour_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)
  pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + scale_x_continuous(limits = c(-1.018, -1.004)) 
  pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + scale_y_continuous(limits = c(-0.01, 0.01)) 
  
  #=================================
  #Add SEMLi
  #=================================
  seml2 = data.frame(x = 0.0, y = 0.0, z = 0.0);
  seml2 = NCtoSYS(seml2, SEML2$GAMMA, SEML2$C1);
  pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_point(data = seml2, aes(x= xEM, y = yEM), size = 4) 
  
  #=================================
  # Title
  #=================================
  pp_pts_x0SEM_y0SEM_eP   = pp_pts_x0SEM_y0SEM_eP + ggtitle(paste0("t = ", toString(t0/CST_SEM_PERIOD_EM), "T"))
  
  #=================================
  # Save
  #=================================
  filename = paste0(getwd(), "/ProjMap/VIDEO/", "pp_pts_x0SEM_y0SEM_eP_t0_", sprintf("%2.4f", t0), '.pdf')
  ggsave(pp_pts_x0SEM_y0SEM_eP, width = xSize, height = ySize, file = filename)
}