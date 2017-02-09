######################################################################################
# R script used in EML2_TO_SEML.R for the video
# of a projection map (connections between EML2 and SEML1,2)
######################################################################################

#===============================================================================
# Define prefix for continuation datafiles
#===============================================================================
FILE_PREFIX_CONT = switch(DATA_SOURCE, 
                          "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, 
                                           "cont_atf_order_", ORDER, "_dest_", 
                                           LIB_POINT_SEM),
                          "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, 
                                                 "cont_atf_order_", ORDER, 
                                                 "_dest_", LIB_POINT_SEM) 
                          )
#===============================================================================
# Define limit of projection
#===============================================================================
#hard
projection_lim_max = 1e-2;
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
#NCEM to EM
#===============================================================================
proj_map$x0_CMU_EM <- -CST_GAMMA_LIB_EM * (proj_map$x0_CMU_NCEM-CST_C1_LIB_EM)
proj_map$y0_CMU_EM <- -CST_GAMMA_LIB_EM * (proj_map$y0_CMU_NCEM)
proj_map$z0_CMU_EM <- +CST_GAMMA_LIB_EM * (proj_map$z0_CMU_NCEM)

#===============================================================================
#L2
#===============================================================================
eml2 = data.frame(x0_CMU_NCEM = 0.0, y0_CMU_NCEM = 0.0, z0_CMU_NCEM = 0.0);
eml2 = NCtoSYS(seml2, CST_GAMMA_LIB_EM, CST_C1_LIB_EM);



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
  #=============================================================================
  # Select the data
  #=============================================================================
  proj_map_tem = proj_map[which(proj_map$t0_CMU_EM == t0),] 
  
  #=============================================================================
  # Select the continuation data if they exist
  #=============================================================================
  ratio_desired = t0/CST_SEM_PERIOD_EM;
  FILE_SUFFIX_CONT      = paste0("_t0_", ratio_desired);
  proj_map_cont = get_proj_map_cont(FILE_PREFIX_CONT, FILE_SUFFIX_CONT)
  
  #NC to EM
  proj_map_cont$x0_CMU_EM <- -CST_GAMMA_LIB_EM * (proj_map_cont$x0_CMU_NCEM-CST_C1_LIB_EM)
  proj_map_cont$y0_CMU_EM <- -CST_GAMMA_LIB_EM * (proj_map_cont$y0_CMU_NCEM)
  proj_map_cont$z0_CMU_EM <- +CST_GAMMA_LIB_EM * (proj_map_cont$z0_CMU_NCEM)
  
  
  #=============================================================================
  # Plot
  #=============================================================================
  #Raw data
  pp_tiles_x0EM_y0EM_eP = plotdf_point(proj_map_tem, "x0_CMU_EM", "y0_CMU_EM", expression("x"[0]^"em"), expression("y"[0]^"em"), "pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 3)
  pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + scale_colour_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)
  
  #Continuation
  if(nrow(proj_map_cont) >0)
  {
    pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + geom_path(data = proj_map_cont, aes(x0_CMU_EM,  y0_CMU_EM) , color = "red", size = 1)
    
  }
  
  #Position of L2
  pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + geom_point(data = eml2, aes(x= xEM, y = yEM), size = 4) 
  
  pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + scale_x_continuous(limits = c(-1.18, -1.12)) 
  pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + scale_y_continuous(limits = c(-0.1, 0.1)) 
  
  #=============================================================================
  # Title
  #=============================================================================
  pp_tiles_x0EM_y0EM_eP   = pp_tiles_x0EM_y0EM_eP + ggtitle(paste0("t = ", toString(t0/CST_SEM_PERIOD_EM), "T"))
  
  #=============================================================================
  # Save
  #=============================================================================
  #filename = paste0(getwd(), "/ProjMap/VIDEO/", "pp_tiles_x0EM_y0EM_eP_t0_", sprintf("%2.4f", t0), '.pdf')
  filename = paste0(getwd(), "/ProjMap/VIDEO/", "pp_tiles_x0EM_y0EM_eP_t0_", sprintf("%i", index), '.png')
  ggsave(pp_tiles_x0EM_y0EM_eP, width = xSize, height = ySize, file = filename)
  
  index = index +1;
}