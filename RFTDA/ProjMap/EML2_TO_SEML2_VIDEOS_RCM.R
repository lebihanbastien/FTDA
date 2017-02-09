######################################################################################
# R script used in EML2_TO_SEML.R for the video
# of a projection map (connections between EML2 and SEML1,2)
######################################################################################
#===============================================================================
# Define limit of projection
#===============================================================================
#hard
projection_lim_max = 2e-3;
projection_lim_mid = 5e-4;

#===============================================================================
# Define limit of projection for colors
#===============================================================================
projection_color_lim = c(0, 1e-2);

#===============================================================================
#Keep only values with a projection distance below projection_lim_max.
#===============================================================================
proj_map = proj_map_source[which(proj_map_source$pmin_dist_SEM <= projection_lim_max),] 

#===============================================================================
# Functions for titles 
#===============================================================================
title1 <- function(t0)
{
  str = paste0("t = ", toString(t0/CST_SEM_PERIOD_EM), "T")
  return(str)
}

title2 <- function(s3)
{
  str = paste0("s3 = ", toString(s3))
  return(str)
}

title3 <- function(s1)
{
  str = paste0("s1 = ", toString(s1))
  return(str)
}

#===============================================================================
# Functions for filename 
#===============================================================================
filename1 <- function(t0, index_pdf)
{
  str = paste0(getwd(), "/ProjMap/VIDEO/", "pp_tiles_s1EM_s3EM_eP_t0_", toString(index_pdf), '.pdf')
  return(str)
}

filename2 <- function(s3, index_pdf)
{
  str = paste0(getwd(), "/ProjMap/VIDEO/", "pp_tiles_t0EM_s1EM_eP_s3_", toString(index_pdf), '.pdf')
  return(str)
}


filename3 <- function(s1, index_pdf)
{
  str = paste0(getwd(), "/ProjMap/VIDEO/", "pp_tiles_t0EM_s3EM_eP_s1_", toString(index_pdf), '.pdf')
  return(str)
}

#===============================================================================
# User-defined video type 
#===============================================================================
splash = paste0("The different possibilities of video are the following:\n",
                "1. s1_CMU_EM vs s3_CMU_EM as t0_CMU_EM increases;\n",
                "2. t0_CMU_EM vs s1_CMU_EM as s3_CMU_EM increases;\n",
                "3. t0_CMU_EM vs s3_CMU_EM as s1_CMU_EM increases;\n",
                "Enter the corresponding value (1, 2, or 3):")
video_type = readline(prompt=splash)



#===============================================================================
# Select the parameters for the rest of the computation, 
# depending on the user choice
#===============================================================================
if(video_type == 1)
{
  #=============================================================================
  # s1_CMU_EM vs s3_CMU_EM as t0_CMU_EM increases
  #=============================================================================
  # Select all unique values of time
  all_vec    = proj_map$t0_CMU_EM
  unique_vec = sort(unique(proj_map$t0_CMU_EM))
  # x and y values
  xx = "s1_CMU_EM"
  yy = "s3_CMU_EM"
  # Corresponding strings to display on plots
  xxs = expression("s"[1])
  yys = expression("s"[3])
  # Corresponding limits
  xxl = c(-35, 35)
  yyl = c(-35, 35)
  # Corresponding title function
  titlef = title1
  # Corresponding filename function
  filenamef = filename1
}else if(video_type == 2)
{
  #=============================================================================
  # t0_CMU_EM vs s1_CMU_EM as s3_CMU_EM increases
  #=============================================================================
  # Select all unique values of time
  all_vec    = proj_map$s3_CMU_EM
  unique_vec = sort(unique(proj_map$s3_CMU_EM))
  # x and y values
  xx = "t0_CMU_EM"
  yy = "s1_CMU_EM"
  # Corresponding strings to display on plots
  xxs = expression("t"[0])
  yys = expression("s"[1])
  # Corresponding limits
  xxl = c(0, 1)
  yyl = c(-35, 35)
  # Corresponding title function
  titlef = title2
  # Corresponding filename function
  filenamef = filename2
}else if(video_type == 3)
{
  #=============================================================================
  # t0_CMU_EM vs s3_CMU_EM as s1_CMU_EM
  #=============================================================================
  # Select all unique values of time
  all_vec    = proj_map$s1_CMU_EM
  unique_vec = sort(unique(proj_map$s1_CMU_EM))
  # x and y values
  xx = "t0_CMU_EM"
  yy = "s3_CMU_EM"
  # Corresponding strings to display on plots
  xxs = expression("t"[0])
  yys = expression("s"[3])
  # Corresponding limits
  xxl = c(0, 1)
  yyl = c(-35, 35)
  # Corresponding title function
  titlef = title3
  # Corresponding filename function
  filenamef = filename3
  
}else{
  #=============================================================================
  # otherwise...
  #=============================================================================
  stop("unknown video-type. Must be 1, 2 or 3")
}


#===============================================================================
# Loop
#===============================================================================
index_pdf = 0
for(zz in unique_vec)
{
  #=================================
  # Select the data
  #=================================
  proj_map_sub = proj_map[which(all_vec == zz),] 

  #=================================
  # Plot
  #=================================
  pp_tiles_eP = plotdf_tile_1(proj_map_sub, xx, yy, xxs, yys, "pmin_dist_SEM", "Projection \ndistance", FALSE, colorLimits = projection_color_lim, na.value = "white")
  pp_tiles_eP = pp_tiles_eP + scale_x_continuous(limits = xxl) 
  pp_tiles_eP = pp_tiles_eP + scale_y_continuous(limits = yyl) 
  #pp_tiles_eP =  pp_tiles_eP + theme_bw()
  
  #=================================
  # Title
  #=================================
  pp_tiles_eP   = pp_tiles_eP + ggtitle(titlef(zz))
  
  #=================================
  # Save
  #=================================
  filename = filenamef(zz, index_pdf)
  ggsave(pp_tiles_eP, width = xSize, height = ySize, file = filename)
  index_pdf = index_pdf+1
}