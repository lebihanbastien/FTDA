################################################################################
# R script to handle a projection map (connections between EML2 and SEML1,2)
# continuation case.
# 
# In this script, each continuation files is loaded, and we look for a certain
# phase (phase_pk) at the poincaré section x = cst. 
#
# WARNING: EML2_TO_SEML2.R must be loaded once before.
#
################################################################################

#===============================================================================
# Which type of data?
#===============================================================================
# Prefix
FILE_PREFIX_CONT      = paste0(ftincppdafolder, "plot/QBCP/EM/L2/Serv/", "cont_atf_order_", LOCAL_ORDER, "_dest_", LIB_POINT_SEM)
FILE_PREFIX_CONT_TRAJ = paste0(ftincppdafolder, "plot/QBCP/EM/L2/Serv/", "cont_atf_traj_order_", LOCAL_ORDER, "_dest_", LIB_POINT_SEM)

#Select phase angle
phase_pk = 0.5;#pi/2;

#===============================================================================
# Function for spline interpolation, column wise
#===============================================================================
fun <- function(dfc, theta_v, theta_c)
{
  sx = spline(theta_v, dfc, xout = theta_c)
  return(sx$y)
}


#===============================================================================
# Loop on specified time
#===============================================================================
proj_map_multi_cont = data.frame()
pk_map_multi_cont   = data.frame()
pk_map_multi_r    = data.frame()

for(ratio_desired in seq(0, 1.0, 0.005))
{
  #===========================================================================
  # Desired time
  #=========================================================================== 
  time_desired  = ratio_desired*CST_SEM_PERIOD_EM;
  
  #=========================================================================== 
  # Suffix from time_desired
  #=========================================================================== 
  FILE_SUFFIX_CONT  = paste0("_t0_", ratio_desired);
  
  #=========================================================================== 
  # Get data from file
  #=========================================================================== 
  proj_map_cont = get_proj_map_cont(FILE_PREFIX_CONT, FILE_SUFFIX_CONT)
  
  if(nrow(proj_map_cont) > 0)
  {
    #=========================================================================== 
    # Post process: additional columns
    #===========================================================================
    proj_map_cont$r0_CMU_EM = proj_map_cont$t0_CMU_EM / CST_SEM_PERIOD_EM
    proj_map_cont$phase     = proj_map_cont$thetae_NCSEM %% (2*pi)
    proj_map_cont$phase_pi  = proj_map_cont$phase/(pi)
    
    #=========================================================================== 
    # Post process: select the point closest to a certain phase angle at Pk map
    #===========================================================================
    #phase_index =  which.min(abs(proj_map_cont$phase-phase_pk))
    #pk_map_cont = proj_map_cont[which(proj_map_cont$phase == proj_map_cont$phase[phase_index]),]
    phase_eps     =  0.5;
    phase_indices = which(abs(proj_map_cont$phase - phase_pk) < phase_eps)
    pk_map_cont   = proj_map_cont[phase_indices,]
    pk_map_multi_cont = rbind(pk_map_multi_cont, pk_map_cont) 
    
    #===========================================================================
    # Parameters for the interpolation
    #===========================================================================
    rnmin = 1;
    rnmax = as.numeric(length(row.names(proj_map_cont)));
    rnl = 5;
    
    #===========================================================================
    # Loop on the possible thetae_NCSEM, such that thetae_NCSEM %% (2*pi) close to phase_pk
    #===========================================================================
    for (ic in seq(-15, 15, 1))
    {
      #=========================================================================     
      # Current theta_c, such that theta_c %% (2*pi) == phase_pk
      #=========================================================================
      theta_c  = phase_pk + 2*pi*ic;
      
      #=========================================================================
      # If there exist a minimum close to the good phase, we go on
      #=========================================================================
      argmin    = which(abs(proj_map_cont$thetae_NCSEM - theta_c) == min(abs(proj_map_cont$thetae_NCSEM - theta_c)))
      proj_map_min = proj_map_cont[argmin,]
      
      if(abs(proj_map_min$phase[1]- phase_pk) < phase_eps)
      {
        #=======================================================================
        # Loop on all the solutions found
        #=======================================================================
        for(id in seq(1, length(rownames(proj_map_cont))-1))
        {
          #=====================================================================
          # If we find an interval that contains phase_pk, we interpolate around 
          # this value to find the exact thetae_NCSEM %% (2*pi) == phase_pk
          #=====================================================================
          condition = (proj_map_cont$thetae_NCSEM[id] - theta_c)*(proj_map_cont$thetae_NCSEM[id+1] - theta_c) < 0
          condition = condition & abs(proj_map_cont$phase[id]- phase_pk) < phase_eps
          if(condition)
          {
            #We create a temp df
            pk_map_c = proj_map_cont[c(id, id+1),] 
            
            #We take the values around the current position
            rn = as.numeric(rownames(pk_map_c));
            low  = max(rnmin, rn - rnl)
            high = min(rnmax, rn + rnl)
            rnv = seq(low, high, 1)
            
            #Interpolation
            pk_map_r = proj_map_cont[rnv,]
            pk_map_line = colwise(fun)(pk_map_r, pk_map_r$thetae_NCSEM, theta_c)
            
            #Add to the data pk_map_multi_r
            pk_map_multi_r = rbind(pk_map_multi_r, pk_map_line)
          }
        }
      }
    }
    
    #=========================================================================== 
    # Add it to the gigantic data
    #=========================================================================== 
    proj_map_multi_cont = rbind(proj_map_multi_cont, proj_map_cont)
  }
}

#=============================================================================
# HEURISTIC POST-PROCESS: 
# 1. we get rid of the solutions that did not reach the Poincaré section
# 2. we get rid (for now) of the solution for which ye > 0, which should 
# correspond to solution for which the crossing has been computed at the second 
# flyby... Not cool!
#=============================================================================
# If te = 0, there is no crossing
proj_map_multi_cont = proj_map_multi_cont[which(proj_map_multi_cont$te_NCSEM != 0), ] 

# Only y < 0
proj_map_multi_cont = proj_map_multi_cont[which(proj_map_multi_cont$ye_CMS_NCSEM < 0), ] 


#=============================================================================
# PLOTS
#=============================================================================
pp_Pk_s1_s3 = plotdf_point(proj_map_multi_cont, "s1_CMU_EM", "s3_CMU_EM", "s1_CMU_EM", "s3_CMU_EM", "r0_CMU_EM", "r0", isColorFac = "false")
pp_Pk_s1_s3 = pp_Pk_s1_s3 + geom_point(data = pk_map_multi_r, aes(s1_CMU_EM, s3_CMU_EM), color = "green", size = 3)
pp_Pk_s1_s3

pp_Pk_x0_y0 = plotdf_point(proj_map_multi_cont, "x0_CMU_NCEM", "y0_CMU_NCEM", "x0_CMU_NCEM", "y0_CMU_NCEM", "r0_CMU_EM", "r0", isColorFac = "false")
pp_Pk_x0_y0

pp_Pk_x0_t = plotdf_point(proj_map_multi_cont, "phase", "x0_CMU_NCEM", "phase", "x0_CMU_NCEM", "r0_CMU_EM", "r0", isColorFac = "false",  pointSize = 2)
pp_Pk_x0_t = pp_Pk_x0_t + geom_point(data = pk_map_multi_cont, aes(phase, x0_CMU_NCEM), color = "black", size = 3)
pp_Pk_x0_t = pp_Pk_x0_t + geom_point(data = pk_map_multi_r, aes(phase, x0_CMU_NCEM), color = "green", size = 3)
pp_Pk_x0_t


pp_Pk_ye_pye = plotdf_point(proj_map_multi_cont, "ye_CMS_NCSEM", "pye_CMS_NCSEM", "ye_CMS_NCSEM", "pye_CMS_NCSEM", "r0_CMU_EM", "r0", isColorFac = "false",  pointSize = 2)
pp_Pk_ye_pye = pp_Pk_ye_pye + geom_point(data = pk_map_multi_cont, aes(ye_CMS_NCSEM, pye_CMS_NCSEM), color = "black", size = 3)
pp_Pk_ye_pye = pp_Pk_ye_pye + geom_point(data = pk_map_multi_r,    aes(ye_CMS_NCSEM, pye_CMS_NCSEM), color = "green", size = 3)
pp_Pk_ye_pye

xx = "x0_CMU_NCEM"
yy = "y0_CMU_NCEM"
pp_pts_phase_x0y0 = plotdf_point(pk_map_multi_r, xx, yy, xx, yy, "r0_CMU_EM", "r0", isColorFac = "false", pointSize = 3)
pp_pts_phase_x0y0

xx = "s1_CMU_EM"
yy = "s3_CMU_EM"
pp_pts_phase_s1s3 = plotdf_point(pk_map_multi_r, xx, yy, xx, yy, "r0_CMU_EM", "r0", isColorFac = "false", pointSize = 3)
pp_pts_phase_s1s3

xx = "ye_CMS_NCSEM"
yy = "pye_CMS_NCSEM"
pp_pts_phase_yepye = plotdf_point(pk_map_multi_r, xx, yy, xx, yy, "r0_CMU_EM", "r0", isColorFac = "false", pointSize = 3)
pp_pts_phase_yepye
