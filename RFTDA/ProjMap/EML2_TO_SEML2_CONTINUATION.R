################################################################################
# R script to handle a projection map (connections between EML2 and SEML1,2)
# continuation case
#
# WARNING: EML2_TO_SEML2.R must be loaded once before.
#
# This scripts basically reads datafiles and updates the following dataframes:
#
# - proj_map_cont : solutions of a continuation procedures (s1, s3, tf...)
# - proj_map_cont_traj : the corresponding trajectories 
#
################################################################################

#===============================================================================
# Which specific time?
#===============================================================================
ratio_desired = 0.52# 0.065; # x SEML.us_em.T
time_desired  = ratio_desired*CST_SEM_PERIOD_EM;

#Family
FAMILY = "";#_fam2";

#===============================================================================
# Which type of data?
#===============================================================================
# Prefix
FILE_PREFIX_CONT      = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "cont_atf_order_", ORDER, "_dest_", LIB_POINT_SEM),
                         "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "cont_atf_order_", ORDER, "_dest_", LIB_POINT_SEM) )
FILE_PREFIX_CONT_TRAJ = switch(DATA_SOURCE, "LOCAL" = paste0(ftincppdafolder, FILE_SUBFOLDER, "cont_atf_traj_order_", ORDER, "_dest_", LIB_POINT_SEM),
                          "FROM_SERVER" = paste0(ftincppdafolder, FILE_SUBFOLDER, "cont_atf_traj_order_", ORDER, "_dest_", LIB_POINT_SEM) )

# Suffix
FILE_SUFFIX_CONT      = paste0("_t0_", ratio_desired);
FILE_SUFFIX_CONT_TRAJ = paste0("_t0_", ratio_desired);

#===============================================================================
# Get data from file
#===============================================================================

#-------------------------------------------------------------------------------
# Routine to actually get proj_map_cont from file
#-------------------------------------------------------------------------------
get_proj_map_cont <- function(FILE_PREFIX_CONT, FILE_SUFFIX_CONT, FAMILY="")
{
  #-----------------------------------------------------------------------------
  #Name of the data file
  #-----------------------------------------------------------------------------
  filename_cont = paste0(FILE_PREFIX_CONT, FILE_SUFFIX_CONT, FAMILY, ".txt");
  
  #-----------------------------------------------------------------------------
  # 3 possibilities:
  #
  #  - The datafile exists, e.g.:
  #         ".../Serv/cont_atf_order_20_dest_L2_t0_0.03.txt"
  #         ".../Serv/cont_atf_order_20_dest_L2_t0_0.03_fam2.txt"
  #
  #  - The datafile is here but split in two files, of the form:
  #         ".../Serv/cont_atf_order_20_dest_L2_t0_0.03_leg1.txt"
  #         ".../Serv/cont_atf_order_20_dest_L2_t0_0.03_leg2.txt"
  #
  #  - The datafile does not exist.
  #
  #-----------------------------------------------------------------------------
  if (file.exists(filename_cont))
  {
      #Read table
      proj_map_cont = read.table(file = filename_cont, header = TRUE)
  }else #we try the split leg1/leg2 option
  {
    #Try with leg2 notation
    filename_cont = paste0(FILE_PREFIX_CONT, FILE_SUFFIX_CONT, FAMILY, "_leg2.txt");
    if (file.exists(filename_cont))
    {
      #Read table
      proj_map_cont = read.table(file = filename_cont, header = TRUE)
      #Reverse table
      proj_map_cont = proj_map_cont[rev(rownames(proj_map_cont)),]
      
      #Try with leg1 notation
      filename_cont = paste0(FILE_PREFIX_CONT, FILE_SUFFIX_CONT, FAMILY, "_leg1.txt");
      if (file.exists(filename_cont))
      {
        #Read table
        proj_map_cont_leg1 = read.table(file = filename_cont, header = TRUE)
        #Get rid of first line, which is redundant with leg1
        proj_map_cont_leg1 = tail(proj_map_cont_leg1, nrow(proj_map_cont_leg1)-1)
        #Rbind with leg1
        proj_map_cont = rbind(proj_map_cont, proj_map_cont_leg1)
      }
      
    }else #no data
    {
      proj_map_cont = data.frame()
    }
  }
  
  #-----------------------------------------------------------------------------
  #Small postprocess (should be changed in the C++ code...)
  #-----------------------------------------------------------------------------
  proj_map_cont$thetae_NCSEM = proj_map_cont$te_NCSEM
  
  #-----------------------------------------------------------------------------
  #Rename of the rownames
  #-----------------------------------------------------------------------------
  if(nrow(proj_map_cont) > 0)
  {
    row.names(proj_map_cont) = seq(1, length(row.names(proj_map_cont)))
  }
    
  #-----------------------------------------------------------------------------
  #Return the dataframe
  #-----------------------------------------------------------------------------
  
  return(proj_map_cont)
}

#-------------------------------------------------------------------------------
# Continuation
#-------------------------------------------------------------------------------
proj_map_cont = get_proj_map_cont(FILE_PREFIX_CONT, FILE_SUFFIX_CONT, FAMILY)

#-------------------------------------------------------------------------------
# Routine to actually get proj_map_cont_traj from file
#-------------------------------------------------------------------------------
get_proj_map_cont_traj <- function(FILE_PREFIX_CONT_TRAJ, FILE_SUFFIX_CONT_TRAJ, FAMILY="")
{
  #-----------------------------------------------------------------------------
  # Column names
  #-----------------------------------------------------------------------------
  names = c("label",  "t_CMU_SEM", "x_CMS_NCSEM", "y_CMS_NCSEM", "z_CMS_NCSEM", "px_CMS_NCSEM", "py_CMS_NCSEM", "pz_CMS_NCSEM");
  
  #-----------------------------------------------------------------------------
  # Build the filename 
  #-----------------------------------------------------------------------------
  filepre_cont_traj  = paste0(FILE_PREFIX_CONT_TRAJ, FILE_SUFFIX_CONT_TRAJ, FAMILY);
  filename_cont_traj = paste0(filepre_cont_traj, ".bin");
  
  #-----------------------------------------------------------------------------
  # 3 possibilities:
  #
  #  - The datafile exists, e.g.:
  #         ".../Serv/cont_atf_order_20_dest_L2_t0_0.03.bin"
  #         ".../Serv/cont_atf_order_20_dest_L2_t0_0.03_fam2.bin"
  #
  #  - The datafile is here but split in two files, of the form:
  #         ".../Serv/cont_atf_order_20_dest_L2_t0_0.03_leg1.bin"
  #         ".../Serv/cont_atf_order_20_dest_L2_t0_0.03_leg2.bin"
  #
  #  - The datafile does not exist.
  #-----------------------------------------------------------------------------
  if (file.exists(filename_cont_traj))
  {
    #Read table
    proj_map_cont_traj = dffbinary(filename_cont_traj, 8, names)
    
  }else
  {
    #Try with leg2 notation
    filename_cont_traj = paste0(filepre_cont_traj, "_leg2.bin");
    if (file.exists(filename_cont_traj))
    {
      #Read table
      proj_map_cont_traj = dffbinary(filename_cont_traj, 8, names)
      #Reverse table
      proj_map_cont_traj = proj_map_cont_traj[rev(rownames(proj_map_cont_traj)),]
      #Change the labels
      proj_map_cont_traj$label = max(proj_map_cont_traj$label) - proj_map_cont_traj$label
      
      #Try with leg1 notation
      filename_cont_traj = paste0(filepre_cont_traj, "_leg1.bin");
      if (file.exists(filename_cont_traj))
      {
        #Read table
        proj_map_cont_leg1 = dffbinary(filename_cont_traj, 8, names)
        #Get rid of first line, which is redundant with leg1
        proj_map_cont_leg1 = tail(proj_map_cont_leg1, nrow(proj_map_cont_leg1)-1)
        #Change the labels
        proj_map_cont_leg1$label = max(proj_map_cont_traj$label) + proj_map_cont_leg1$label
        #Rbind with leg1
        proj_map_cont_traj = rbind(proj_map_cont_traj, proj_map_cont_leg1)
      }
      
    }else
    {
      proj_map_cont_traj = data.frame()
    }
  }
  return(proj_map_cont_traj)
}

#-------------------------------------------------------------------------------
# Continuation trajectories
#-------------------------------------------------------------------------------
proj_map_cont_traj = get_proj_map_cont_traj(FILE_PREFIX_CONT_TRAJ, FILE_SUFFIX_CONT_TRAJ, FAMILY)

#===============================================================================
# Postprocess for continuation
#===============================================================================
if(nrow(proj_map_cont_traj) > 0 && nrow(proj_map_cont) >0)
{
  #-----------------------------------------------------------------------------
  # For proj_map_cont_traj
  #-----------------------------------------------------------------------------
  #Small postprocess
  proj_map_cont_traj$x_CMS_SEM = -CST_GAMMA_LIB_SEM*(proj_map_cont_traj$x_CMS_NCSEM - CST_C1_LIB_SEM)
  proj_map_cont_traj$y_CMS_SEM = -CST_GAMMA_LIB_SEM*(proj_map_cont_traj$y_CMS_NCSEM - 0)
  proj_map_cont_traj$z_CMS_SEM = +CST_GAMMA_LIB_SEM*(proj_map_cont_traj$z_CMS_NCSEM - 0)
  
  #Label only
  proj_map_cont_label = proj_map_cont_traj[which(proj_map_cont_traj$t_CMU_SEM == min(proj_map_cont_traj$t_CMU_SEM)),]
  
  #-----------------------------------------------------------------------------
  # For proj_map_cont
  #-----------------------------------------------------------------------------
  proj_map_cont$x0_CMS_SEM = -CST_GAMMA_LIB_SEM*(proj_map_cont$x0_CMS_NCSEM - CST_C1_LIB_SEM)
  proj_map_cont$y0_CMS_SEM = -CST_GAMMA_LIB_SEM*(proj_map_cont$y0_CMS_NCSEM - 0)
  proj_map_cont$z0_CMS_SEM = +CST_GAMMA_LIB_SEM*(proj_map_cont$z0_CMS_NCSEM - 0)
  
  proj_map_cont$t0_CMU_EMT = proj_map_cont$t0_CMU_EM/CST_SEM_PERIOD_EM
  
  #Copy label into proj_map_cont
  proj_map_cont$label = proj_map_cont_label$label
}

#===============================================================================
# Plots (requires to load the EML2_TO_SEML2_POSTPROCESS.R script before)
#===============================================================================
source("ProjMap/EML2_TO_SEML2_POSTPROCESS.R")
source("ProjMap/EML2_TO_SEML2_PLOTS.R")

