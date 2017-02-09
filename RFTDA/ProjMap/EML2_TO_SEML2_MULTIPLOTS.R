######################################################################################
# R script used in EML2_TO_SEML.R for the multplotting 
# of a projection map (connections between EML2 and SEML1,2)
######################################################################################
#===============================================================================
# Different plots as a function of s1 & s3
#===============================================================================
ppl     = list();
ttm_l   = list();
index = 1;
for(i in 1:length(s3_CMU_EM_vec))
{
  #Test: only some values of s3
  bool = s3_CMU_EM_vec[i] > -10 & s3_CMU_EM_vec[i] < 10 
  if(bool)
  {
    #Select the values in dataset
    ttm_l[[index]] = proj_map[which(proj_map$s3_CMU_EM == s3_CMU_EM_vec[i]),] 
    
    #Test: enough points to be plotted
    if(length(ttm_l[[index]]$s3_CMU_EM) > 2)
    {
      #------------------------------------------------
      # Plot : tiles in the s1_CMU_EM/s3_CMU_EM space
      #------------------------------------------------
      ppl[[index]]   = plotdf_tile_n(ttm_l[[index]], "t0_CMU_EM", "s1_CMU_EM", "t0", "s1 (EML2)", "pmin_dist_SEM", "Projection \ndistance", projection_lim_mid, FALSE, colorLimits = projection_color_lim)
      ppl[[index]]   = ppl[[index]] + scale_y_continuous(limits = c(-35, 35), breaks = seq(-34,34,4))  
      ppl[[index]]   = ppl[[index]] + ggtitle(paste0("s3 = ", toString(s3_CMU_EM_vec[i])))
      
      #------------------------------------------------
      #Index
      #------------------------------------------------
      index = index+1
    }
  }
}
#Actual plot in multiplot format
multiplot(plotlist = ppl, cols = floor(sqrt(length(s3_CMU_EM_vec))))
stop()


#===============================================================================
# Different plots as a function of the initial time
#===============================================================================
ppl     = list();
ttm_l   = list();
index = 1;
for(i in 1:length(t0_CMU_EM_vec))
{
  ttm_l[[i]] = proj_map[which(proj_map$t0_CMU_EM == t0_CMU_EM_vec[i]),] 
  
  if(length(ttm_l[[i]]$t0_CMU_EM) > 2)
  {
    #------------------------------------------------
    # Plot : tiles in the s1_CMU_EM/s3_CMU_EM space
    #------------------------------------------------
    ppl[[i]]   = plotdf_tile_n(ttm_l[[i]], "s1_CMU_EM", "s3_CMU_EM", "s1_CMU_EM", "s3_CMU_EM", "pmin_dist_SEM", "pmin_dist_SEM", projection_lim_mid, FALSE)
    ppl[[i]]   = ppl[[i]] + ggtitle(paste0("t = ", toString(t0_CMU_EM_vec[i]/CST_SEM_PERIOD_EM), "T"))
    
    #------------------------------------------------
    # Plot : points in the x0/y0 space
    #------------------------------------------------
    #     ppl[[i]] = plotdf_point(ttm_l[[i]], "x0_CMU_EM", "y0_CMU_EM", "x0", "y0","pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 3)
    #     ppl[[i]] = ppl[[i]] + scale_colour_gradient2("pmin_dist_SEM", space="Lab", midpoint = projection_lim_mid, mid = "white", high = muted("blue"), guide = FALSE)
    #     ppl[[i]] = ppl[[i]] + ggtitle(paste0("t = ", toString(t0_CMU_EM_vec[i]/CST_SEM_PERIOD_EM), "T"))
    
    #------------------------------------------------
    # Plot : points in the x0_CM_SEM/y0_CM_SEM space
    #------------------------------------------------
    #     ppl[[index]] = plotdf_point(ttm_l[[i]], "xf_CM_SEM", "yf_CM_SEM", "xf_CM_SEM", "yf_CM_SEM","pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 3)
    #     ppl[[index]] = ppl[[index]] + scale_colour_gradient2("pmin_dist_SEM", space="Lab", midpoint = projection_lim_mid, mid = "white", high = muted("blue"), guide = FALSE)
    #     #Add SEMLi
    #     ppl[[index]] = ppl[[index]] + geom_point(data = seml2, aes(x= xEM, y= yEM), size = 4) 
    #     #Title
    #     ppl[[index]]  = ppl[[index]] + ggtitle(paste0("t = ", toString(t0_CMU_EM_vec[i]/CST_SEM_PERIOD_EM), "T"))
    
    #------------------------------------------------
    #Index
    #------------------------------------------------
    index = index+1
  }
}
#Actual plot in multiplot format
multiplot(plotlist = ppl, cols = 6)
stop()

#===============================================================================
# Best trajectories
#===============================================================================
#Number of solutions
nosol = switch(DATA_SOURCE, "LOCAL" = 20, "FROM_SERVER" = 50)
#Number of solutions on each plot
nos = floor(nosol/5);
pptrajl     = list();
for(k in seq(1, floor(nosol/nos)))
{
  pptrajl[[k]] = ggplot();
  for(i in seq((k-1)*nos,k*nos-1))
  {
    #Select one solution
    proj_map_sol_min = proj_map_sol[which(proj_map_sol$label == i),]
    
    #Trajectory
    n = length(proj_map_sol_min$label)-2
    imapcs_min_traj = as.data.frame(proj_map_sol_min[1:n,], drop=false)
    
    #Final point
    imapcs_min_yf = as.data.frame(proj_map_sol_min[n+1,], drop=false)
    #Projection point
    imapcs_min_yp = as.data.frame(proj_map_sol_min[n+2,], drop=false)
    
    
    #Initial EML2 orbit
    pptrajl[[k]] = pptrajl[[k]] + geom_path(data = imapcs_min_traj, aes(x= x_orb_eml_SEM, y = y_orb_eml_SEM), color = "black") 
    
    #Manifold leg
    pptrajl[[k]] = pptrajl[[k]] + geom_path(data = imapcs_min_traj, aes(x= x_man_SEM, y = y_man_SEM), color = "red") 
    
    #Final SEML2 orbit
    pptrajl[[k]] = pptrajl[[k]] + geom_path(data = imapcs_min_traj, aes(x= x_orb_seml_SEM, y = y_orb_seml_SEM), color = "blue") 
    
    #Final point
    pptrajl[[k]] = pptrajl[[k]] + geom_point(data = imapcs_min_yf, aes(x= x_man_SEM, y = y_man_SEM), size = 2, color = "red")
    pptrajl[[k]] = pptrajl[[k]] + geom_segment(data = imapcs_min_yf, 
                                               aes(x = x_man_SEM, y = y_man_SEM, xend = x_man_SEM + 2e-1*px_man_SEM, yend = y_man_SEM + 2e-1*py_man_SEM), 
                                               arrow = arrow(length = unit(0.1,"cm")), color = "red")
    #Projection point
    pptrajl[[k]] = pptrajl[[k]] + geom_point(data = imapcs_min_yp, aes(x= x_man_SEM, y = y_man_SEM), size = 2, color = "blue")
    pptrajl[[k]] = pptrajl[[k]] + geom_segment(data = imapcs_min_yp, 
                                               aes(x = x_man_SEM, y = y_man_SEM, xend = x_man_SEM + 2e-1*px_man_SEM, yend = y_man_SEM + 2e-1*py_man_SEM), 
                                               arrow = arrow(length = unit(0.1,"cm")), color = "blue")
    #SEML2
    pptrajl[[k]] = pptrajl[[k]] + geom_point(data = seml2, aes(x= xEM, y = yEM), size = 4) 
    
    #Title
    pptrajl[[k]]  = pptrajl[[k]] + ggtitle(paste0("t = ", toString(imapcs_min_traj$t0_CMU_EM[1]/SEMperiod("SEM")), "T"))
    
    #Display
    pptrajl[[k]]  = pptrajl[[k]] + coord_fixed()
    pptrajl[[k]]
  }
}
multiplot(plotlist = pptrajl, cols = 5)
