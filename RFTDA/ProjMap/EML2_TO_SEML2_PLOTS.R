######################################################################################
# R script used in EML2_TO_SEML.R for the plotting 
# of a projection map (connections between EML2 and SEML1,2)
######################################################################################
#===============================================================================
# CONTINUATION: Only some solutions
#===============================================================================
proj_map_cont_traj_one =  proj_map_cont_traj[which(proj_map_cont_traj$label %% 5 == 0),]
pp_path_sol = ggplot() + geom_path(data = proj_map_cont_traj_one, aes(x = x_CMS_SEM, y = y_CMS_SEM, colour = label, group = label), size = 1)
pp_path_sol = pp_path_sol + scale_colour_gradient(space="Lab", high = "black", low = "white", guide = FALSE)
pp_path_sol = pp_path_sol + custom_theme 
#Add SEMLi
seml2 = data.frame(x0_CMU_EM = 0.0, y0_CMU_EM = 0.0, z0_CMU_EM = 0.0);
seml2 = NCtoSYS(seml2, CST_GAMMA_LIB_SEM, CST_C1_LIB_SEM);
pp_path_sol = pp_path_sol + geom_point(data = seml2, aes(x= xEM, y = yEM), size = 4) 
#Labels
pp_path_sol = pp_path_sol + labs(x = "X (SEM)", y = "Y (SEM)")
pp_path_sol


#===============================================================================
# CONTINUATION: ALL found solutions as a "continuum"
#===============================================================================
proj_map_cont_traj_one =  proj_map_cont_traj#[which(proj_map_cont_traj$label < 2),]
pp_path_cont = ggplot() + geom_path(data = proj_map_cont_traj_one, aes(x = x_CMS_SEM, y = y_CMS_SEM, colour = label, group = factor(label)), size = 1.5)
pp_path_cont = pp_path_cont + scale_colour_gradient2("label", space="Lab", midpoint = max(proj_map_cont_traj_one$label)/2, mid = "white", high = muted("blue"))
pp_path_cont = pp_path_cont + custom_theme 
#Add SEMLi
seml2 = data.frame(x0_CMU_EM = 0.0, y0_CMU_EM = 0.0, z0_CMU_EM = 0.0);
seml2 = NCtoSYS(seml2, CST_GAMMA_LIB_SEM, CST_C1_LIB_SEM);
pp_path_cont = pp_path_cont + geom_point(data = seml2, aes(x= xEM, y = yEM), size = 4) 
#Labels
pp_path_cont = pp_path_cont + labs(x = "X (SEM)", y = "Y (SEM)")
pp_path_cont

#===============================================================================
# TILES
#===============================================================================
#------------------------------------------------
# Plot : tiles (pmin_dist_SEM) in the s1_CMU_EM/s3_CMU_EM space
#------------------------------------------------
pp_tiles_s1EM_s3EM_eP = plotdf_tile_n(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", "s1 (EML2)", "s3 (EML2)", "pmin_dist_SEM", "Projection \ndistance", projection_lim_mid, TRUE, colorLimits = projection_color_lim)
pp_tiles_s1EM_s3EM_eP = pp_tiles_s1EM_s3EM_eP + scale_x_continuous(limits = c(-35, 35), breaks = seq(-34,34,4)) 
pp_tiles_s1EM_s3EM_eP = pp_tiles_s1EM_s3EM_eP + scale_y_continuous(limits = c(-35, 35), breaks = seq(-34,34,4))  
pp_tiles_s1EM_s3EM_eP

#------------------------------------------------
# Plot : tiles (sf_CM_SEM) in the s1_CMU_EM/s3_CMU_EM space
#------------------------------------------------
pp_tiles_s1EM_s3EM_rf = plotdf_tile_n(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", "s1 (EML2)", "s3 (EML2)", "sf_CM_SEM", "Radius at final point", 0.4, TRUE)
pp_tiles_s1EM_s3EM_rf = pp_tiles_s1EM_s3EM_rf + scale_x_continuous(limits = c(-35, 35), breaks = seq(-34,34,4))
pp_tiles_s1EM_s3EM_rf = pp_tiles_s1EM_s3EM_rf + scale_y_continuous(limits = c(-35, 35), breaks = seq(-34,34,4)) 
pp_tiles_s1EM_s3EM_rf

#------------------------------------------------
# Plot : tiles (pmin_dist_SEM) in the s1_CMU_EM/s3_CMU_EM space, with continutation
#------------------------------------------------
pp_tiles_s1EM_s3EM_eP_cont = plotdf_tile_n(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", "s1 (EML2)", "s3 (EML2)", "pmin_dist_SEM", "Projection \ndistance", projection_lim_mid, TRUE, colorLimits = projection_color_lim)
pp_tiles_s1EM_s3EM_eP_cont = pp_tiles_s1EM_s3EM_eP_cont + scale_x_continuous(limits = c(-35, 35), breaks = seq(-34,34,4))
pp_tiles_s1EM_s3EM_eP_cont = pp_tiles_s1EM_s3EM_eP_cont + scale_y_continuous(limits = c(-35, 35), breaks = seq(-34,34,4)) 
pp_tiles_s1EM_s3EM_eP_cont

#Adding some specific points
s1 = c(-19.6734010192183, -17.1501764344739, 18.0368297229752, 8.33528384435544, 16.4521052628656)#, -3.276492877366309e+01, -3.043854834129248e+01, +3.218682800989279e+01, +3.435177734715551e+01, +2.287726449350953e+01, +2.092962158553479e+01)
s3 = c(1.11039067757786, 11.8945781333513, -32.5007644319156, -1.19257336959027, 15.5416792554699)#, +2.693265287540921e+01, +2.879570444016836e-01, +3.303849815659568e+01, -1.033936728130126e+01, -2.731602874014531e+01, -3.446251071019438e+01)
some_points = data.frame(s1_CMU_EM = s1, s3_CMU_EM = s3);
#pp_tiles_s1EM_s3EM_eP_cont = pp_tiles_s1EM_s3EM_eP_cont + geom_point(data = some_points, aes(s1, s3), color = "red", size = 5)
pp_tiles_s1EM_s3EM_eP_cont

#Adding some continuation results
pp_tiles_s1EM_s3EM_eP_cont = pp_tiles_s1EM_s3EM_eP_cont + geom_path(data = proj_map_cont, aes(s1_CMU_EM, s3_CMU_EM), color = "black", size = 1)
pp_tiles_s1EM_s3EM_eP_cont

#------------------------------------------------
# Plot : tiles (tof_EM) in the s1_CMU_EM/s3_CMU_EM space
#------------------------------------------------
pp_tiles_s1EM_s3EM_tofEM = plotdf_tile_n(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", "s1_CMU_EM", "s3_CMU_EM", "tof_EM", "tof_EM", 7, TRUE)
pp_tiles_s1EM_s3EM_tofEM

#------------------------------------------------
# Plot : tiles (dv_at_projection_SI) in the s1_CMU_EM/s3_CMU_EM space
#------------------------------------------------
pp_tiles_s1EM_s3EM_DV = plotdf_tile_n(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", "s1_CMU_EM", "s3_CMU_EM", "dv_at_projection_SI", "dv_at_projection_SI", 0.04, TRUE)
pp_tiles_s1EM_s3EM_DV

#===============================================================================
# POINTS: EML2 focus
#===============================================================================
#------------------------------------------------
# Plot : points (t0_CMU_EMT) in the s1_CMU_EM/s3_CMU_EM space
#------------------------------------------------
pp_pts_s1EM_s3EM_t0EM = plotdf_point(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", "s1_CMU_EM", "s3_CMU_EM", "t0_CMU_EMT", "t0_CMU_EMT", 0, pointSize = 3)
pp_pts_s1EM_s3EM_t0EM = pp_pts_s1EM_s3EM_t0EM + scale_colour_gradient2("t0_CMU_EMT", space="Lab", midpoint = 0.5, mid = "white", high = muted("blue"))
pp_pts_s1EM_s3EM_t0EM

#------------------------------------------------
# Plot : points, pmin_dist_SEM vs s1_CMU_EM
#------------------------------------------------
pp_pts_s1EM_eP = plotdf_point(proj_map_tem, "s1_CMU_EM", "pmin_dist_SEM", "s1_CMU_EM", "pmin_dist_SEM", pointSize = 3)
pp_pts_s1EM_eP

#------------------------------------------------
# Plot : points, pmin_dist_SEM vs s3_CMU_EM
#------------------------------------------------
pp_pts_s3EM_eP = plotdf_point(proj_map_tem, "s3_CMU_EM", "pmin_dist_SEM", "s3_CMU_EM", "pmin_dist_SEM", pointSize = 3)
pp_pts_s3EM_eP

#------------------------------------------------
# Plot : points, dx_CM_SEM, dy_CM_SEM... vs s3_CMU_EM
#------------------------------------------------
pp_pts_s3EM_dX = plotdf_point(proj_map_tem, "s3_CMU_EM", "dx_CM_SEM", "s3_CMU_EM", "dx_CM_SEM", pointSize = 3)
pp_pts_s3EM_dX = pp_pts_s3EM_dX + geom_point(data = proj_map_tem, aes(s3_CMU_EM, dy_CM_SEM),  color = "red", size = 3)
pp_pts_s3EM_dX = pp_pts_s3EM_dX + geom_point(data = proj_map_tem, aes(s3_CMU_EM, dpx_CM_SEM), color = "blue", size = 3)
pp_pts_s3EM_dX = pp_pts_s3EM_dX + geom_point(data = proj_map_tem, aes(s3_CMU_EM, dpy_CM_SEM), color = "green", size = 3)
pp_pts_s3EM_dX

#------------------------------------------------
# Plot : tiles (pmin_dist_SEM) in the s1_CMU_EM/s1_CM_SEM space
#------------------------------------------------
pp_pts_s1EM_s1SEM = plotdf_point(proj_map_tem, "s1_CMU_EM", "s1_CM_SEM", "s1_CMU_EM", "s1_CM_SEM", pointSize = 3)
pp_pts_s1EM_s1SEM


#------------------------------------------------
# Plot : points (pmin_dist_SEM) in the x0_CMU_NCEM/y0_CMU_NCEM space
#------------------------------------------------
pp_tiles_x0EM_y0EM_eP = plotdf_point(proj_map_tem, "x0_CMU_NCEM", "y0_CMU_NCEM", "x0_CMU_NCEM", "y0_CMU_NCEM","pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 3)
pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + scale_colour_gradient2("pmin_dist_SEM", space="Lab", midpoint = projection_lim_mid, mid = "white", high = muted("blue"))
pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CM_SEM - SNMAX) < 0.001),], aes(x0_CMU_NCEM, y0_CMU_NCEM), color = "red", size = 5)
pp_tiles_x0EM_y0EM_eP

#Adding some continuation results
pp_tiles_x0EM_y0EM_eP = pp_tiles_x0EM_y0EM_eP + geom_path(data = proj_map_cont, aes(x0_CMU_NCEM,  y0_CMU_NCEM) , color = "black", size = 1)
pp_tiles_x0EM_y0EM_eP

#===============================================================================
# POINTS: SEML2 focus
#===============================================================================

#------------------------------------------------
# Plot : points (pmin_dist_SEM) in the s1_CM_SEM/s3_CM_SEM space
#------------------------------------------------
pp_tiles_s1SEM_s3SEM_eP = plotdf_point(proj_map_tem, "s1_CM_SEM", "s3_CM_SEM", "s1_CM_SEM", "s3_CM_SEM","pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 3)
pp_tiles_s1SEM_s3SEM_eP = pp_tiles_s1SEM_s3SEM_eP +scale_colour_gradient2("pmin_dist_SEM", space="Lab", midpoint = projection_lim_mid, mid = "white", high = muted("blue"))
#Saturation of the constraint on s1_CM_SEM/s3_CM_SEM
pp_tiles_s1SEM_s3SEM_eP = pp_tiles_s1SEM_s3SEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CM_SEM - SNMAX) < 0.001),], aes(s1_CM_SEM, s3_CM_SEM), color = "red", size = 5)
#Saturation of the constraint on s1_CMU_EM/s3_CMU_EM
pp_tiles_s1SEM_s3SEM_eP = pp_tiles_s1SEM_s3SEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CMU_EM - GSIZEMAX) < 0.2),], aes(s1_CM_SEM, s3_CM_SEM), color = "green", size = 5)
pp_tiles_s1SEM_s3SEM_eP

#Adding some continuation results
pp_tiles_s1SEM_s3SEM_eP = pp_tiles_s1SEM_s3SEM_eP + geom_path(data = proj_map_cont, aes(s1_CMS_SEM, s3_CMS_SEM), color = "black", size = 1)
pp_tiles_s1SEM_s3SEM_eP



#------------------------------------------------
# Plot : points (pmin_dist_SEM) in the x0_CM_SEM/y0_CM_SEM space
#------------------------------------------------
pp_pts_x0SEM_y0SEM_eP = plotdf_point(proj_map_tem, "x0_CM_SEM", "y0_CM_SEM", "x0_CM_SEM", "y0_CM_SEM","pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 3)
pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + scale_colour_gradient2("pmin_dist_SEM", space="Lab", midpoint = projection_lim_mid, mid = "white", high = muted("blue"))
pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CM_SEM - SNMAX)    < 0.001),], aes(x0_CM_SEM, y0_CM_SEM), color = "red", size = 2)
pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$nf_CM_NCSEM - YNMAX) < 0.001),], aes(x0_CM_SEM, y0_CM_SEM), color = "magenta", size = 2)


#Add the points of arrival at SEML2
pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_point(data = proj_map_tem, aes(x= xf_CM_SEM, y = yf_CM_SEM, color = pmin_dist_SEM), size = 3)
pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CM_SEM - SNMAX) < 0.001),], aes(xf_CM_SEM, yf_CM_SEM), color = "red", size = 2)
pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CMU_EM - GSIZEMAX) < 0.2),], aes(xf_CM_SEM, yf_CM_SEM), color = "green", size = 2)
pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$nf_CM_NCSEM - YNMAX) < 0.001),], aes(xf_CM_SEM, yf_CM_SEM), color = "magenta", size = 2)


#Add SEMLi
seml2 = data.frame(x = 0.0, y = 0.0, z = 0.0);
seml2 = NCtoSYS(seml2, CST_GAMMA_LIB_SEM, CST_C1_LIB_SEM);
pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_point(data = seml2, aes(x= xEM, y = yEM), size = 4) 

#Adding some continuation results
#pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + geom_path(data = proj_map_cont, aes(x0_CMS_SEM,  y0_CMS_SEM) , color = "black", size = 1)
#pp_pts_x0SEM_y0SEM_eP


#Scaling
# pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + scale_x_continuous(breaks = seq(-1.012, -1.0025, 0.001), limits = c(-1.012, -1.002)) 
# pp_pts_x0SEM_y0SEM_eP = pp_pts_x0SEM_y0SEM_eP + scale_y_continuous(limits = c(-0.005, +0.005)) 
pp_pts_x0SEM_y0SEM_eP

#===============================================================================
# POINTS: pmin_dist_SEM vs the world
#===============================================================================
#------------------------------------------------
# Plot : points, pmin_dist_SEM vs tof_EM
#------------------------------------------------
pp_pts_tofEM_eP = plotdf_point(proj_map_tem, "tof_EM", "pmin_dist_SEM", "tof_EM", "pmin_dist_SEM", pointSize = 3)
#Saturation of the constraint on s1_CM_SEM/s3_CM_SEM
pp_pts_tofEM_eP = pp_pts_tofEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CM_SEM - SNMAX) < 0.001),], aes(tof_EM, pmin_dist_SEM), color = "red", size = 3)
#Saturation of the constraint on s1_CMU_EM/s3_CMU_EM
pp_pts_tofEM_eP = pp_pts_tofEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CMU_EM - GSIZEMAX) < 1),], aes(tof_EM, pmin_dist_SEM), color = "green", size = 3)
#Saturation of the constraint on yprojection
pp_pts_tofEM_eP = pp_pts_tofEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$nf_CM_NCSEM - YNMAX) < 0.001),], aes(tof_EM, pmin_dist_SEM), color = "blue", size = 3)
pp_pts_tofEM_eP

#------------------------------------------------
# Plot : points, pmin_dist_SEM vs t0_CMU_EMT
#------------------------------------------------
pp_pts_t0EMT_eP = plotdf_point(proj_map_tem, "t0_CMU_EMT", "pmin_dist_SEM", "t0_CMU_EM/T", "pmin_dist_SEM", pointSize = 3)
pp_pts_t0EMT_eP

#------------------------------------------------
# Plot : points, pmin_dist_SEM vs tf_man_SEMT
#------------------------------------------------
pp_pts_tfEMT_eP = plotdf_point(proj_map_tem, "tf_man_SEMT", "pmin_dist_SEM", "tf_man_SEM/T", "pmin_dist_SEM", pointSize = 3)
pp_pts_tfEMT_eP

#------------------------------------------------
# Plot : points, pmin_dist_SEM vs ns_CM_SEM
#------------------------------------------------
pp_pts_nsSEM_eP = plotdf_point(proj_map_tem, "ns_CM_SEM", "pmin_dist_SEM", pointSize = 3)
pp_pts_nsSEM_eP = pp_pts_nsSEM_eP + geom_point(data = proj_map_tem[which(abs(proj_map_tem$ns_CM_SEM - SNMAX) < 0.001),], aes(ns_CM_SEM, pmin_dist_SEM), color = "red", size = 5)
pp_pts_nsSEM_eP

#------------------------------------------------
# Plot : points, pmin_dist_SEM vs nf_CM_NCSEM
#------------------------------------------------
pp_pts_nfNCSEM_eP = plotdf_point(proj_map_tem, "nf_CM_NCSEM", "pmin_dist_SEM", "nf_CM_NCSEM", "pmin_dist_SEM", pointSize = 3)
pp_pts_nfNCSEM_eP

#===============================================================================
# TRAJECTORY
#===============================================================================
#--------------------------------
#Select one solution
#--------------------------------
proj_map_sol_min = proj_map_sol[which(proj_map_sol$label == 0),]
#proj_map_sol_min = proj_map_sol[which(proj_map_sol$pmin_dist_SEM == min(proj_map_sol$pmin_dist_SEM)),]

#--------------------------------
# Post process
#--------------------------------
#Trajectory: we get rid of the last two values that contain other data
#(probably not the best choice, but works for now!)
n = length(proj_map_sol_min$label)-2
imapcs_min_traj = as.data.frame(proj_map_sol_min[1:n,], drop=false)
#Final point
imapcs_min_yf = as.data.frame(proj_map_sol_min[n+1,], drop=false)
#Projection point
imapcs_min_yp = as.data.frame(proj_map_sol_min[n+2,], drop=false)

#--------------------------------
# Plot it
#--------------------------------
#Initial EML2 orbit
pp_path_traj = plotdf_path(imapcs_min_traj, "x_orb_eml_SEM", "y_orb_eml_SEM", "Xsem", "Ysem")

#Manifold leg
pp_path_traj = pp_path_traj + geom_path(data = imapcs_min_traj, aes(x= x_man_SEM, y = y_man_SEM), color = "red") 

#Final SEML2 orbit
pp_path_traj = pp_path_traj + geom_path(data = imapcs_min_traj, aes(x= x_orb_seml_SEM, y = y_orb_seml_SEM), color = "blue") 

#Final point
pp_path_traj = pp_path_traj + geom_point(data = imapcs_min_yf, aes(x= x_man_SEM, y = y_man_SEM), size = 2, color = "red")
pp_path_traj = pp_path_traj + geom_segment(data = imapcs_min_yf, 
                                           aes(x = x_man_SEM, y= y_man_SEM, xend = x_man_SEM + 2e-1*px_man_SEM, yend = y_man_SEM + 2e-1*py_man_SEM), 
                                           arrow = arrow(length = unit(0.1,"cm")), color = "red")
#Projection point
pp_path_traj = pp_path_traj + geom_point(data = imapcs_min_yp, aes(x= x_man_SEM, y = y_man_SEM), size = 2, color = "blue")
pp_path_traj = pp_path_traj + geom_segment(data = imapcs_min_yp, 
                                           aes(x = x_man_SEM, y = y_man_SEM, xend = x_man_SEM + 2e-1*px_man_SEM, yend = y_man_SEM + 2e-1*py_man_SEM), 
                                           arrow = arrow(length = unit(0.1,"cm")), color = "blue")
#SEML2
pp_path_traj = pp_path_traj + geom_point(data = seml2, aes(x= xEM, y = yEM), size = 4) 

#Title
pp_path_traj  = pp_path_traj + ggtitle(paste0("t = ", toString(imapcs_min_traj$t0_CMU_EM[1]/SEMperiod("SEM")), "T"))

#Display
pp_path_traj  = pp_path_traj + coord_fixed()
pp_path_traj



#--------------------------------
#Energy
#--------------------------------
pp_path_energy  = plotdf_path(imapcs_min_traj, "t_orb_eml_SEM", "h_orb_eml_SEM")
pp_path_energy  = pp_path_energy + geom_path(data = imapcs_min_traj, aes(x= t0_CMU_EM, y = h_man_SEM), color = "red") 
pp_path_energy  = pp_path_energy + geom_path(data = imapcs_min_traj, aes(x= t_orb_seml_SEM, y = h_orb_seml_SEM), color = "blue") 
pp_path_energy  = pp_path_energy + geom_point(data = imapcs_min_yp,  aes(x= t0_CMU_EM, y = h_man_SEM), size = 2, color = "blue")
pp_path_energy


#------------------------------------------------------------------------------------------------------------
#Interpolation
#------------------------------------------------------------------------------------------------------------
# #Add an arbitrary value
# # xs1 = c(proj_map_tem$s1_CMU_EM, 35, 40)
# # xs3 = c(proj_map_tem$s3_CMU_EM, -3.5, 5.5)
# # Plot : points (pmin_dist_SEM) in the s1_CMU_EM/s3_CMU_EM space
# ppePmp = plotdf_point(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", "s1_CMU_EM", "s3_CMU_EM", pointSize = 3)
# 
# 
# #Interpolation (spline)
# splineData <- data.frame(
#   spline(xs1, xs3),
#   method = "spline()"
# )
# #Corresponding function
# splineF = splinefun(xs1, xs3)
# 
# #Interpolation (linear)
# approxData <- data.frame(
#   approx(xs1, xs3),
#   method = "approx()"
# )
# #Corresponding function
# approxF = approxfun(xs1, xs3)
# 
# 
# #Plot again, with the interpolation
# ppePmp = ppePmp + geom_line(dat=splineData, aes(x0_CMU_EM, y0_CMU_EM), color = "red")
# ppePmp = ppePmp + geom_line(dat=approxData, aes(x0_CMU_EM, y0_CMU_EM), color = "blue")
# ppePmp
# 
# #------------------------------------------------------------------------------------------------------------
# #Write interpolation data in file
# #------------------------------------------------------------------------------------------------------------
# dfint = data.frame(approx(xs1, xs3))
# write.table(dfint, file = "data_interpolation.txt", sep = " ", row.names = FALSE, col.names = FALSE)
