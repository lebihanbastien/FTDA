######################################################################################
# R script used in EML2_TO_SEML.R for the plotting 
# of a projection map (connections between EML2 and SEML1,2)
######################################################################################


#===============================================================================
# CONTINUATION: Only certain values
#===============================================================================
maxlabel = Inf
proj_map_cont_traj_in =  proj_map_cont_traj[which(proj_map_cont_traj$label < maxlabel),]
proj_map_cont_in      =  proj_map_cont[which(proj_map_cont$label < maxlabel),]

#===============================================================================
# CONTINUATION: Only one solution
#===============================================================================
proj_map_cont_traj_one =  proj_map_cont_traj_in[which(proj_map_cont_traj_in$label < 2),]
pp_path_sol = ggplot() + geom_path(data = proj_map_cont_traj_one, aes(x = x_CMS_SEM, y = y_CMS_SEM, colour = factor(label), group = factor(label)), size = 1.5)
pp_path_sol = pp_path_sol + scale_colour_brewer(palette = "Set1", direction = -1, guide = FALSE) #scale_colour_manual(values = brewer.pal(2,"Set1"), guide = FALSE)


#Add SEMLi
seml2 = data.frame(x0_CMU_NCEM = 0.0, y0_CMU_NCEM = 0.0, z0_CMU_NCEM = 0.0);
seml2 = NCtoSYS(seml2, CST_GAMMA_LIB_SEM, CST_C1_LIB_SEM);
pp_path_sol = pp_path_sol + geom_point(data = seml2, aes(x= xEM, y = yEM), size = 3) 

#Add Earth
earth = data.frame(x = -1.0, y = 0.0, z0_CMU_NCEM = 0.0);
pp_path_sol = pp_path_sol + geom_point(data = earth, aes(x= x, y = y), size = 4) 

#Labels
pp_path_sol = pp_path_sol + labs(x = "$X$", y = "$Y$")

#Limits
pp_path_sol = pp_path_sol + scale_x_continuous(limits = c(-1.013, -0.997))
pp_path_sol = pp_path_sol + scale_y_continuous(limits = c(-0.011, 0.011)) 

#Theme
pp_path_sol = pp_path_sol + custom_bw_theme
pp_path_sol = pp_path_sol+ coord_fixed(ratio=1)

#Save
pp_path_sol = pp_path_sol + labs(x = "$X$", y = "$Y$")
ggplot2tikz(pp_path_sol, width = xSize, height = ySize, file = paste0(filepre_cont_tra, "_tex.tex"))
pp_path_sol = pp_path_sol + labs(x = "X", y = "Y")
ggsave(pp_path_sol, width = xSize, height = xSize,  bg = "transparent",  file = paste0(filepre_cont_tra, ".pdf")) #Save in pdf

pp_path_sol
#stop()


#===============================================================================
# CONTINUATION: Only some solutions
#===============================================================================
proj_map_cont_traj_one =  proj_map_cont_traj_in[which(proj_map_cont_traj_in$label %% 6 == 0),]
pp_path_some = ggplot() + geom_path(data = proj_map_cont_traj_one, aes(x = x_CMS_SEM, y = y_CMS_SEM, colour = label, group = label), colour = "black", size = 0.8)
#pp_path_some = pp_path_some + scale_colour_gradient(space="Lab", high = "black", low = "white", guide = FALSE)

#Add SEMLi
seml2 = data.frame(x0_CMU_NCEM = 0.0, y0_CMU_NCEM = 0.0, z0_CMU_NCEM = 0.0);
seml2 = NCtoSYS(seml2, CST_GAMMA_LIB_SEM, CST_C1_LIB_SEM);
pp_path_some = pp_path_some + geom_point(data = seml2, aes(x= xEM, y = yEM), size = 3) 

#Add Earth
earth = data.frame(x = -1.0, y = 0.0, z0_CMU_NCEM = 0.0);
pp_path_some = pp_path_some + geom_point(data = earth, aes(x= x, y = y), size = 4) 

#Labels
pp_path_some = pp_path_some + labs(x = "$X$", y = "$Y$")

#Limits
#pp_path_some = pp_path_some + scale_x_continuous(limits = c(-1.013, -0.997))
#pp_path_some = pp_path_some + scale_y_continuous(limits = c(-0.012, 0.013)) 

#Theme
pp_path_some = pp_path_some + custom_theme
pp_path_some = pp_path_some+ coord_fixed(ratio=1)

#Save
pp_path_some = pp_path_some + labs(x = "$X$", y = "$Y$")
ggplot2tikz(pp_path_some, width = xSize, height = ySize, file = paste0(filepre_cont_tra, "_", toString(maxlabel), "_tex.tex"))
pp_path_some = pp_path_some + labs(x = "X", y = "Y")
ggsave(pp_path_some, width = xSize, height = xSize,  bg = "transparent",  file = paste0(filepre_cont_tra, "_", toString(maxlabel), ".pdf")) #Save in pdf
ggsave(pp_path_some, width = xSize, height = xSize,  bg = "transparent",  file = paste0(filepre_cont_tra, "_", toString(maxlabel), ".png")) #Save in png

pp_path_some

#===============================================================================
# TILES
#===============================================================================
#------------------------------------------------
# Plot : tiles (pmin_dist_SEM) in the s1_CMU_EM/s3_CMU_EM space, with continutation
#------------------------------------------------
ppEM = plotdf_tile_n(proj_map_tem, "s1_CMU_EM", "s3_CMU_EM", expression("s"[1] ("EML"[2])), expression("s"[3] ("EML"[2])), "pmin_dist_SEM", "Projection \ndistance", projection_lim_mid, TRUE, colorLimits = projection_color_lim)
ppEM = ppEM + scale_x_continuous(limits = c(-35, 35), breaks = seq(-36,36,6))
ppEM = ppEM + scale_y_continuous(limits = c(-35, 35), breaks = seq(-36,36,6)) 
ppEM

#Adding some specific points
s1 = c(-19.6734010192183, -17.1501764344739, 18.0368297229752, 8.33528384435544, 16.4521052628656)#, -3.276492877366309e+01, -3.043854834129248e+01, +3.218682800989279e+01, +3.435177734715551e+01, +2.287726449350953e+01, +2.092962158553479e+01)
s3 = c(1.11039067757786, 11.8945781333513, -32.5007644319156, -1.19257336959027, 15.5416792554699)#, +2.693265287540921e+01, +2.879570444016836e-01, +3.303849815659568e+01, -1.033936728130126e+01, -2.731602874014531e+01, -3.446251071019438e+01)
some_points = data.frame(s1_CMU_EM = s1, s3_CMU_EM = s3);
#ppEM = ppEM + geom_point(data = some_points, aes(s1, s3), color = "white", size = 5)
#ppEM = ppEM + geom_point(data = some_points, aes(s1, s3), color = "black", size = 3)

#Adding some continuation results
ppEM = ppEM + geom_path(data = proj_map_cont_in, aes(s1_CMU_EM, s3_CMU_EM), color = "white", size = 3)
ppEM = ppEM + geom_path(data = proj_map_cont_in, aes(s1_CMU_EM, s3_CMU_EM), color = "black", size = 1)
ppEM


#Colors
ppEM = ppEM + scale_colour_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)
ppEM = ppEM + scale_fill_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)
#Theme
ppEM = ppEM + custom_theme

#Save
ggsave(ppEM, width = xSize, height = xSize,  bg = "transparent",  file = paste0(filepre, "_", toString(maxlabel), ".png")) #Save png
ggsave(ppEM, width = xSize, height = xSize,  bg = "transparent",  file = paste0(filepre, "_", toString(maxlabel), ".pdf")) #Save in pdf

#===============================================================================
# POINTS
#===============================================================================
#------------------------------------------------
# Plot : points (pmin_dist_SEM) in the x0_CMU_NCEM/y0_CMU_NCEM space
#------------------------------------------------
ppx0NCEM = plotdf_point(proj_map_tem, "x0_CMU_NCEM", "y0_CMU_NCEM", "x0_CMU_NCEM", "y0_CMU_NCEM","pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 1)

#Colors
ppx0NCEM = ppx0NCEM + scale_colour_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)
ppx0NCEM = ppx0NCEM + scale_fill_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)
#Theme
ppx0NCEM = ppx0NCEM + custom_theme

#Add EMLi
eml2 = data.frame(x0_CMU_NCEM = 0.0, y0_CMU_NCEM = 0.0, z0_CMU_NCEM = 0.0);
ppx0NCEM = ppx0NCEM + geom_point(data = eml2, aes(x= x0_CMU_NCEM, y = y0_CMU_NCEM), size = 3) 

#Add Moon
#moon = data.frame(x = -1.0, y = 0.0, z0_CMU_NCEM = 0.0);
#ppx0NCEM = ppx0NCEM + geom_point(data = moon, aes(x= x, y = y), size = 4) 

ppx0NCEM


#------------------------------------------------
# Plot : points (pmin_dist_SEM) in the x0_CMU_EM/y0_CMU_EM space
#------------------------------------------------
ppx0EM = plotdf_point(proj_map_tem, "x0_CMU_EM", "y0_CMU_EM", "x (EM)", "y (EM)","pmin_dist_SEM", "pmin_dist_SEM", 0, pointSize = 1)

#Colors
ppx0EM = ppx0EM + scale_colour_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)
ppx0EM = ppx0EM + scale_fill_gradient("pmin_dist_SEM", space="Lab", high = "white", low = muted("blue"), limits = projection_color_lim, guide = FALSE)

#Theme
ppx0EM = ppx0EM + custom_theme
ppx0EM = ppx0EM+ coord_fixed(ratio=1)

#Add EMLi
eml2 = data.frame(x = 0.0, y = 0.0, z = 0.0);
eml2 = NCtoSYS(eml2, CST_GAMMA_LIB_EM, CST_C1_LIB_EM);
ppx0EM = ppx0EM + geom_point(data = eml2, aes(x= xEM, y = yEM), size = 3) 

#Add Moon
moon = data.frame(x = -1.0, y = 0.0, z = 0.0);
moon = NCtoSYS(moon, CST_GAMMA_LIB_EM, CST_C1_LIB_EM);
ppx0EM = ppx0EM + geom_point(data = moon, aes(x= xEM, y = yEM), size = 4) 

ppx0EM

#Save
ggsave(ppx0EM, width = xSize, height = xSize,  bg = "transparent",  file = paste0(filepre, "_EM.png")) #Save png
