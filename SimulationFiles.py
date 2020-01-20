#!/usr/bin/python3 python

import numpy as np
import pandas as pd
import os
import sys
import PyGrams as PyGrm


ValuesBCC = (['test/Dyka_MPM/CurveConstat.txt', 'NULL'],
             ['NULL', 'NULL'],
             ['NULL', 'NULL' ],
             ['NULL', 'NULL'])

PyGrm.Write_GramsBox('ejemplo', 'GID', 'mesh_ejemplo.msh', 'V', ValuesBCC)

PyGrm.Write_GramsSolid2D('ejemplo','MPM_Mesh.msh')

# file = open(sys.argv[1], 'r')

# Manning_value = float(sys.argv[2]) #500
# turbulence = float(sys.argv[3]) #0.5

# # Read files
# topography = pd.read_excel(file.name,
#                            sheet_name=0,
#                            header = None,
#                            skiprows = 1,
#                            dtype = {'0': float,
#                                     '1': float,
#                                     '2': float,
#                            } ,
#                            decimal=',')


# source = pd.read_excel(file.name,
#                        sheet_name=1,
#                        header = None,
#                        skiprows = 1,
#                        dtype = {'0': float,
#                                 '1': float,
#                                 '2': float,
#                        } ,
#                        decimal=',')


# os.system('echo "-----------------------------------" ')
# os.system('echo "------- Files lecture [Done]-------" ')
# os.system('echo "-----------------------------------" ')



# # Creating vectors data

# topo_x = topography[0].values
# topo_y = topography[1].values
# topo_height = topography[2].values

# source = source.sort_values(by=2,axis=1,ascending=False)
# source = source[source[2] > 0.0]
# source_x = source[0].values
# source_y = source[1].values
# source_thickness = source[2].values



# # Name of the simulation
# name_simulation = 'Goldau'


# # Remove previous data

# os.system('rm %s.top' %(name_simulation))
# os.system('rm %s.pts' %(name_simulation))
# os.system('rm %s.master.dat' %(name_simulation))
# os.system('rm %s.MASTER.dat' %(name_simulation))
# os.system('rm %s.dat' %(name_simulation))
# os.system('rm %s.sh' %(name_simulation))


# # Creating files


# # Define constants (.top)
# ictop_topo = 11
# deltax_topo = 25

# SPH.topo_data(name_simulation,ictop_topo,deltax_topo,
#               topo_x,topo_y,topo_height)


# ############################################################


# # Define constants (.pts)
# grid_spacing_pts = 25
# deltax_pts = 15
# deltay_pts = 15
# facthsml_pts = 2

# SPH.pts_data(name_simulation,grid_spacing_pts,
#              deltax_pts,deltay_pts,facthsml_pts,
#              source_x,source_y,source_thickness)


# ############################################################

# # Define constants (.MASTER.dat & .master.dat)

# if_sph_master = 1 
# if_gfl_master = 0 
# if_tgf_master = 0
# SPH_problem_type_master = 1
# SPH_t_integ_Alg_master = 4
# dt_master = 0.1
# time_end_master = 600
# maxtimesteps_master = 1000000
# print_step_master = 25
# t_save_step_master = 25
# t_plot_step_master = 25
# dt_sph_master = 0.1
# ic_adapt_master = 1
# Ntime_curves_master = 0
# max_pts_master = 6


# SPH.master_data(name_simulation,if_sph_master,if_gfl_master,if_tgf_master,
#                 SPH_problem_type_master,SPH_t_integ_Alg_master,dt_master,
#                 time_end_master,maxtimesteps_master,print_step_master,
#                 t_save_step_master,t_plot_step_master,dt_sph_master,
#                 ic_adapt_master,Ntime_curves_master,max_pts_master)




# ############################################################

# # Generate file autorun (.sh)

# SPH.autorun_gen(name_simulation)


# ############################################################

# # Define constants (.dat)

# ic_SW_Alg_dat = 0

# nhist_dat = 0

# ndimn_dat = 2

# ic_soil_dat = 1
# ic_water_dat = 0
# ic_vps_dat = 0
# ic_abs_dat = 0
# ic_array_dat = (ic_soil_dat,ic_water_dat,ic_vps_dat,ic_abs_dat)

# Soil_unkno_dat = 7
# h_inf_SW_dat = 0.1

# pa_sph_dat = 2
# nnps_dat = 2
# sle_dat = 2
# skf_dat = 1

# sum_den = True
# av_vel = True 
# virt_part = False
# nor_dens = False

# C1 = 9.81 # Gravity
# C2 = 2000 # Density
# C3 = Manning_value # Manning
# C4 = 0.
# C5 = 7 # Rheology model
# C6 = 0.0 # Tauy0
# C7 = 0.0 # constK
# C8 = 0.0 # visco
# C9 = turbulence # tanfi8 
# C10 = 0.001 # hfrict0
# C11 = -999
# C12 = turbulence
# C13 = 0.0
# C14 = 0.0
# C15 = 0.001

# C16 = 0.0
# C17 = 2000.0
# C18 = 1000.0
# C19 = 1.0
# C20 = 0
# C21 = 1
# C22 = 2e-3
# C23 = 0
# C24 = 0
# C25 = 0
# C26 = 0
# C27 = 1e8
# C28 = 1
# C29 = 0
# C30 = -999

# C31 = 0
# C32 = 0
# C33 = 0
# C34 = 0
# C35 = 0
# C36 = 0
# C37 = 0
# C38 = 0
# C39 = 0
# C40 = 0
# C41 = 0
# C42 = 1 # Solid law (0/1)
# C43 = 7000 # T_end law
# C44 = 3 # Exponential factor
# C45 = -999


# C46 = 0
# C47 = 0.50
# C48 = 0
# C49 = 0
# C50 = 0
# C51 = 0
# C52 = 0
# C53 = 0
# C54 = 0
# C55 = 0
# C56 = 0
# C57 = 0
# C58 = 0
# C59 = 0
# C60 = 0

# C1_C15_array_dat = (C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15)
# C16_C30_array_dat = (C16,C17,C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,
#                      C29,C30)
# C31_C45_array_dat = (C31,C32,C33,C34,C35,C36,C37,C38,C39,C40,C41,C42,C43,
#                      C44,C45)
# C46_C60_array_dat = (C46,C47,C48,C49,C50,C51,C52,C53,C54,C55,C56,C57,C58,
#                      C59,C60)


# K0_activated_dat = 0

# icpwp_dat = 0

# coarse_mesh_dat = 0 # coarse mesh saving utility

# control_points = 0

# GID_hs_dat = 1
# GID_disp_dat = 1
# GID_vel_dat = 0
# GID_Pw_dat = 0
# GID_eros_dat = 0
# GID_Z_dat = 0
# GID_hrel_dat = 0
# GID_hw_dat = 0
# GID_eta_dat = 0
# GID_10_dat = 0
# GID_hs_hw_dat = 0
# GID_dumm_dat = 0

# GID_filter_dat = (GID_hs_dat,GID_disp_dat,GID_vel_dat,GID_Pw_dat,
#                   GID_Pw_dat,GID_eros_dat,GID_Z_dat,GID_hrel_dat,
#                   GID_hw_dat,GID_eta_dat,GID_10_dat,GID_hs_hw_dat,
#                   GID_dumm_dat)

# T_change_to_W = 1.e+12


# # Data execution

# SPH.dat_data(name_simulation,ic_SW_Alg_dat,nhist_dat,ndimn_dat,
#              ic_array_dat,Soil_unkno_dat,h_inf_SW_dat,pa_sph_dat,
#              nnps_dat,sle_dat,skf_dat,sum_den,av_vel,virt_part,
#              nor_dens,C1_C15_array_dat,C16_C30_array_dat,
#              C31_C45_array_dat,C46_C60_array_dat,K0_activated_dat,
#              icpwp_dat,coarse_mesh_dat,control_points,GID_filter_dat,
#              T_change_to_W)

# ############################################################

# # Run simulation

# os.system('./%s.sh' %(name_simulation))


# ############################################################


# # Move results to a new folder

# os.system('mkdir Simulation_%f_%f'
#           %(Manning_value,turbulence))
# os.system('mv %s.dat ./Simulation_%f_%f'
#           %(name_simulation,Manning_value,turbulence))
# os.system('mv %s.restart.dat ./Simulation_%f_%f'
#           %(name_simulation,Manning_value,turbulence))
# os.system('mv *.res ./Simulation_%f_%f'
#           %(Manning_value,turbulence))
# os.system('mv *.chk ./Simulation_%f_%f'
#           %(Manning_value,turbulence))
# os.system('mv *.msh ./Simulation_%f_%f'
#           %(Manning_value,turbulence))


# os.system('echo "-----------------------------------" ')
# os.system('echo "----- Simulation finished ---------" ')
# os.system('echo "-----------------------------------" ')
