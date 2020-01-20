#!/usr/bin/env python

import numpy as np
import os

# Creating topo file

def topo_data(name_simulation,ictop_topo,deltax_topo,topo_x,topo_y,topo_height):

    topo_file = open('%s.top' %(name_simulation),'w')

    topo_file.write('ictop \t %s \n' %(name_simulation))
    topo_file.write('%d \n'%(ictop_topo))

    topo_file.write('np \t deltx \n')
    topo_file.write('%d \t %d \n' %(topo_x.size-1,deltax_topo))

    topo_file.write('Topo_x \t Topo_y \t Topo_zx coordinate \n')
    for i in range(0,topo_x.size-1):
        topo_file.write('%.2f \t %.2f \t %.2f \n'
                        %(topo_x[i],topo_y[i],topo_height[i]))
        
    topo_file.write('terrain \n')
    topo_file.write('%d \n' %(0))
    topo_file.close()
    os.system('echo "-----------------------------------" ')
    os.system('echo "---------- %s.top [Done]-----------" ' %
              (name_simulation))
    os.system('echo "-----------------------------------" ')
    
##############################################################################################################


    
# Creating pts file

def pts_data(name_simulation,grid_spacing_pts,deltax_pts,deltay_pts,
             facthsml_pts,source_x,source_y,source_thickness):

    source_file = open('%s.pts' %(name_simulation),'w')
    source_file.write('npoin source \t grid spacing \t deltax \t delty \t facthsml \n')
    source_file.write('%d \t %d \t %d \t %d \t %d \n'
                      %(source_x.size-1,grid_spacing_pts,deltax_pts,deltay_pts,facthsml_pts))
    source_file.write('---  X   ---------  Y   -------  h  ----- \n')
    for i in range(0,source_x.size-1):
        source_file.write('%.2f \t %.2f \t %.2f \n' 
                          %(source_x[i],source_y[i],source_thickness[i]))
    source_file.close()
    
    os.system('echo "-----------------------------------" ')
    os.system('echo "---------- %s.pts [Done]-----------" ' %
              (name_simulation))
    os.system('echo "-----------------------------------" ')
    
########################################################################


    
# Creating MASTER and master files

def master_data(name_simulation,if_sph_master,if_gfl_master,if_tgf_master,
                SPH_problem_type_master,SPH_t_integ_Alg_master,dt_master,
                time_end_master,maxtimesteps_master,print_step_master,
                t_save_step_master,t_plot_step_master,dt_sph_master,
                ic_adapt_master,Ntime_curves_master,max_pts_master):

    master = open('%s.master.dat' %(name_simulation),'w')

    # Number of comment lines
    master.write('%d \n' 
                 %(1)) 
    master.write('%s \n' 
                 %(name_simulation))

    # Boolean (1/0) Calculus modules actives
    master.write('if_sph \t if_gfl \t if_tgf \n') 
    master.write('%d \t %d \t %d \n' 
                 %(if_sph_master,if_gfl_master,if_tgf_master)) 

    # SPH parameters
    master.write('SPH_problem_type \t SPH_t_integ_Alg \n')
    master.write('%d \t %d \n' 
                 %(SPH_problem_type_master,SPH_t_integ_Alg_master))

    # SPH problem name
    master.write('sph problem name \n')
    master.write('%s \n' 
                 %(name_simulation))

    # SPH time parameters I
    master.write('dt \t time_end \t maxtimesteps \n')
    master.write('%.2f \t %.2f \t %d \n' 
                 %(dt_master,time_end_master,maxtimesteps_master))

    # Output results
    master.write('print_step \t save_step \t plot_step \n')
    master.write('%d \t %d \t %d \n' 
                 %(print_step_master,t_save_step_master,t_plot_step_master))

    # SPH time parameters II
    master.write('dt_sph \t ic_adapt \n')
    master.write('%.2f \t %d \n' 
                 %(dt_sph_master,ic_adapt_master))

    # Other parameters
    master.write('Ntime curves \t max pts in them \n')
    master.write('%d \t %d \n' 
                 %(Ntime_curves_master,max_pts_master))

    # End parameters
    master.write('dt \t time_end \t maxtimesteps \n')
    master.write('%d \t %.2f \t %d \n' 
                 %(-1,2.00,100000))

    # Close file
    master.close()
     
    os.system('echo "-----------------------------------" ')
    os.system('echo "------ %s.master.dat [Done]--------" ' %
              (name_simulation))
    os.system('echo "-----------------------------------" ')

    # Duplicate file
    os.system('cp %s.master.dat %s.MASTER.dat' 
              %(name_simulation,name_simulation))
    
    os.system('echo "-----------------------------------" ')
    os.system('echo "------ %s.MASTER.dat [Done]--------" ' %
              (name_simulation))
    os.system('echo "-----------------------------------" ')
    
#######################################################################



# Creating DAT file

def dat_data(name_simulation,ic_SW_Alg_dat,nhist_dat,ndimn_dat,
             ic_array_dat,Soil_unkno_dat,h_inf_SW_dat,pa_sph_dat,
             nnps_dat,sle_dat,skf_dat,sum_den,av_vel,virt_part,
             nor_dens,C1_C15_array_dat,C16_C30_array_dat,
             C31_C45_array_dat,C46_C60_array_dat,K0_activated_dat,
             icpwp_dat,coarse_mesh_dat,control_points,GID_filter_dat,
             T_change_to_W):
    
    dat = open('%s.dat' %(name_simulation),'w')

    # Comment lines
    dat.write('%d \n' 
              %(1))
    dat.write('---- %s \n' 
              %(name_simulation))

    # 
    dat.write('ic_SW_Alg\n')
    dat.write('%d \n' 
              %(ic_SW_Alg_dat))

    #
    dat.write('nhist \n')
    dat.write('%d \n' 
              %(nhist_dat))

    # Number of dimmensions
    dat.write('ndimn \n')
    dat.write('%d \n' 
              %(ndimn_dat))

    # Boolean variable phases
    dat.write('ic_soil \t ic_water \t ic_vps \t ic_abs \n')
    dat.write('%d \t %d \t %d \t %d \n' 
              %(ic_array_dat[0],
                ic_array_dat[1],
                ic_array_dat[2],
                ic_array_dat[3]))

    # 
    dat.write('Soil unkno \t h_inf_SW \n')
    dat.write('%d \t %.2f \n' 
              %(Soil_unkno_dat,h_inf_SW_dat))

    # Pts file name
    dat.write('pts file name \n')
    dat.write('%s \n'
              %(name_simulation)) 

    # SPH algorithm options
    dat.write('pa_sph \t nnps \t sle \t skf \n')
    dat.write('%d \t %d \t %d \t %d \n' 
              %(pa_sph_dat,nnps_dat,sle_dat,skf_dat))

    # SPH boolean algorithm 
    dat.write('sum_den \t av_vel \t virt_part \t nor_dens \n')
    dat.write('%.1s \t %.1s \t %.1s \t %.1s \n ' 
              %(sum_den,av_vel,virt_part,nor_dens))       

    # Physical parameters I
    dat.write('C1 \t C2 \t C3 \t C4 \t C5 \t C6 \t C7 \t C8 \t C9 \t C10 \t C11 \t C12 \t C13 \t C14 \t C15 \n')
    dat.write('%.2f \t %.2f \t %.2f \t %.2f \t %d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2e \t %d \t %.2f \t %.2f \t %.2f \t %.2e \n'
              %(C1_C15_array_dat[0],C1_C15_array_dat[1],C1_C15_array_dat[2],
                C1_C15_array_dat[3],C1_C15_array_dat[4],C1_C15_array_dat[5],
                C1_C15_array_dat[6],C1_C15_array_dat[7],C1_C15_array_dat[8],
                C1_C15_array_dat[9],C1_C15_array_dat[10],C1_C15_array_dat[11],
                C1_C15_array_dat[12],C1_C15_array_dat[13],C1_C15_array_dat[14]))  

    # Physical parameters II
    dat.write('C16 \t C17 \t C18 \t C19 \t C20 \t C21 \t C22 \t C23 \t C24 \t C25 \t C26 \t C27 \t C28 \t C29 \t C30 \n')
    dat.write('%.2f \t %.2f \t %.2f \t %.2f \t %d \t %d \t %.2e \t %d \t %d \t %d \t %d \t %.2e \t %d \t %d \t %d \n'
              %(C16_C30_array_dat[0],C16_C30_array_dat[1],C16_C30_array_dat[2],
                C16_C30_array_dat[3],C16_C30_array_dat[4],C16_C30_array_dat[5],
                C16_C30_array_dat[6],C16_C30_array_dat[7],C16_C30_array_dat[8],
                C16_C30_array_dat[9],C16_C30_array_dat[10],C16_C30_array_dat[11],
                C16_C30_array_dat[12],C16_C30_array_dat[13],C16_C30_array_dat[14]))    

    # Physical parameters III
    dat.write('C31 \t C32 \t C33 \t C34 \t C35 \t C36 \t C37 \t C38 \t C39 \t C40 \t C41 \t C42 \t C43 \t C44 \t C45 \n')
    dat.write('%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \n ' 
              %(C31_C45_array_dat[0],C31_C45_array_dat[1],C31_C45_array_dat[2],
                C31_C45_array_dat[3],C31_C45_array_dat[4],C31_C45_array_dat[5],
                C31_C45_array_dat[6],C31_C45_array_dat[7],C31_C45_array_dat[8],
                C31_C45_array_dat[9],C31_C45_array_dat[10],C31_C45_array_dat[11],
                C31_C45_array_dat[12],C31_C45_array_dat[13],C31_C45_array_dat[14])) 

    # Physical parameters IV
    dat.write('c46 \t c47 \t c48 \t c49 \t c50 \t c51 \t c52 \t c53 \t c54 \t c55 \t c56 \t c57 \t c58 \t c59 \t c60 \n')
    dat.write('%d \t %.2f \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \n ' 
              %(C46_C60_array_dat[0],C46_C60_array_dat[1],C46_C60_array_dat[2],
                C46_C60_array_dat[3],C46_C60_array_dat[4],C46_C60_array_dat[5],
                C46_C60_array_dat[6],C46_C60_array_dat[7],C46_C60_array_dat[8],
                C46_C60_array_dat[9],C46_C60_array_dat[10],C46_C60_array_dat[11],
                C46_C60_array_dat[12],C46_C60_array_dat[13],C46_C60_array_dat[14]))
    

    # Physical parameter V
    dat.write('K0 activated \n')
    dat.write('%d \n' 
              %(K0_activated_dat))

    # Numerical parameter I
    dat.write('icpwp \n') 
    dat.write('%d \n' 
              %(icpwp_dat))

    # Numerical parameter II
    dat.write('coarse mesh saving utility \n')
    dat.write('%d \n' 
              %(coarse_mesh_dat))

    # Numerical parameter III
    dat.write('control points \n')
    dat.write('%d \n' 
              %(control_points))

    # GID filter   
    dat.write('h_s \t disp \t vel \t  P_w \t eros \t Z \t hrel \t hw \t eta \t 10 \t hs+hw \t dumm \n')
    dat.write('%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \n'
              %(GID_filter_dat[0],GID_filter_dat[1],GID_filter_dat[2],
                GID_filter_dat[3],GID_filter_dat[4],GID_filter_dat[5],
                GID_filter_dat[6],GID_filter_dat[7],GID_filter_dat[8],
                GID_filter_dat[9],GID_filter_dat[10],GID_filter_dat[11]))

    # 
    dat.write('T_change_to_W \n')
    dat.write('%.2e \n' 
              %(T_change_to_W))

    # Close file
    dat.close()
    
    os.system('echo "-----------------------------------" ')
    os.system('echo "---------- %s.dat [Done]-----------" ' %
              (name_simulation))
    os.system('echo "-----------------------------------" ')

def autorun_gen(name_simulation):
    autorun = open('%s.sh' %(name_simulation),'w')

    autorun.write('#!/usr/bin/expect -f \n')
    autorun.write('set timeout -1 \n')
    autorun.write('spawn ./SPH_FEM_v1010208 \n')
    autorun.write('expect "Input GLOBAL problem name?" \n')
    autorun.write('send "%s \r" \n' %(name_simulation))
    autorun.write('expect "Input the problem name (TOPO mesh) ?" \n')
    autorun.write('send "%s \r" \n' %(name_simulation))
    autorun.write('expect "FORTRAN PAUSE: enter <return> or <ctrl>d to continue>" \n')
#    autorun.write('send "\r" \n')
    autorun.write('send "^d" \n')
    autorun.close()
    
    os.system('echo "-----------------------------------" ')
    os.system('echo "---------- %s.sh [Done]-----------" ' %
              (name_simulation))
    os.system('echo "-----------------------------------" ')
    
    os.system('chmod a+x %s.sh' %(name_simulation)) # Execution permission
