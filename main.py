# -*- encoding=utf-8 -*-

import math, os, errno, pickle, time, shutil
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# own
from dem_to_pf import *
from pf_to_dem import *
from ic import *
from tools import *
from Parameters import *

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Functions

def run_moose(dict_user, dict_sample):
    '''
    Prepare and run moose simulation.
    '''
    # from dem to pf
    # compute mass
    compute_mass(dict_user, dict_sample) # from tools.py
    tic_dem_to_pf = time.perf_counter() # compute dem_to_pf performances
    tic_tempo = time.perf_counter() # compute performances
    move_phasefield(dict_user, dict_sample) # in dem_to_pf.py
    tac_tempo = time.perf_counter() # compute performances
    dict_user['move_pf'] = dict_user['move_pf'] + tac_tempo-tic_tempo 
    # check mass
    compute_mass_loss(dict_user, dict_sample, 'L_loss_move_pf') # from tools.py
    
    # write pf file
    write_array_txt(dict_sample, 'eta_1', dict_sample['eta_1_map']) # from dem_to_pf.py
    write_array_txt(dict_sample, 'eta_2', dict_sample['eta_2_map']) # from dem_to_pf.py

    # from dem to pf
    # compute mass
    compute_mass(dict_user, dict_sample) # from tools.py
    # compute diffusivity map
    tic_tempo = time.perf_counter() # compute performances
    compute_kc(dict_user, dict_sample) # in dem_to_pf.py
    tac_tempo = time.perf_counter() # compute performances
    dict_user['comp_kc'] = dict_user['comp_kc'] + tac_tempo-tic_tempo 
    # check mass
    compute_mass_loss(dict_user, dict_sample, 'L_loss_kc') # from tools.py
    # compute activity map
    tic_tempo = time.perf_counter() # compute performances
    compute_as(dict_user, dict_sample) # in dem_to_pf.py
    tac_tempo = time.perf_counter() # compute performances
    dict_user['comp_as'] = dict_user['comp_as'] + tac_tempo-tic_tempo 
    
    # generate .i file
    tic_tempo = time.perf_counter() # compute performances
    write_i(dict_user, dict_sample) # in dem_to_pf.py
    tac_tempo = time.perf_counter() # compute performances
    tac_dem_to_pf = time.perf_counter() # compute dem_to_pf performances
    dict_user['L_t_dem_to_pf'].append(tac_dem_to_pf-tic_dem_to_pf)
    dict_user['write_i'] = dict_user['write_i'] + tac_tempo-tic_tempo 
    
    # pf
    # compute mass
    compute_mass(dict_user, dict_sample) # from tools.py
    print('Running PF')
    tic_pf = time.perf_counter() # compute pf performances
    os.system('mpiexec -n '+str(dict_user['n_proc'])+' ~/projects/moose/modules/phase_field/phase_field-opt -i pf.i')
    tac_pf = time.perf_counter() # compute pf performances
    dict_user['L_t_pf'].append(tac_pf-tic_pf)
    dict_user['solve_pf'] = dict_user['solve_pf'] + tac_pf-tic_pf 
    
    # from pf to dem
    tic_tempo = time.perf_counter() # compute performances
    last_j_str = sort_files(dict_user, dict_sample) # in pf_to_dem.py
    tac_dem_to_pf = time.perf_counter() # compute dem_to_pf performances
    tac_tempo = time.perf_counter() # compute performances
    dict_user['sort_pf'] = dict_user['sort_pf'] + tac_tempo-tic_tempo 

    print('Reading data')
    tic_pf_to_dem = time.perf_counter() # compute pf_to_dem performances
    read_vtk(dict_user, dict_sample, last_j_str) # in pf_to_dem.py
    tac_pf_to_dem = time.perf_counter() # compute pf_to_dem performances
    dict_user['L_t_pf_to_dem_2'].append(tac_pf_to_dem-tic_pf_to_dem)
    dict_user['read_pf'] = dict_user['read_pf'] + tac_pf_to_dem-tic_pf_to_dem 
    # check mass
    compute_mass_loss(dict_user, dict_sample, 'L_loss_pf') # from tools.py
    
# ------------------------------------------------------------------------------------------------------------------------------------------ #

def run_yade(dict_user, dict_sample):
    '''
    Prepare and run yade simulation.
    '''
    # from pf to dem
    tic_pf_to_dem = time.perf_counter() # compute pf_to_dem performances
    compute_levelset(dict_user, dict_sample) # from pf_to_dem.py
    tac_pf_to_dem = time.perf_counter() # compute pf_to_dem performances
    dict_user['L_t_pf_to_dem_1'].append(tac_pf_to_dem-tic_pf_to_dem)
    dict_user['comp_ls'] = dict_user['comp_ls'] + tac_pf_to_dem-tic_pf_to_dem 
    
    # transmit data
    tic_tempo = time.perf_counter() # compute performances
    dict_save = {
    'radius': dict_user['radius'],
    'E': dict_user['E'],
    'Poisson': dict_user['Poisson'],
    'n_ite_max': dict_user['n_ite_max'],
    'i_DEMPF_ite': dict_sample['i_DEMPF_ite'],
    'L_pos_w': dict_user['L_pos_w']
    }
    with open('data/main_to_dem.data', 'wb') as handle:
        pickle.dump(dict_save, handle, protocol=pickle.HIGHEST_PROTOCOL)
    tac_tempo = time.perf_counter() # compute performances
    dict_user['save_dem'] = dict_user['save_dem'] + tac_tempo-tic_tempo 
    
    # dem
    print('Running DEM')
    tic_dem = time.perf_counter() # compute dem performances
    os.system('yadedaily -j '+str(dict_user['n_proc'])+' -x -n dem_base.py')
    #os.system('yadedaily -j '+str(dict_user['n_proc'])+' dem_base.py')
    #os.system('~/AlexSM/myYade/install/bin/yade-2024-01-31.git-c315081 -j '+str(dict_user['n_proc'])+' -x -n dem_base.py')
    tac_dem = time.perf_counter() # compute dem performances
    dict_user['L_t_dem'].append(tac_dem-tic_dem)
    dict_user['solve_dem'] = dict_user['solve_dem'] + tac_dem-tic_dem 
    
    # load data
    tic_tempo = time.perf_counter() # compute performances
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    tac_tempo = time.perf_counter() # compute performances
    dict_user['read_dem'] = dict_user['read_dem'] + tac_tempo-tic_tempo 
        
# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Plan
    
# get parameters
dict_user = get_parameters() # from Parameters.py
dict_sample = {}

# folders
create_folder('vtk') # from tools.py
create_folder('plot') # from tools.py
if 'configuration_eta' in dict_user['L_figures'] or\
   'configuration_c' in dict_user['L_figures']:
    create_folder('plot/configuration')
create_folder('data') # from tools.py
create_folder('input') # from tools.py
create_folder('dict') # from tools.py

# compute performances
tic = time.perf_counter()

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Create initial condition

if dict_user['Shape'] == 'Sphere':
    create_1_sphere(dict_user, dict_sample) # from ic.py
create_solute(dict_user, dict_sample) # from ic.py

# compute tracker
plot_sum_mean_etai_c(dict_user, dict_sample) # from tools.py

# check if the mesh map is inside the database
check_mesh_database(dict_user, dict_sample) # from tools.py

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Performances

dict_user['move_pf'] = 0
dict_user['comp_con_vol'] = 0
dict_user['comp_kc'] = 0
dict_user['comp_as'] = 0
dict_user['write_i'] = 0
dict_user['solve_pf'] = 0
dict_user['sort_pf'] = 0
dict_user['read_pf'] = 0
dict_user['comp_ls'] = 0
dict_user['save_dem'] = 0
dict_user['solve_dem'] = 0
dict_user['read_dem'] = 0
dict_user['plot_s_m_etai_c'] = 0
dict_user['plot_perf'] = 0

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# PFDEM iteration

dict_sample['i_DEMPF_ite'] = 0
while dict_sample['i_DEMPF_ite'] < dict_user['n_DEMPF_ite']:
    dict_sample['i_DEMPF_ite'] = dict_sample['i_DEMPF_ite'] + 1
    print('\nStep',dict_sample['i_DEMPF_ite'],'/',dict_user['n_DEMPF_ite'],'\n')

    # print configuration c and eta_1
    plot_slices(dict_user, dict_sample) # from tools

    # DEM
    run_yade(dict_user, dict_sample)

    # DEM->PF, PF, PF->DEM
    run_moose(dict_user, dict_sample)

    # Evolution of sum and mean of etai + c
    tic_tempo = time.perf_counter() # compute performances
    plot_sum_mean_etai_c(dict_user, dict_sample) # from tools.py
    tac_tempo = time.perf_counter() # compute performances
    dict_user['plot_s_m_etai_c'] = dict_user['plot_s_m_etai_c'] + tac_tempo-tic_tempo 
    
    # plot performances
    tic_tempo = time.perf_counter() # compute performances
    plot_performances(dict_user, dict_sample) # from tools.py
    tac_tempo = time.perf_counter() # compute performances
    dict_user['plot_perf'] = dict_user['plot_perf'] + tac_tempo-tic_tempo 

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# close simulation

# save
with open('dict/dict_user', 'wb') as handle:
    pickle.dump(dict_user, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('dict/dict_sample', 'wb') as handle:
    pickle.dump(dict_sample, handle, protocol=pickle.HIGHEST_PROTOCOL)

# compute performances
tac = time.perf_counter()
hours = (tac-tic)//(60*60)
minutes = (tac-tic - hours*60*60)//(60)
seconds = int(tac-tic - hours*60*60 - minutes*60)
print("\nSimulation time : "+str(hours)+" hours "+str(minutes)+" minutes "+str(seconds)+" seconds")
print('Simulation ends')

# save mesh database 
save_mesh_database(dict_user, dict_sample) # from tools.py

# sort files
reduce_n_vtk_files(dict_user, dict_sample) # from tools.py

#output
print('\n')
print('mass csv:', dict_user['L_sum_mass'][0], '->', dict_user['L_sum_mass'][-1])
print('=', 100*abs(dict_user['L_sum_mass'][0]-dict_user['L_sum_mass'][-1])/dict_user['L_sum_mass'][0], '%')

