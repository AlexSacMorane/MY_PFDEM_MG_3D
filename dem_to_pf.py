# -*- encoding=utf-8 -*-

import pickle, math, os, shutil
from pathlib import Path
from scipy.ndimage import binary_dilation, label
import numpy as np
import matplotlib.pyplot as plt

# own
from tools import *
from pf_to_dem import *

# -----------------------------------------------------------------------------#

def move_phasefield(dict_user, dict_sample):
    '''
    Move phase field maps by interpolation.
    '''
    print('Updating phase field maps')

    # load data
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    displacement = dict_save['displacement']

    # tracker
    dict_user['L_displacement'].append(dict_save['displacement'])
    if 'displacement' in dict_user['L_figures']:
      plot_displacement(dict_user, dict_sample) # from tools.py

    # here eta_1 does not move 

    # loading old variables
    eta_2_map = dict_sample['eta_2_map']
    # updating phase map
    eta_2_map_new = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
    # iteration on x
    for i_x in range(len(dict_sample['x_L'])):
        x = dict_sample['x_L'][i_x]
        i_x_old = 0
        # eta 1 
        if displacement[0] < 0:
            if x-displacement[0] <= dict_sample['x_L'][-1]:
                # look for window
                while not (dict_sample['x_L'][i_x_old] <= x-displacement[0] and x-displacement[0] <= dict_sample['x_L'][i_x_old+1]):
                    i_x_old = i_x_old + 1
                # interpolate
                eta_2_map_new[i_x, :, :] = (eta_2_map[(i_x_old+1), :, :] - eta_2_map[i_x_old, :, :])/(dict_sample['x_L'][i_x_old+1] - dict_sample['x_L'][i_x_old])*\
                                           (x-displacement[0] - dict_sample['x_L'][i_x_old]) + eta_2_map[i_x_old, :, :]
        elif displacement[0] > 0:
            if dict_sample['x_L'][0] <= x-displacement[0]:
                # look for window
                while not (dict_sample['x_L'][i_x_old] <= x-displacement[0] and x-displacement[0] <= dict_sample['x_L'][i_x_old+1]):
                    i_x_old = i_x_old + 1
                # interpolate
                eta_2_map_new[i_x, :, :] = (eta_2_map[(i_x_old+1), :, :] - eta_2_map[i_x_old, :, :])/(dict_sample['x_L'][i_x_old+1] - dict_sample['x_L'][i_x_old])*\
                                           (x-displacement[0] - dict_sample['x_L'][i_x_old]) + eta_2_map[i_x_old, :, :]
        else :
            eta_2_map_new = eta_2_map
    
    # loading old variables
    eta_2_map = eta_2_map_new.copy()
    # updating phase map
    eta_2_map_new = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
    # iteration on y
    for i_y in range(len(dict_sample['y_L'])):
        y = dict_sample['y_L'][i_y]
        i_y_old = 0
        # eta 1 
        if displacement[1] < 0:
            if y-displacement[1] <= dict_sample['y_L'][-1]:
                # look for window
                while not (dict_sample['y_L'][i_y_old] <= y-displacement[1] and y-displacement[1] <= dict_sample['y_L'][i_y_old+1]):
                    i_y_old = i_y_old + 1
                # interpolate
                eta_2_map_new[:, i_y, :] = (eta_2_map[:, (i_y_old+1), :] - eta_2_map[:, i_y_old, :])/(dict_sample['y_L'][i_y_old+1] - dict_sample['y_L'][i_y_old])*\
                                           (y-displacement[1] - dict_sample['y_L'][i_y_old]) + eta_2_map[:, i_y_old, :]
        elif displacement[1] > 0:
            if dict_sample['y_L'][0] <= y-displacement[1]:
                # look for window
                while not (dict_sample['y_L'][i_y_old] <= y-displacement[1] and y-displacement[1] <= dict_sample['y_L'][i_y_old+1]):
                    i_y_old = i_y_old + 1
                # interpolate
                eta_2_map_new[:, i_y, :] = (eta_2_map[:, (i_y_old+1), :] - eta_2_map[:, i_y_old, :])/(dict_sample['y_L'][i_y_old+1] - dict_sample['y_L'][i_y_old])*\
                                           (y-displacement[1] - dict_sample['y_L'][i_y_old]) + eta_2_map[:, i_y_old, :]
        else :
            eta_2_map_new = eta_2_map

    # loading old variables
    eta_2_map = eta_2_map_new.copy()
    # updating phase map
    eta_2_map_new = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
    # iteration on z
    for i_z in range(len(dict_sample['z_L'])):
        z = dict_sample['z_L'][i_z]
        i_z_old = 0
        # eta 1 
        if displacement[2] < 0:
            if z-displacement[2] <= dict_sample['z_L'][-1]:
                # look for window
                while not (dict_sample['z_L'][i_z_old] <= z-displacement[2] and z-displacement[2] <= dict_sample['z_L'][i_z_old+1]):
                    i_z_old = i_z_old + 1
                # interpolate
                eta_2_map_new[:, :, i_z] = (eta_2_map[:, :, (i_z_old+1)] - eta_2_map[:, :, i_z_old])/(dict_sample['z_L'][i_z_old+1] - dict_sample['z_L'][i_z_old])*\
                                           (z-displacement[2] - dict_sample['z_L'][i_z_old]) + eta_2_map[:, :, i_z_old]
        elif displacement[2] > 0:
            if dict_sample['z_L'][0] <= z-displacement[2]:
                # look for window
                while not (dict_sample['z_L'][i_z_old] <= z-displacement[2] and z-displacement[2] <= dict_sample['z_L'][i_z_old+1]):
                    i_z_old = i_z_old + 1
                # interpolate
                eta_2_map_new[:, :, i_z] = (eta_2_map[:, :, (i_z_old+1)] - eta_2_map[:, :, i_z_old])/(dict_sample['z_L'][i_z_old+1] - dict_sample['z_L'][i_z_old])*\
                                           (z-displacement[2] - dict_sample['z_L'][i_z_old]) + eta_2_map[:, :, i_z_old]
        else :
            eta_2_map_new = eta_2_map

    # The solute map is updated
    # the solute is push out/in of the grain
    # this is done in compute_kc() from dem_to_pf.py called later

    # update variables
    dict_sample['eta_2_map'] = eta_2_map_new

# -----------------------------------------------------------------------------#

def compute_contact(dict_user, dict_sample):
  '''
  Compute the contact characteristics:
    - box
    - maximum surface
    - volume
  '''
  # box initialization
  x_box_min = None
  x_box_max = None
  y_box_min = None
  y_box_max = None
  z_box_min = None
  z_box_max = None
  # volume
  vol_contact = 0
  # surface 
  surf_contact = 0

  # iterate on mesh
  for i_x in range(len(dict_sample['x_L'])):
    for i_y in range(len(dict_sample['y_L'])):
      # initialyze surface at this z
      surf_contact_z = 0
      for i_z in range(len(dict_sample['z_L'])):
        # contact detection
        if dict_sample['eta_1_map'][i_x, i_y, i_z] > dict_user['eta_contact_box_detection'] and\
           dict_sample['eta_2_map'][i_x, i_y, i_z] > dict_user['eta_contact_box_detection']:
          # compute box dimensions
          if x_box_min == None:
            x_box_min = dict_sample['x_L'][i_x]
            x_box_max = dict_sample['x_L'][i_x]
            y_box_min = dict_sample['y_L'][i_y]
            y_box_max = dict_sample['y_L'][i_y]
            z_box_min = dict_sample['z_L'][i_z]
            z_box_max = dict_sample['z_L'][i_z]
          else :
            if dict_sample['x_L'][i_x] < x_box_min:
              x_box_min = dict_sample['x_L'][i_x]
            if x_box_max < dict_sample['x_L'][i_x]:
              x_box_max = dict_sample['x_L'][i_x]
            if dict_sample['y_L'][i_y] < y_box_min:
              y_box_min = dict_sample['y_L'][i_y]
            if y_box_max < dict_sample['y_L'][i_y]:
              y_box_max = dict_sample['y_L'][i_y]
            if dict_sample['z_L'][i_z] < z_box_min:
              z_box_min = dict_sample['z_L'][i_z]
            if z_box_max < dict_sample['z_L'][i_z]:
              z_box_max = dict_sample['z_L'][i_z]
          # compute volume contact
          vol_contact = vol_contact + (dict_sample['x_L'][1]-dict_sample['x_L'][0])*\
                                      (dict_sample['y_L'][1]-dict_sample['y_L'][0])*\
                                      (dict_sample['z_L'][1]-dict_sample['z_L'][0])
          # compute surface at this z
          surf_contact_z = surf_contact_z + (dict_sample['x_L'][1]-dict_sample['x_L'][0])*\
                                            (dict_sample['y_L'][1]-dict_sample['y_L'][0])
      # compare surface at this z with the maximum registered
      if surf_contact < surf_contact_z:
         surf_contact = surf_contact_z  
  # if no contact detected
  if x_box_min == None:
    x_box_min = 0
    x_box_max = 0
    y_box_min = 0
    y_box_max = 0
    z_box_min = 0
    z_box_max = 0
  # save
  dict_sample['contact_box'] = [x_box_min, x_box_max, y_box_min, y_box_max, z_box_min, z_box_max]
  dict_sample['vol_contact'] = vol_contact
  dict_sample['surf_contact'] = surf_contact
  dict_user['L_contact_box_x'].append(x_box_max-x_box_min)
  dict_user['L_contact_box_y'].append(y_box_max-y_box_min)
  dict_user['L_contact_box_z'].append(z_box_max-z_box_min)
  dict_user['L_contact_volume'].append(vol_contact)
  dict_user['L_contact_surface'].append(surf_contact)

# -----------------------------------------------------------------------------#

def compute_as(dict_user, dict_sample):
    '''
    Compute activity of solid.
    '''
    # load data
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    normal_force = dict_save['normal_force']

    # init
    dict_sample['as_map'] = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
    
    # iterate on mesh
    for i_x in range(len(dict_sample['x_L'])):
        for i_y in range(len(dict_sample['y_L'])):
            for i_z in range(len(dict_sample['z_L'])):
                # contact detection
                if dict_sample['eta_1_map'][i_x, i_y, i_z] > dict_user['eta_contact_box_detection'] and\
                   dict_sample['eta_2_map'][i_x, i_y, i_z] > dict_user['eta_contact_box_detection']:
                    # determine pressure
                    P = normal_force/dict_sample['surf_contact'] # Pa
                else :
                    # determine pressure
                    P = 0 # Pa
                # save in the map
                dict_sample['as_map'][i_x, i_y, i_z] = math.exp(P*dict_user['V_m']/(dict_user['R_cst']*dict_user['temperature']))
    # save
    dict_user['L_contact_as'].append(math.exp(normal_force/dict_sample['surf_contact']*dict_user['V_m']/(dict_user['R_cst']*dict_user['temperature'])))
    dict_user['L_contact_pressure'].append(normal_force/dict_sample['surf_contact'])

    # plot 
    plot_as_pressure(dict_user, dict_sample) # from tools.py

    # write as
    write_array_txt(dict_sample, 'as', dict_sample['as_map'])

#-------------------------------------------------------------------------------

def compute_kc(dict_user, dict_sample):
    '''
    Compute the diffusion coefficient of the solute.
    Then write a .txt file needed for MOOSE simulation.

    This .txt file represent the diffusion coefficient map.
    '''
    # compute
    kc_map = np.array(np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z'])), dtype = bool)
    kc_pore_map =  np.array(np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z'])), dtype = bool)
    # iterate on x, y and z
    for i_z in range(len(dict_sample['z_L'])):
        for i_y in range(len(dict_sample['y_L'])):
            for i_x in range(len(dict_sample['x_L'])):
                if dict_sample['eta_1_map'][i_x, i_y, i_z] < dict_user['eta_contact_box_detection'] and\
                   dict_sample['eta_2_map'][i_x, i_y, i_z] < dict_user['eta_contact_box_detection']: # out of the grain
                    kc_map[i_x, i_y, i_z] = True
                    kc_pore_map[i_x, i_y, i_z] = True
                elif dict_sample['eta_1_map'][i_x, i_y, i_z] > dict_user['eta_contact_box_detection'] and\
                     dict_sample['eta_2_map'][i_x, i_y, i_z] > dict_user['eta_contact_box_detection']: # in the contact
                    kc_map[i_x, i_y, i_z] = True
                    kc_pore_map[i_x, i_y, i_z] = False
                else : # in the grain
                    kc_map[i_x, i_y, i_z] = False
                    kc_pore_map[i_x, i_y, i_z] = False

    # dilation
    dilated_M = binary_dilation(kc_map, dict_user['struct_element'])

    #compute the map of the solute diffusion coefficient
    kc_map = dict_user['D_solute']*dilated_M + 99*dict_user['D_solute']*kc_pore_map
    dict_sample['kc_map'] = kc_map

    # write
    write_array_txt(dict_sample, 'kc', dict_sample['kc_map'])

    # compute the number of grain detected in kc_map
    invert_dilated_M = np.invert(dilated_M)
    labelled_image, num_features = label(invert_dilated_M)
    dict_user['L_grain_kc_map'].append(num_features)

    # plot 
    if 'n_grain_kc_map' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
        ax1.plot(dict_user['L_grain_kc_map'])
        ax1.set_title('Number of grains detected (-)',fontsize=20)
        fig.tight_layout()
        fig.savefig('plot/n_grain_detected.png')
        plt.close(fig)

    # loading old variable
    c_map = dict_sample['c_map']
    # updating solute map
    c_map_new = c_map.copy()

    # iterate on the mesh
    for i_z in range(len(dict_sample['y_L'])):
        for i_y in range(len(dict_sample['y_L'])):
            for i_x in range(len(dict_sample['x_L'])):
                if not dilated_M[i_x, i_y, i_z]: 
                    c_map_new[i_x, i_y, i_z] = dict_user['C_eq']
                    
    # save data
    dict_sample['c_map'] = c_map_new

    # write txt for the solute concentration map
    write_array_txt(dict_sample, 'c', dict_sample['c_map'])

#-------------------------------------------------------------------------------

def write_array_txt(dict_sample, namefile, data_array):
    '''
    Write a .txt file needed for MOOSE simulation.

    This .txt represents the map of a numpy array.
    '''
    file_to_write = open('data/'+namefile+'.txt','w')
    # x
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # y
    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # z
    file_to_write.write('AXIS Z\n')
    line = ''
    for z in dict_sample['z_L']:
        line = line + str(z)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # data
    file_to_write.write('DATA\n')
    for k in range(len(dict_sample['z_L'])):
        for j in range(len(dict_sample['y_L'])):
            for i in range(len(dict_sample['x_L'])):
                file_to_write.write(str(data_array[i,j,k])+'\n')
    # close
    file_to_write.close()

#-------------------------------------------------------------------------------

def write_i(dict_user, dict_sample):
  '''
  Create the .i file to run MOOSE simulation.

  The file is generated from a template nammed PF_ACS_base.i
  '''
  file_to_write = open('pf.i','w')
  file_to_read = open('pf_base.i','r')
  lines = file_to_read.readlines()
  file_to_read.close()

  j = 0
  for line in lines :
    j = j + 1
    if j == 4:
      line = line[:-1] + ' ' + str(len(dict_sample['x_L'])-1)+'\n'
    elif j == 5:
      line = line[:-1] + ' ' + str(len(dict_sample['y_L'])-1)+'\n'
    elif j == 6:
      line = line[:-1] + ' ' + str(len(dict_sample['z_L'])-1)+'\n'
    elif j == 7:
      line = line[:-1] + ' ' + str(min(dict_sample['x_L']))+'\n'
    elif j == 8:
      line = line[:-1] + ' ' + str(max(dict_sample['x_L']))+'\n'
    elif j == 9:
      line = line[:-1] + ' ' + str(min(dict_sample['y_L']))+'\n'
    elif j == 10:
      line = line[:-1] + ' ' + str(max(dict_sample['y_L']))+'\n'
    elif j == 11:
      line = line[:-1] + ' ' + str(min(dict_sample['z_L']))+'\n'
    elif j == 12:
      line = line[:-1] + ' ' + str(max(dict_sample['z_L']))+'\n'
    elif j == 88 or j==94:
      line = line[:-1] + ' ' + str(1/dict_user['V_m'])+'\n'
    elif j == 108:
      line = line[:-1] + "'"+str(dict_user['Mobility_eff'])+' '+str(dict_user['kappa_eta'])+" 1'\n"
    elif j == 130:
      line = line[:-1] + ' ' + str(dict_user['Energy_barrier'])+"'\n"
    elif j == 143:
      line = line[:-1] + "'" + str(dict_user['C_eq']) + ' ' + str(dict_user['k_diss']) + ' ' + str(dict_user['k_prec']) + "'\n"
    elif j == 204 or j == 205 or j == 207 or j == 208:
      line = line[:-1] + ' ' + str(dict_user['crit_res']) +'\n'
    elif j == 211:
      line = line[:-1] + ' ' + str(dict_user['dt_PF']*dict_user['n_t_PF']) +'\n'
    elif j == 215:
      line = line[:-1] + ' ' + str(dict_user['dt_PF']) +'\n'
    file_to_write.write(line)

  file_to_write.close()
    
#-------------------------------------------------------------------------------

def sort_dem_files(dict_user, dict_sample):
  '''
  Sort the files from the YADE simulation.
  '''
  # rename files
  os.rename('vtk/ite_PFDEM_'+str(dict_sample['i_DEMPF_ite'])+'_lsBodies.0.vtm',\
            'vtk/ite_PFDEM_'+str(dict_sample['i_DEMPF_ite'])+'.vtm')
