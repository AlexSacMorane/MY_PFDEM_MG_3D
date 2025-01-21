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
    L_displacement = dict_save['L_displacement']

    # tracker
    dict_user['L_L_displacement'].append(dict_save['L_displacement'])
    if 'displacement' in dict_user['L_figures']:
        plot_displacement(dict_user, dict_sample) # from tools.py

    # iterate on grains
    for i_grain in range(0, len(dict_sample['L_etai_map'])):
        # read displacement
        displacement = L_displacement[i_grain]
        # print
        print('grain', i_grain, ':', displacement)

        # TRANSLATION on x
        # loading old variables
        eta_i_map = dict_sample['L_etai_map'][i_grain]
        # updating phase map
        eta_i_map_new = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
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
                    eta_i_map_new[i_x, :, :] = (eta_i_map[(i_x_old+1), :, :] - eta_i_map[i_x_old, :, :])/(dict_sample['x_L'][i_x_old+1] - dict_sample['x_L'][i_x_old])*\
                                              (x-displacement[0] - dict_sample['x_L'][i_x_old]) + eta_i_map[i_x_old, :, :]
            elif displacement[0] > 0:
                if dict_sample['x_L'][0] <= x-displacement[0]:
                    # look for window
                    while not (dict_sample['x_L'][i_x_old] <= x-displacement[0] and x-displacement[0] <= dict_sample['x_L'][i_x_old+1]):
                        i_x_old = i_x_old + 1
                    # interpolate
                    eta_i_map_new[i_x, :, :] = (eta_i_map[(i_x_old+1), :, :] - eta_i_map[i_x_old, :, :])/(dict_sample['x_L'][i_x_old+1] - dict_sample['x_L'][i_x_old])*\
                                              (x-displacement[0] - dict_sample['x_L'][i_x_old]) + eta_i_map[i_x_old, :, :]
            else :
                eta_i_map_new = eta_i_map

        # TRANSLATION on y
        # loading old variables
        eta_i_map = eta_i_map_new.copy()
        # updating phase map
        eta_i_map_new = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
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
                    eta_i_map_new[:, i_y, :] = (eta_i_map[:, (i_y_old+1), :] - eta_i_map[:, i_y_old, :])/(dict_sample['y_L'][i_y_old+1] - dict_sample['y_L'][i_y_old])*\
                                              (y-displacement[1] - dict_sample['y_L'][i_y_old]) + eta_i_map[:, i_y_old, :]
            elif displacement[1] > 0:
                if dict_sample['y_L'][0] <= y-displacement[1]:
                    # look for window
                    while not (dict_sample['y_L'][i_y_old] <= y-displacement[1] and y-displacement[1] <= dict_sample['y_L'][i_y_old+1]):
                        i_y_old = i_y_old + 1
                    # interpolate
                    eta_i_map_new[:, i_y, :] = (eta_i_map[:, (i_y_old+1), :] - eta_i_map[:, i_y_old, :])/(dict_sample['y_L'][i_y_old+1] - dict_sample['y_L'][i_y_old])*\
                                              (y-displacement[1] - dict_sample['y_L'][i_y_old]) + eta_i_map[:, i_y_old, :]
            else :
                eta_i_map_new = eta_i_map

        # TRANSLATION on z
        # loading old variables
        eta_i_map = eta_i_map_new.copy()
        # updating phase map
        eta_i_map_new = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
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
                    eta_i_map_new[:, :, i_z] = (eta_i_map[:, :, (i_z_old+1)] - eta_i_map[:, :, i_z_old])/(dict_sample['z_L'][i_z_old+1] - dict_sample['z_L'][i_z_old])*\
                                              (z-displacement[2] - dict_sample['z_L'][i_z_old]) + eta_i_map[:, :, i_z_old]
            elif displacement[2] > 0:
                if dict_sample['z_L'][0] <= z-displacement[2]:
                    # look for window
                    while not (dict_sample['z_L'][i_z_old] <= z-displacement[2] and z-displacement[2] <= dict_sample['z_L'][i_z_old+1]):
                        i_z_old = i_z_old + 1
                    # interpolate
                    eta_i_map_new[:, :, i_z] = (eta_i_map[:, :, (i_z_old+1)] - eta_i_map[:, :, i_z_old])/(dict_sample['z_L'][i_z_old+1] - dict_sample['z_L'][i_z_old])*\
                                              (z-displacement[2] - dict_sample['z_L'][i_z_old]) + eta_i_map[:, :, i_z_old]
            else :
                eta_i_map_new = eta_i_map

        # ROTATION
        if displacement[3] != 0:
            # compute center of mass
            center_x = 0
            center_y = 0
            center_z = 0
            counter = 0
            # iterate on x
            for i_x in range(len(dict_sample['x_L'])):
                # iterate on y
                for i_y in range(len(dict_sample['y_L'])):
                    # iterate on z
                    for i_z in range(len(dict_sample['z_L'])):
                        # criterion to verify the point is inside the grain
                        if dict_user['eta_contact_box_detection'] < eta_i_map_new[i_x, i_y, i_z]:
                            center_x = center_x + dict_sample['x_L'][i_x]
                            center_y = center_y + dict_sample['y_L'][i_y]
                            center_z = center_z + dict_sample['z_L'][i_z]
                            counter = counter + 1
            # compute the center of mass
            center_x = center_x/counter
            center_y = center_y/counter
            center_z = center_z/counter
            center = np.array([center_x, center_y, center_z])

            # compute matrice of rotation
            # cf french wikipedia "quaternions et rotation dans l'espace"
            a = displacement[3]
            b = displacement[4]
            c = displacement[5]
            d = displacement[6]
            M_rot = np.array([[a*a+b*b-c*c-d*d,     2*b*c-2*a*d,     2*a*c+2*b*d],
                              [    2*a*d+2*b*c, a*a-b*b+c*c-d*d,     2*c*d-2*a*b],
                              [    2*b*d-2*a*c,     2*a*b+2*c*d, a*a-b*b-c*c+d*d]])
            M_rot_inv = np.linalg.inv(M_rot)

            # loading old variables
            eta_i_map = eta_i_map_new.copy()
            # updating phase map
            eta_i_map_new = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
            # iteration on x
            for i_x in range(len(dict_sample['x_L'])):
                # iteration on y
                for i_y in range(len(dict_sample['y_L'])):
                    # iteration on z
                    for i_z in range(len(dict_sample['z_L'])):
                        # create vector of the node
                        p = np.array([dict_sample['x_L'][i_x], dict_sample['y_L'][i_y], dict_sample['z_L'][i_z]])
                        # remove the center of the grain
                        pp = p - center
                        # applied the invert rotation
                        pp = np.dot(M_rot_inv, pp)
                        # applied center
                        pp = pp + center
                        # initialization
                        found = True
                        # look for the vector in the x-axis
                        if dict_sample['x_L'][0] <= pp[0] and pp[0] <= dict_sample['x_L'][-1]:
                            i_x_old = 0
                            while not (dict_sample['x_L'][i_x_old] <= pp[0] and pp[0] <= dict_sample['x_L'][i_x_old+1]):
                                i_x_old = i_x_old + 1
                        else :
                            found = False
                        # look for the vector in the y-axis
                        if dict_sample['y_L'][0] <= pp[1] and pp[1] <= dict_sample['y_L'][-1]:
                            i_y_old = 0
                            while not (dict_sample['y_L'][i_y_old] <= pp[1] and pp[1] <= dict_sample['y_L'][i_y_old+1]):
                                i_y_old = i_y_old + 1
                        else :
                            found = False
                        # look for the vector in the z-axis
                        if dict_sample['z_L'][0] <= pp[2] and pp[2] <= dict_sample['z_L'][-1]:
                            i_z_old = 0
                            while not (dict_sample['z_L'][i_z_old] <= pp[2] and pp[2] <= dict_sample['z_L'][i_z_old+1]):
                                i_z_old = i_z_old + 1
                        else :
                            found = False
                        # triple interpolation if point found
                        if found :
                            # points
                            p000 = np.array([  dict_sample['x_L'][i_x_old],   dict_sample['y_L'][i_y_old],   dict_sample['z_L'][i_z_old]])
                            p100 = np.array([dict_sample['x_L'][i_x_old+1],   dict_sample['y_L'][i_y_old],   dict_sample['z_L'][i_z_old]])
                            p010 = np.array([  dict_sample['x_L'][i_x_old], dict_sample['y_L'][i_y_old+1],   dict_sample['z_L'][i_z_old]])
                            p001 = np.array([  dict_sample['x_L'][i_x_old],   dict_sample['y_L'][i_y_old], dict_sample['z_L'][i_z_old+1]])
                            p101 = np.array([dict_sample['x_L'][i_x_old+1],   dict_sample['y_L'][i_y_old], dict_sample['z_L'][i_z_old+1]])
                            p011 = np.array([  dict_sample['x_L'][i_x_old], dict_sample['y_L'][i_y_old+1], dict_sample['z_L'][i_z_old+1]])
                            p110 = np.array([dict_sample['x_L'][i_x_old+1], dict_sample['y_L'][i_y_old+1],   dict_sample['z_L'][i_z_old]])
                            p111 = np.array([dict_sample['x_L'][i_x_old+1], dict_sample['y_L'][i_y_old+1], dict_sample['z_L'][i_z_old+1]])
                            # values
                            q000 = eta_i_map[  i_x_old,   i_y_old,   i_z_old]
                            q100 = eta_i_map[i_x_old+1,   i_y_old,   i_z_old]
                            q010 = eta_i_map[  i_x_old, i_y_old+1,   i_z_old]
                            q001 = eta_i_map[  i_x_old,   i_y_old, i_z_old+1]
                            q101 = eta_i_map[i_x_old+1,   i_y_old, i_z_old+1]
                            q011 = eta_i_map[  i_x_old, i_y_old+1, i_z_old+1]
                            q110 = eta_i_map[i_x_old+1, i_y_old+1,   i_z_old]
                            q111 = eta_i_map[i_x_old+1, i_y_old+1, i_z_old+1]

                            # interpolation following the z-axis
                            # points
                            p00 = np.array([  dict_sample['x_L'][i_x_old],   dict_sample['y_L'][i_y_old]])
                            p10 = np.array([dict_sample['x_L'][i_x_old+1],   dict_sample['y_L'][i_y_old]])
                            p01 = np.array([  dict_sample['x_L'][i_x_old], dict_sample['y_L'][i_y_old+1]])
                            p11 = np.array([dict_sample['x_L'][i_x_old+1], dict_sample['y_L'][i_y_old+1]])
                            # values
                            q00 = (q000*(p001[2]-pp[2]) + q001*(pp[2]-p000[2]))/(p001[2]-p000[2])
                            q10 = (q100*(p101[2]-pp[2]) + q101*(pp[2]-p100[2]))/(p101[2]-p100[2])
                            q01 = (q010*(p011[2]-pp[2]) + q011*(pp[2]-p010[2]))/(p011[2]-p010[2])
                            q11 = (q110*(p111[2]-pp[2]) + q111*(pp[2]-p110[2]))/(p111[2]-p110[2])

                            # interpolation following the y-axis
                            # points
                            p0 = np.array([  dict_sample['x_L'][i_x_old]])
                            p1 = np.array([dict_sample['x_L'][i_x_old+1]])
                            # values
                            q0 = (q00*(p01[1]-pp[1]) + q01*(pp[1]-p00[1]))/(p01[1]-p00[1])
                            q1 = (q10*(p11[1]-pp[1]) + q11*(pp[1]-p10[1]))/(p11[1]-p10[1])

                            # interpolation following the x-axis
                            eta_i_map_new[i_x, i_y, i_z] = (q0*(p1[0]-pp[0]) + q1*(pp[0]-p0[0]))/(p1[0]-p0[0])

                        else :
                            eta_i_map_new[i_x, i_y, i_z] = 0

        # update variables
        dict_sample['L_etai_map'][i_grain] = eta_i_map_new

# -----------------------------------------------------------------------------#

def compute_contact(dict_user, dict_sample):
  '''
  Compute the contact characteristics:
    - box
    - maximum surface
    - volume
  '''
  # load data
  with open('data/dem_to_main.data', 'rb') as handle:
      dict_save = pickle.load(handle)
  L_contact = dict_save['L_contact']
  # initialization
  dict_sample['L_contact_box'] = []
  dict_sample['L_vol_contact'] = []
  dict_sample['L_surf_contact'] = []
  # iterate on contacts
  for contact in L_contact:
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
          if dict_sample['L_etai_map'][contact[0]][i_x, i_y, i_z] > dict_user['eta_contact_box_detection'] and\
             dict_sample['L_etai_map'][contact[1]][i_x, i_y, i_z] > dict_user['eta_contact_box_detection']:
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
    dict_sample['L_contact_box'].append([x_box_min, x_box_max, y_box_min, y_box_max, z_box_min, z_box_max])
    dict_sample['L_vol_contact'].append(vol_contact)
    dict_sample['L_surf_contact'].append(surf_contact)

  # save
  # iterate on potential contact
  ij = 0
  for i_grain in range(len(dict_sample['L_etai_map'])-1):
      for j_grain in range(i_grain+1, len(dict_sample['L_etai_map'])):
          i_contact = 0
          contact_found = L_contact[i_contact][0:2] == [i_grain, j_grain]
          while not contact_found and i_contact < len(L_contact)-1:
                i_contact = i_contact + 1
                contact_found = L_contact[i_contact][0:2] == [i_grain, j_grain]
          if dict_sample['i_DEMPF_ite'] == 1:
              if contact_found:
                  dict_user['L_L_contact_box_x'].append([dict_sample['L_contact_box'][i_contact][1]-dict_sample['L_contact_box'][i_contact][0]])
                  dict_user['L_L_contact_box_y'].append([dict_sample['L_contact_box'][i_contact][3]-dict_sample['L_contact_box'][i_contact][2]])
                  dict_user['L_L_contact_box_z'].append([dict_sample['L_contact_box'][i_contact][5]-dict_sample['L_contact_box'][i_contact][4]])
                  dict_user['L_L_contact_volume'].append([dict_sample['L_vol_contact'][i_contact]])
                  dict_user['L_L_contact_surface'].append([dict_sample['L_surf_contact'][i_contact]])
              else:
                  dict_user['L_L_contact_box_x'].append([0])
                  dict_user['L_L_contact_box_y'].append([0])
                  dict_user['L_L_contact_box_z'].append([0])
                  dict_user['L_L_contact_volume'].append([0])
                  dict_user['L_L_contact_surface'].append([0])
          else :
              if contact_found:
                  dict_user['L_L_contact_box_x'][ij].append(dict_sample['L_contact_box'][i_contact][1]-dict_sample['L_contact_box'][i_contact][0])
                  dict_user['L_L_contact_box_y'][ij].append(dict_sample['L_contact_box'][i_contact][3]-dict_sample['L_contact_box'][i_contact][2])
                  dict_user['L_L_contact_box_z'][ij].append(dict_sample['L_contact_box'][i_contact][5]-dict_sample['L_contact_box'][i_contact][4])
                  dict_user['L_L_contact_volume'][ij].append(dict_sample['L_vol_contact'][i_contact])
                  dict_user['L_L_contact_surface'][ij].append(dict_sample['L_surf_contact'][i_contact])
              else:
                  dict_user['L_L_contact_box_x'][ij].append(0)
                  dict_user['L_L_contact_box_y'][ij].append(0)
                  dict_user['L_L_contact_box_z'][ij].append(0)
                  dict_user['L_L_contact_volume'][ij].append(0)
                  dict_user['L_L_contact_surface'][ij].append(0)
          ij = ij + 1

# -----------------------------------------------------------------------------#

def compute_as(dict_user, dict_sample):
    '''
    Compute activity of solid.
    '''
    # load data
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    L_contact = dict_save['L_contact']

    # init
    dict_sample['as_map'] = np.ones((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
    L_pressure_tempo = []
    L_as_tempo = []
    # iterate on mesh
    for i_x in range(len(dict_sample['x_L'])):
        for i_y in range(len(dict_sample['y_L'])):
            for i_z in range(len(dict_sample['z_L'])):
                p_saved = False
                for i_contact in range(len(L_contact)):
                  contact = L_contact[i_contact]
                  # contact detection
                  if dict_sample['L_etai_map'][contact[0]][i_x, i_y, i_z] > dict_user['eta_contact_box_detection'] and\
                     dict_sample['L_etai_map'][contact[1]][i_x, i_y, i_z] > dict_user['eta_contact_box_detection']:
                      # determine pressure
                      P = contact[3][2]/dict_sample['L_surf_contact'][i_contact] # Pa
                      # tempo save
                      if not p_saved :
                        L_pressure_tempo.append(P)
                        L_as_tempo.append(math.exp(P*dict_user['V_m']/(dict_user['R_cst']*dict_user['temperature'])))
                        p_saved = True
                      # save in the map
                      # do not erase data
                      if dict_sample['as_map'][i_x, i_y, i_z] == 1:
                          dict_sample['as_map'][i_x, i_y, i_z] = math.exp(P*dict_user['V_m']/(dict_user['R_cst']*dict_user['temperature']))
    
    # save
    # iterate on potential contact
    ij = 0
    for i_grain in range(len(dict_sample['L_etai_map'])-1):
        for j_grain in range(i_grain+1, len(dict_sample['L_etai_map'])):
            i_contact = 0
            contact_found = L_contact[i_contact][0:2] == [i_grain, j_grain]
            while not contact_found and i_contact < len(L_contact)-1:
                i_contact = i_contact + 1
                contact_found = L_contact[i_contact][0:2] == [i_grain, j_grain]
            if dict_sample['i_DEMPF_ite'] == 1:
                if contact_found:
                    dict_user['L_L_contact_pressure'].append([L_pressure_tempo[i_contact]])
                    dict_user['L_L_contact_as'].append([L_as_tempo[i_contact]])
                else:
                    dict_user['L_L_contact_pressure'].append([0])
                    dict_user['L_L_contact_as'].append([1])
            else :
                if contact_found:
                    dict_user['L_L_contact_pressure'][ij].append(L_pressure_tempo[i_contact])
                    dict_user['L_L_contact_as'][ij].append(L_as_tempo[i_contact])
                else:
                    dict_user['L_L_contact_pressure'][ij].append(0)
                    dict_user['L_L_contact_as'][ij].append(1)
            ij = ij + 1

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
                # count the number of eta > eta_criterion
                c_eta_crit = 0
                for i_grain in range(len(dict_sample['L_etai_map'])):
                   if dict_sample['L_etai_map'][i_grain][i_x, i_y, i_z] > dict_user['eta_contact_box_detection']:
                      c_eta_crit = c_eta_crit + 1
                # compute coefficient of diffusion
                if c_eta_crit == 0: # out of the grain
                    kc_map[i_x, i_y, i_z] = True
                    kc_pore_map[i_x, i_y, i_z] = True
                elif c_eta_crit >= 2: # in the contact
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

    # HERE MUST BE MODIFIED
    # Move the solute to connserve the mass

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
    elif j == 17:
      line = ''
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line + '\t[./eta'+str(i_grain+1)+']\n'+\
                      '\t\torder = FIRST\n'+\
                      '\t\tfamily = LAGRANGE\n'+\
                      '\t\t[./InitialCondition]\n'+\
                      '\t\t\ttype = FunctionIC\n'+\
                      '\t\t\tfunction = eta'+str(i_grain+1)+'_txt\n'+\
                      '\t\t[../]\n'+\
                      '\t[../]\n'
    elif j == 27:
      line = ''
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line + '\t# Order parameter eta'+str(i_grain+1)+'\n'+\
                      '\t[./deta'+str(i_grain+1)+'dt]\n'+\
                      '\t\ttype = TimeDerivative\n'+\
                      '\t\tvariable = eta'+str(i_grain+1)+'\n'+\
                      '\t[../]\n'+\
                      '\t[./ACBulk'+str(i_grain+1)+']\n'+\
                      '\t\ttype = AllenCahn\n'+\
                      '\t\tvariable = eta'+str(i_grain+1)+'\n'+\
                      '\t\tmob_name = L\n'+\
                      '\t\tf_name = F_total\n'+\
                      '\t[../]\n'+\
                      '\t[./ACInterface'+str(i_grain+1)+']\n'+\
                      '\t\ttype = ACInterface\n'+\
                      '\t\tvariable = eta'+str(i_grain+1)+'\n'+\
                      '\t\tmob_name = L\n'+\
                      "\t\tkappa_name = 'kappa_eta'\n"+\
                      '\t[../]\n'
    elif j == 33:
      line = ''
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line + '\t[./eta'+str(i_grain+1)+'_c]\n'+\
                      '\t\ttype = CoefCoupledTimeDerivative\n'+\
                      "\t\tv = 'eta"+str(i_grain+1)+"'\n"+\
                      '\t\tvariable = c\n'+\
                      '\t\tcoef = '+str(1/dict_user['V_m'])+'\n'+\
                      '\t[../]\n'
    elif j == 46:
      line = line[:-1] + "'"+str(dict_user['Mobility_eff'])+' '+str(dict_user['kappa_eta'])+" 1'\n"
    elif j == 66:
      line = line[:-1] + "'"
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line +'eta'+str(i_grain+1)+' '
      line = line[:-1] + "'\n"
    elif j == 68:
      line = line[:-1] + "'" + str(dict_user['Energy_barrier'])+"'\n"
    elif j == 69:
      line = line[:-1] + "'"
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line +'h*(eta'+str(i_grain+1)+'^2*(1-eta'+str(i_grain+1)+')^2)+'
      line = line[:-1] + "'\n"
    elif j == 78:
      line = line[:-1] + "'"
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line +'eta'+str(i_grain+1)+' '
      line = line + "c'\n"
    elif j == 81:
      line = line[:-1] + "'" + str(dict_user['C_eq']) + ' ' + str(dict_user['k_diss']) + ' ' + str(dict_user['k_prec']) + "'\n"
    elif j == 82:
      line = line[:-1] + "'if(c<c_eq*as,k_diss,k_prec)*as*(1-c/(c_eq*as))*("
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line +'3*eta'+str(i_grain+1)+'^2-2*eta'+str(i_grain+1)+'^3+'
      line = line[:-1] +")'\n"
    elif j == 91:
      line = line[:-1] + "'"
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line +'eta'+str(i_grain+1)+' '
      line = line + "c'\n"
    elif j == 92:
      line = line[:-1] + "'F("
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line +'eta'+str(i_grain+1)+','
      line = line[:-1] + ") Ed("
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line +'eta'+str(i_grain+1)+','
      line = line + "c)'\n"
    elif j == 101:
      line = ''
      for i_grain in range(len(dict_sample['L_etai_map'])):
        line = line + '\t[eta'+str(i_grain+1)+'_txt]\n'+\
                      '\t\ttype = PiecewiseMultilinear\n'+\
                      "\t\tdata_file = data/eta_"+str(i_grain+1)+".txt\n"+\
                      '\t[]\n'
    elif j == 135 or j == 136 or j == 138 or j == 139:
      line = line[:-1] + ' ' + str(dict_user['crit_res']) +'\n'
    elif j == 142:
      line = line[:-1] + ' ' + str(dict_user['dt_PF']*dict_user['n_t_PF']) +'\n'
    elif j == 146:
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
