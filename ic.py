# -*- encoding=utf-8 -*-

import numpy as np
import math, skfmm, pickle

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def create_spheres(dict_user, dict_sample):
    '''
    Create initial conditions with spheres.
    Mesh and phase field maps are generated
    '''    
    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Create initial mesh
    print("Creating initial mesh")

    x_L = np.linspace(dict_user['x_min'], dict_user['x_max'], dict_user['n_mesh_x'])
    y_L = np.linspace(dict_user['y_min'], dict_user['y_max'], dict_user['n_mesh_y'])
    z_L = np.linspace(dict_user['z_min'], dict_user['z_max'], dict_user['n_mesh_z'])

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # iterate on grains     
    print("Creating initial phase field maps")

    L_etai_map = []
    for i_grain in range(len(dict_user['L_pos_g'])):
        # position of the grains
        pos_i = dict_user['L_pos_g'][i_grain]

        # Create initial phase map  
        eta_i_map = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
    
        # iteration on x
        for i_x in range(len(x_L)):
            x = x_L[i_x]
            # iteration on y
            for i_y in range(len(y_L)):
                y = y_L[i_y]
                # iteration on z
                for i_z in range(len(z_L)):
                    z = z_L[i_z]

                    # distance to grains
                    d_node_to_g = np.linalg.norm(np.array([x,y,z])-np.array(pos_i))

                    # compute phase variable
                    if d_node_to_g <= dict_user['radius']-dict_user['w_int']/2 :
                        eta_i_map[i_x, i_y, i_z] = 1
                    elif dict_user['radius']-dict_user['w_int']/2 < d_node_to_g and d_node_to_g < dict_user['radius']+dict_user['w_int']/2:
                        eta_i_map[i_x, i_y, i_z] = 0.5*(1+math.cos(math.pi*(d_node_to_g-dict_user['radius']+dict_user['w_int']/2)/dict_user['w_int']))
                    elif dict_user['radius']+dict_user['w_int']/2 <= d_node_to_g :
                        eta_i_map[i_x, i_y, i_z] = 0
        
        # save
        L_etai_map.append(eta_i_map)

    # save dict
    dict_sample['L_etai_map'] = L_etai_map
    dict_sample['x_L'] = x_L
    dict_sample['y_L'] = y_L
    dict_sample['z_L'] = z_L

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def load_microstructure(dict_user, dict_sample):
    '''
    Load microstructure as initial conditions.
    Mesh and phase field maps are generated
    '''    
    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Load data
    with open('level_set_part0.data', 'rb') as handle:
        dict_save = pickle.load(handle)

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Create initial mesh
    print("Creating initial mesh")

    x_L = dict_save['L_x']
    y_L = dict_save['L_y']
    z_L = dict_save['L_z']

    # update dict
    dict_user['n_mesh_x'] = len(x_L)
    dict_user['n_mesh_y'] = len(y_L)
    dict_user['n_mesh_z'] = len(z_L)
    dict_user['w_int'] = ((x_L[1]-x_L[0])+(y_L[1]-y_L[0])+(z_L[1]-z_L[0]))/3*dict_user['n_int']

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Create walls

    dict_user['L_pos_w'] = dict_save['L_pos_w']
    
    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Load data
    print("Creating initial phase field maps")
    L_etai_map = []
    
    for i_data in [1,2]:
        
        with open('level_set_part'+str(i_data)+'.data', 'rb') as handle:
            dict_save = pickle.load(handle)

        # ------------------------------------------------------------------------------------------------------------------------------------------ #
        # iterate on grains     
        
        for i_grain in range(len(dict_save['L_sdf_i_map'])):
            # Create initial phase map  
            eta_i_map = np.zeros((len(x_L), len(y_L), len(z_L)))
        
            # iteration on x
            for i_x in range(len(x_L)):
                # iteration on y
                for i_y in range(len(y_L)):
                    # iteration on z
                    for i_z in range(len(z_L)):

                        # distance to grains
                        sdf = dict_save['L_sdf_i_map'][i_grain][i_x, i_y, i_z]

                        # compute phase variable
                        if sdf > dict_user['w_int']/2 :
                            eta_i_map[i_x, i_y, i_z] = 0
                        elif sdf < -dict_user['w_int']/2 :
                            eta_i_map[i_x, i_y, i_z] = 1
                        else :
                            eta_i_map[i_x, i_y, i_z] = 0.5*(1+math.cos(math.pi*(sdf+dict_user['w_int']/2)/dict_user['w_int']))
            
            # save
            L_etai_map.append(eta_i_map)

    # save dict
    dict_sample['L_etai_map'] = L_etai_map
    dict_sample['x_L'] = x_L
    dict_sample['y_L'] = y_L
    dict_sample['z_L'] = z_L
 
# ------------------------------------------------------------------------------------------------------------------------------------------ #

def create_solute(dict_user, dict_sample):
    '''
    Create the map of the solute distribution.
    '''
    c_map = np.zeros((dict_user['n_mesh_x'], dict_user['n_mesh_y'], dict_user['n_mesh_z']))
    for i_x in range(len(dict_sample['x_L'])):
        for i_y in range(len(dict_sample['y_L'])):
            for i_z in range(len(dict_sample['y_L'])):
                c_map[i_x, i_y, i_z] = dict_user['C_eq'] # system at the equilibrium initialy
    # save in dict
    dict_sample['c_map'] = c_map
