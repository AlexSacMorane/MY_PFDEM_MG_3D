#------------------------------------------------------------------------------------------------------------------------------------------ #
# Librairies
#------------------------------------------------------------------------------------------------------------------------------------------#

import numpy as np

#------------------------------------------------------------------------------------------------------------------------------------------ #
# Parameters
#------------------------------------------------------------------------------------------------------------------------------------------#

def get_parameters():
    '''
    Define the parameters used in the simulation.
    '''
    #---------------------------------------------------------------------#
    # Norrmalization
    n_dist = 100*1e-6 # m
    n_time = 24*60*60 # s
    n_mol = 0.73*1e3 * n_dist**3 # mol

    #---------------------------------------------------------------------#
    # PFDEM

    n_DEMPF_ite = 2 # number of PFDEM iterations
    n_proc = 4 # number of processors used
    j_total = 0 # index global of results
    n_max_vtk_files = None # maximum number of vtk files (can be None to save all files)

    # Select Figures to plot
    # Available:
    # n_grain_kc_map, sum_etai_c, configuration_eta, configuration_c, mean_etai_c, mass_loss, performances
    # dem, all_dem, overlap, normal_force, yade_vtk
    # contact_box, contact_volume, contact_surface, as, pressure
    # displacement
    L_figures = ['mean_etai_c',\
                 'overlap', \
                 'contact_volume', 'contact_surface', 'as',\
                 'displacement']

    #---------------------------------------------------------------------#
    # Grain description

    # shape of the grain
    # Sphere
    Shape = 'Sphere'

    # the radius of grains
    radius = 100*1e-6/n_dist # m/m

    #---------------------------------------------------------------------#
    # DEM (Yade)

    # steady state detection
    n_ite_max = 5000 # maximum number of iteration during a DEM step
    n_steady_state_detection = 100 # window size for the steady state detection
    steady_state_detection = 0.05 # criterion for the steady state detection

    # DEM material parameters
    # Young modulus
    E = 7e13 # (kg m-1 s-2)

    # Poisson ratio
    Poisson = 0.3

    # stiffness
    kn = E*radius
    ks = E*Poisson*radius

    # force applied
    force_applied =  0.05*E*radius**2 # (kg m s-2)

    #---------------------------------------------------------------------#
    # Phase-Field (Moose)

    # mesh
    x_min = -1.2*radius
    x_max =  1.2*radius
    y_min =  x_min
    y_max =  x_max
    n_mesh_xy = 40
    factor_z_xy = 2
    z_min = factor_z_xy*x_min
    z_max = factor_z_xy*x_max
    n_mesh_z = factor_z_xy*n_mesh_xy
    m_size_mesh = ((x_max-x_min)/(n_mesh_xy-1)+(y_max-y_min)/(n_mesh_xy-1)+(z_max-z_min)/(n_mesh_z-1))/3
    check_database = True

    # PF material parameters
    # the energy barrier
    Energy_barrier = 1
    # number of mesh in the interface
    n_int = 6
    # the interface thickness
    w = m_size_mesh*n_int
    # the gradient coefficient
    kappa_eta = Energy_barrier*w*w/9.86
    # the mobility
    Mobility_eff = 1*(100*1e-6/(24*60*60))/(n_dist/n_time) # m.s-1/(m.s-1)

    # temperature
    temperature = 623 # K 
    # molar volume
    V_m = (2.2*1e-5)/(n_dist**3/n_mol) # (m3 mol-1)/(m3 mol-1)
    # constant
    R_cst = (8.32)/(n_dist**2/(n_time**2*n_mol)) # (kg m2 s-2 mol-1 K-1)/(m2 s-2 mol-1)

    # kinetics of dissolution and precipitation
    # it affects the tilting coefficient in Ed
    k_diss = 1*(0.005)/(m_size_mesh) # ed_j = ed_i*m_i/m_j
    k_prec = k_diss # -

    # molar concentration at the equilibrium
    C_eq = (0.73*1e3)/(n_mol/n_dist**3) # (mol m-3)/(mol m-3)

    # diffusion of the solute
    size_film = m_size_mesh*5
    D_solute = (4e-14/2/size_film)/(n_dist*n_dist/n_time) # (m2 s-1)/(m2 s-1)
    n_struct_element = int(round(size_film/m_size_mesh,0))
    struct_element = np.array(np.ones((n_struct_element,n_struct_element,n_struct_element)), dtype=bool) # for dilation

    # the time stepping and duration of one PF simualtion
    dt_PF = (0.01*24*60*60)/n_time # time step
    # n_t_PF*dt_PF gives the total time duration
    #n_t_PF = 200 # number of iterations
    n_t_PF = 10 # number of iterations

    # the criteria on residual
    crit_res = 1e-3
    
    # Contact box detection
    eta_contact_box_detection = 0.1 # value of the phase field searched to determine the contact box

    #---------------------------------------------------------------------#
    # Wall positions
    # x, y, z, orientation

    L_pos_w = [[x_min, 0, 0, 0],
               [x_max, 0, 0, 0],
               [0, y_min, 0, 1],
               [0, y_max, 0, 1],
               [0, 0, z_min, 2]]

    #---------------------------------------------------------------------#
    # Grain positions
    # x, y, z

    L_pos_g = [[0, 0, -0.98*radius],
               [0, 0,  0.98*radius]]

    #---------------------------------------------------------------------#
    # trackers

    L_displacement = []
    L_overlap = []
    L_normal_force = []
    L_contact_box_x = []
    L_contact_box_y = []
    L_contact_box_z = []
    L_contact_volume = []
    L_contact_surface = []
    L_contact_as = []
    L_contact_pressure = []
    L_sum_eta_1 = []
    L_sum_eta_2 = []
    L_sum_c = []
    L_sum_mass = []
    L_m_eta_1 = []
    L_m_eta_2 = []
    L_m_c = []
    L_m_mass = []
    L_t_pf_to_dem_1 = []
    L_t_pf_to_dem_2 = []
    L_t_dem = []
    L_t_dem_to_pf = []
    L_t_pf = []
    L_grain_kc_map = []
    L_loss_move_pf_eta1 = []
    L_loss_move_pf_eta2 = []
    L_loss_move_pf_c = []
    L_loss_move_pf_m = []
    L_loss_kc_eta1 = []
    L_loss_kc_eta2 = []
    L_loss_kc_c = []
    L_loss_kc_m = []
    L_loss_pf_eta1 = []
    L_loss_pf_eta2 = []
    L_loss_pf_c = []
    L_loss_pf_m = []

    #---------------------------------------------------------------------#
    # dictionnary

    dict_user = {
    'n_dist': n_dist,
    'n_time': n_time,
    'n_mol': n_mol,
    'n_DEMPF_ite': n_DEMPF_ite,
    'n_proc': n_proc,
    'j_total': j_total,
    'n_max_vtk_files': n_max_vtk_files,
    'L_figures': L_figures,
    'n_ite_max': n_ite_max,
    'n_steady_state_detection': n_steady_state_detection,
    'steady_state_detection': steady_state_detection,
    'E': E,
    'Poisson': Poisson,
    'kn_dem': kn,
    'ks_dem': ks,
    'force_applied': force_applied,
    'Shape': Shape,
    'radius': radius,
    'check_database': check_database,
    'Energy_barrier': Energy_barrier,
    'n_int': n_int,
    'w_int': w,
    'kappa_eta': kappa_eta,
    'Mobility_eff': Mobility_eff,
    'temperature': temperature,
    'V_m': V_m,
    'R_cst': R_cst,
    'k_diss': k_diss,
    'k_prec': k_prec,
    'C_eq': C_eq,
    'size_film': size_film,
    'D_solute': D_solute,
    'struct_element': struct_element,
    'dt_PF': dt_PF,
    'n_t_PF': n_t_PF,
    'crit_res': crit_res,
    'eta_contact_box_detection': eta_contact_box_detection,
    'L_pos_g': L_pos_g,
    'L_pos_w': L_pos_w,
    'L_displacement': L_displacement,
    'L_overlap': L_overlap,
    'L_normal_force': L_normal_force,
    'L_contact_box_x': L_contact_box_x,
    'L_contact_box_y': L_contact_box_y,
    'L_contact_box_z': L_contact_box_z,
    'L_contact_volume': L_contact_volume,
    'L_contact_surface': L_contact_surface,
    'L_contact_as': L_contact_as,
    'L_contact_pressure': L_contact_pressure,
    'L_sum_eta_1': L_sum_eta_1,
    'L_sum_eta_2': L_sum_eta_2,
    'L_sum_c': L_sum_c,
    'L_sum_mass': L_sum_mass,
    'L_m_eta_1': L_m_eta_1,
    'L_m_eta_2': L_m_eta_2,
    'L_m_c': L_m_c,
    'L_m_mass': L_m_mass,
    'L_t_pf_to_dem_1': L_t_pf_to_dem_1,
    'L_t_pf_to_dem_2': L_t_pf_to_dem_2,
    'L_t_dem': L_t_dem,
    'L_t_dem_to_pf': L_t_dem_to_pf,
    'L_t_pf': L_t_pf,
    'L_grain_kc_map': L_grain_kc_map,
    'L_loss_move_pf_eta1': L_loss_move_pf_eta1,
    'L_loss_move_pf_eta2': L_loss_move_pf_eta2,
    'L_loss_move_pf_c': L_loss_move_pf_c,
    'L_loss_move_pf_m': L_loss_move_pf_m,
    'L_loss_kc_eta1': L_loss_kc_eta1,
    'L_loss_kc_eta2': L_loss_kc_eta2,
    'L_loss_kc_c': L_loss_kc_c,
    'L_loss_kc_m': L_loss_kc_m,
    'L_loss_pf_eta1': L_loss_pf_eta1,
    'L_loss_pf_eta2': L_loss_pf_eta2,
    'L_loss_pf_c': L_loss_pf_c,
    'L_loss_pf_m': L_loss_pf_m,
    'x_min': x_min,
    'x_max': x_max,
    'y_min': y_min,
    'y_max': y_max,
    'z_min': z_min,
    'z_max': z_max,
    'n_mesh_x': n_mesh_xy,
    'n_mesh_y': n_mesh_xy,
    'n_mesh_z': n_mesh_z
    }

    return dict_user
