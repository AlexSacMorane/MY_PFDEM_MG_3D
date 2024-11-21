# -*- encoding=utf-8 -*-

from yade import pack, utils, plot, export
import polyhedra_utils
import pickle
import numpy as np

# Own functions
# -----------------------------------------------------------------------------#

def create_materials():
    '''
    Create materials.
    '''
    O.materials.append(FrictMat(young=E, poisson=Poisson, frictionAngle=atan(0.5), density=2000, label='frictMat'))
    O.materials.append(FrictMat(young=E, poisson=Poisson, frictionAngle=0, density=2000, label='frictlessMat'))

# -----------------------------------------------------------------------------#

def create_grains():
    '''
    Recreate level set from data extrapolated with phase field output.
    '''
    print("Creating level set")

    # grid
    grid = RegularGrid(
        min=(min(L_x), min(L_y), min(L_z)),
        nGP=(len(L_x), len(L_y), len(L_z)),
        spacing=L_x[1]-L_x[0] 
    )   

    # grains
    O.bodies.append(
        levelSetBody(grid=grid,
                     distField=sdf_1_map.tolist(),
                     material=0)
    )
    O.bodies[-1].state.blockedDOFs = 'xyzXYZ'
    O.bodies.append(
        levelSetBody(grid=grid,
                     distField=sdf_2_map.tolist(),
                     material=0)
    )
    O.bodies[-1].state.blockedDOFs = 'xy'

# -----------------------------------------------------------------------------#

def create_walls():
    '''
    Recreate walls.
    '''
    for i_wall in range(len(L_pos_w)):
        O.bodies.append(wall((L_pos_w[i_wall][0], L_pos_w[i_wall][1], L_pos_w[i_wall][2]), L_pos_w[i_wall][3], material=1))

# -----------------------------------------------------------------------------#

def check():
    '''
    Try to detect a steady-state.
    A maximum number of iteration is used.
    '''
    if O.iter == n_ite_max:
        O.pause() # stop DEM simulation
    
# -----------------------------------------------------------------------------#

def compute_dt():
    '''
    Compute the time step used in the DEM step.
    '''
    O.dt = 0.2*SpherePWaveTimeStep(radius=radius, density=2000, young=E)

# -----------------------------------------------------------------------------#

def create_engines():
    '''
    Create engines.
    '''
    O.engines = [
            ForceResetter(),
            InsertionSortCollider([Bo1_LevelSet_Aabb(), Bo1_Wall_Aabb()], verletDist=0.00),
            InteractionLoop(
                    [Ig2_LevelSet_LevelSet_ScGeom(), Ig2_Wall_LevelSet_ScGeom()],
                    [Ip2_FrictMat_FrictMat_FrictPhys()],
                    [Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=False)]),
            NewtonIntegrator(damping=0.1, label='newton', gravity=(0, 0, 0)),
    		PyRunner(command='check()',iterPeriod=1, label='checker')
    ]

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Load data

# from main
with open('data/main_to_dem.data', 'rb') as handle:
    dict_save = pickle.load(handle)
radius = dict_save['radius']
E = dict_save['E']
Poisson = dict_save['Poisson']
n_ite_max = dict_save['n_ite_max']
i_DEMPF_ite = dict_save['i_DEMPF_ite']
L_pos_w = dict_save['L_pos_w']

# from plane interpolation
with open('data/level_sets.data', 'rb') as handle:
    dict_save = pickle.load(handle)
L_x = dict_save['L_x']
L_y = dict_save['L_y']
L_z = dict_save['L_z']
sdf_1_map = dict_save['sdf_1_map']
sdf_2_map = dict_save['sdf_2_map']

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Plan simulation

# materials
create_materials() # from dem.py
# create grains and walls
create_grains() # from dem.py
create_walls() # from dem.py
# Engines
create_engines() # from dem.py
# time step
compute_dt() # from dem.py

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# DEM

O.run()
O.wait()

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Output

displacement = np.array(O.bodies[1].state.pos - O.bodies[1].state.refPos)

# Save data
dict_save = {
'displacement': displacement,
}
with open('data/dem_to_main.data', 'wb') as handle:
    pickle.dump(dict_save, handle, protocol=pickle.HIGHEST_PROTOCOL)
