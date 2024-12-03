# -*- encoding=utf-8 -*-

from yade import pack, utils, plot, export
import pickle
import numpy as np

# -----------------------------------------------------------------------------#
# Functions called once
# -----------------------------------------------------------------------------#

def create_materials():
    '''
    Create materials.
    '''
    O.materials.append(FrictMat(young=E, poisson=Poisson, frictionAngle=atan(0.5), density=density, label='frictMat'))
    O.materials.append(FrictMat(young=E, poisson=Poisson, frictionAngle=0, density=density, label='frictlessMat'))

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
    for i_grain in range(len(L_sdf_i_map)):
        O.bodies.append(
            levelSetBody(grid=grid,
                        distField=L_sdf_i_map[i_grain].tolist(),
                        material=0)
        )
        if i_grain == L_id_fixed: 
            O.bodies[-1].state.blockedDOFs = 'xyzXYZ'
        else :
            O.bodies[-1].state.blockedDOFs = 'XYZ'

# -----------------------------------------------------------------------------#

def create_walls():
    '''
    Recreate walls.
    '''
    for i_wall in range(len(L_pos_w)):
        O.bodies.append(wall((L_pos_w[i_wall][0], L_pos_w[i_wall][1], L_pos_w[i_wall][2]), L_pos_w[i_wall][3], material=1))

# -----------------------------------------------------------------------------#

def create_plots():
    '''
    Create plots during the DEM step.
    '''
    plot.plots = {'iteration': ('unbalForce')}

# -----------------------------------------------------------------------------#

def compute_dt():
    '''
    Compute the time step used in the DEM step.
    '''
    O.dt = 0.2*SpherePWaveTimeStep(radius=radius, density=density, young=E)

# -----------------------------------------------------------------------------#

def create_engines():
    '''
    Create engines.

    Overlap based on the distance

    Ip2:
        kn = given
        ks = given    

    Law2:
        Fn = kn.un
        Fs = ks.us
    '''
    O.engines = [
            #VTKRecorder(recorders=["lsBodies"], fileName='./vtk/ite_PFDEM_'+str(i_DEMPF_ite)+'_', iterPeriod=0, multiblockLS=True, label='initial_export'),
            ForceResetter(),
            PyRunner(command='applied_force()',iterPeriod=1, label='apply_force'),
            InsertionSortCollider([Bo1_LevelSet_Aabb(), Bo1_Wall_Aabb()], verletDist=0.00),
            InteractionLoop(
                    [Ig2_LevelSet_LevelSet_ScGeom(), Ig2_Wall_LevelSet_ScGeom()],
                    [Ip2_FrictMat_FrictMat_FrictPhys(kn=MatchMaker(algo='val', val=kn), ks=MatchMaker(algo='val', val=ks))],
                    [Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=False)]),
            NewtonIntegrator(damping=0.1, label='newton', gravity=(0, 0, 0)),
    		PyRunner(command='add_data()',iterPeriod=1, label='data'),
            PyRunner(command='check()',iterPeriod=1, label='checker')
    ]

    if print_vtk:
        O.engines = O.engines + [VTKRecorder(recorders=["lsBodies"], fileName='./vtk/yade/', iterPeriod=100, multiblockLS=True, label='export')]

# -----------------------------------------------------------------------------#
# Functions called multiple times
# -----------------------------------------------------------------------------#

def applied_force():
    '''
    Apply a constant force on the top grain.
    '''
    for id_applied in L_id_applied:
        O.forces.addF(id_applied, (0, 0, -force_applied))

# -----------------------------------------------------------------------------#

def check():
    '''
    Try to detect a steady-state.
    A maximum number of iteration is used.
    '''
    if O.iter < max(n_ite_max*0.01, n_steady_state_detection):
        return
    window = plot.data['unbalForce'][-n_steady_state_detection:]
    if O.iter > n_ite_max or max(window)<steady_state_detection:
        if print_all_dem:
            plot.plot(noShow=True).savefig('plot/dem/'+str(i_DEMPF_ite)+'.png')
        if print_dem:
            plot.plot(noShow=True).savefig('plot/dem.png')
        O.pause() # stop DEM simulation
    
# -----------------------------------------------------------------------------#

def add_data():
    '''
    Add data to plot :
        - iteration
        - unbalannced force (mean resultant forces / mean contact force)
    '''
    plot.addData(iteration=O.iter, unbalForce=unbalancedForce())
    
# -----------------------------------------------------------------------------#
# Load data
# -----------------------------------------------------------------------------#

# from main
with open('data/main_to_dem.data', 'rb') as handle:
    dict_save = pickle.load(handle)
radius = dict_save['radius']
E = dict_save['E']
Poisson = dict_save['Poisson']
kn = dict_save['kn']
ks = dict_save['ks']
n_ite_max = dict_save['n_ite_max']
i_DEMPF_ite = dict_save['i_DEMPF_ite']
L_pos_w = dict_save['L_pos_w']
force_applied = dict_save['force_applied']
n_steady_state_detection = dict_save['n_steady_state_detection']
steady_state_detection = dict_save['steady_state_detection']
print_all_dem = dict_save['print_all_dem']
print_dem = dict_save['print_dem']
print_vtk = dict_save['print_vtk']
L_id_fixed = dict_save['L_id_fixed']
L_id_applied = dict_save['L_id_applied']
density = 2000

# from plane interpolation
with open('data/level_sets.data', 'rb') as handle:
    dict_save = pickle.load(handle)
L_x = dict_save['L_x']
L_y = dict_save['L_y']
L_z = dict_save['L_z']
L_sdf_i_map = dict_save['L_sdf_i_map']

# -----------------------------------------------------------------------------#
# Plan simulation
# -----------------------------------------------------------------------------#

# materials
create_materials()
# create grains and walls
create_grains()
create_walls() 
# Engines
create_engines()
# time step
compute_dt()
# plot
create_plots()

# -----------------------------------------------------------------------------#
# MAIN DEM
# -----------------------------------------------------------------------------#

O.run()
O.wait()

# -----------------------------------------------------------------------------#
# Output
# -----------------------------------------------------------------------------#

L_displacement = []
for i_grain in range(len(L_sdf_i_map)):
    L_displacement.append(np.array(O.bodies[i_grain].state.pos - O.bodies[i_grain].state.refPos))
L_contact = []
for i in O.interactions:
    # work only on grain-grain contact
    if isinstance(O.bodies[i.id1].shape, LevelSet) and isinstance(O.bodies[i.id2].shape, LevelSet):
        if i.id1 < i.id2:
            id_smaller = i.id1
            id_larger  = i.id2
        else:
            id_smaller = i.id2
            id_larger  = i.id1
        L_contact.append([id_smaller, id_larger, i.geom.penetrationDepth, np.array([i.phys.normalForce[0], i.phys.normalForce[1], i.phys.normalForce[2]])])

# Save data
dict_save = {
'L_displacement': L_displacement,
'L_contact': L_contact
}
with open('data/dem_to_main.data', 'wb') as handle:
    pickle.dump(dict_save, handle, protocol=pickle.HIGHEST_PROTOCOL)
