# -*- encoding=utf-8 -*-

from pathlib import Path
import shutil, os, pickle, math
import numpy as np
import matplotlib.pyplot as plt

# own
from pf_to_dem import *

#------------------------------------------------------------------------------------------------------------------------------------------ #

def create_folder(name):
    '''
    Create a new folder. If it already exists, it is erased.
    '''
    if Path(name).exists():
        shutil.rmtree(name)
    os.mkdir(name)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def reduce_n_vtk_files(dict_user, dict_sample):
    '''
    Reduce the number of vtk files for phase-field and dem.

    Warning ! The pf and dem files are not synchronized...
    '''
    if dict_user['n_max_vtk_files'] != None:
        # Phase Field files

        # compute the frequency
        if dict_user['j_total']-1 > dict_user['n_max_vtk_files']:
            f_save = (dict_user['j_total']-1)/(dict_user['n_max_vtk_files']-1)
        else :
            f_save = 1
        # post proccess index
        i_save = 0

        # iterate on time 
        for iteration in range(dict_user['j_total']):
            iteration_str = index_to_str(iteration) # from pf_to_dem.py 
            if iteration >= f_save*i_save:
                i_save_str = index_to_str(i_save) # from pf_to_dem.py
                # rename .pvtu
                os.rename('vtk/pf_'+iteration_str+'.pvtu','vtk/pf_'+i_save_str+'.pvtu')
                # write .pvtu to save all vtk
                file = open('vtk/pf_'+i_save_str+'.pvtu','w')
                file.write('''<?xml version="1.0"?>
                <VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">
                \t<PUnstructuredGrid GhostLevel="1">
                \t\t<PPointData>
                \t\t\t<PDataArray type="Float64" Name="as"/>
                \t\t\t<PDataArray type="Float64" Name="kc"/>
                \t\t\t<PDataArray type="Float64" Name="eta1"/>
                \t\t\t<PDataArray type="Float64" Name="c"/>
                \t\t</PPointData>
                \t\t<PCellData>
                \t\t\t<PDataArray type="Int32" Name="libmesh_elem_id"/>
                \t\t\t<PDataArray type="Int32" Name="subdomain_id"/>
                \t\t\t<PDataArray type="Int32" Name="processor_id"/>
                \t\t</PCellData>
                \t\t<PPoints>
                \t\t\t<PDataArray type="Float64" Name="Points" NumberOfComponents="3"/>
                \t\t</PPoints>''')
                line = ''
                for i_proc in range(dict_user['n_proc']):
                    line = line + '''\t\t<Piece Source="pf_'''+i_save_str+'''_'''+str(i_proc)+'''.vtu"/>\n'''
                file.write(line)
                file.write('''\t</PUnstructuredGrid>
                </VTKFile>''')
                file.close()
                # rename .vtk
                for i_proc in range(dict_user['n_proc']):
                    os.rename('vtk/pf_'+iteration_str+'_'+str(i_proc)+'.vtu','vtk/pf_'+i_save_str+'_'+str(i_proc)+'.vtu')
                i_save = i_save + 1 
            else:
                # delete files
                os.remove('vtk/pf_'+iteration_str+'.pvtu')
                for i_proc in range(dict_user['n_proc']):
                    os.remove('vtk/pf_'+iteration_str+'_'+str(i_proc)+'.vtu')
        # .e file
        os.remove('vtk/pf_out.e')
        # other files
        j = 0
        j_str = index_to_str(j)
        filepath = Path('vtk/pf_other_'+j_str+'.pvtu')
        while filepath.exists():
            for i_proc in range(dict_user['n_proc']):
                os.remove('vtk/pf_other_'+j_str+'_'+str(i_proc)+'.vtu')
            os.remove('vtk/pf_other_'+j_str+'.pvtu')
            j = j + 1
            j_str = index_to_str(j)
            filepath = Path('vtk/pf_other_'+j_str+'.pvtu')

        # DEM files

        # compute the frequency
        if 2*dict_user['n_DEMPF_ite']-1 > dict_user['n_max_vtk_files']:
            f_save = (2*dict_user['n_DEMPF_ite']-1)/(dict_user['n_max_vtk_files']-1)
        else :
            f_save = 1
        # post proccess index
        i_save = 0

        # iterate on time 
        for iteration in range(2*dict_user['n_DEMPF_ite']):
            iteration_str = str(iteration) # from pf_to_dem.py 
            if iteration >= f_save*i_save:
                i_save_str = str(i_save) # from pf_to_dem.py
                os.rename('vtk/2grains_'+iteration_str+'.vtk', 'vtk/2grains_'+i_save_str+'.vtk')
                os.rename('vtk/grain1_'+iteration_str+'.vtk', 'vtk/grain1_'+i_save_str+'.vtk')
                os.rename('vtk/grain2_'+iteration_str+'.vtk', 'vtk/grain2_'+i_save_str+'.vtk')
                i_save = i_save + 1
            else :
                os.remove('vtk/2grains_'+iteration_str+'.vtk')
                os.remove('vtk/grain1_'+iteration_str+'.vtk')
                os.remove('vtk/grain2_'+iteration_str+'.vtk')

#------------------------------------------------------------------------------------------------------------------------------------------ #

def save_mesh_database(dict_user, dict_sample):
    '''
    Save mesh database.
    '''
    # creating a database
    if not Path('mesh_map.database').exists():
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'n_mesh_x': len(dict_sample['x_L']),
        'n_mesh_y': len(dict_sample['y_L']),
        'n_mesh_z': len(dict_sample['z_L']),
        'L_L_i_XYZ_used': dict_sample['L_L_i_XYZ_used'],
        }
        dict_database = {'Run_1': dict_data}
        with open('mesh_map.database', 'wb') as handle:
                pickle.dump(dict_database, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # updating a database
    else :
        with open('mesh_map.database', 'rb') as handle:
            dict_database = pickle.load(handle)
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'n_mesh_x': len(dict_sample['x_L']),
        'n_mesh_y': len(dict_sample['y_L']),
        'n_mesh_z': len(dict_sample['z_L']),
        'L_L_i_XYZ_used': dict_sample['L_L_i_XYZ_used']
        }   
        mesh_map_known = False
        for i_run in range(1,len(dict_database.keys())+1):
            if dict_database['Run_'+str(int(i_run))]['n_proc'] == dict_data['n_proc'] and \
               dict_database['Run_'+str(int(i_run))]['n_mesh_x'] == dict_data['n_mesh_x'] and \
               dict_database['Run_'+str(int(i_run))]['n_mesh_y'] == dict_data['n_mesh_y'] and \
               dict_database['Run_'+str(int(i_run))]['n_mesh_z'] == dict_data['n_mesh_z']:
                mesh_map_known = True
        # new entry
        if not mesh_map_known: 
            key_entry = 'Run_'+str(int(len(dict_database.keys())+1))
            dict_database[key_entry] = dict_data
            with open('mesh_map.database', 'wb') as handle:
                pickle.dump(dict_database, handle, protocol=pickle.HIGHEST_PROTOCOL)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def check_mesh_database(dict_user, dict_sample):
    '''
    Check mesh database.
    '''
    if Path('mesh_map.database').exists():
        with open('mesh_map.database', 'rb') as handle:
            dict_database = pickle.load(handle)
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'n_mesh_x': len(dict_sample['x_L']),
        'n_mesh_y': len(dict_sample['y_L']),
        'n_mesh_z': len(dict_sample['z_L'])
        }   
        mesh_map_known = False
        for i_run in range(1,len(dict_database.keys())+1):
            if dict_database['Run_'+str(int(i_run))]['n_proc'] == dict_data['n_proc'] and \
               dict_database['Run_'+str(int(i_run))]['n_mesh_x'] == dict_data['n_mesh_x'] and \
               dict_database['Run_'+str(int(i_run))]['n_mesh_y'] == dict_data['n_mesh_y'] and \
               dict_database['Run_'+str(int(i_run))]['n_mesh_z'] == dict_data['n_mesh_z']:
                mesh_map_known = True
                i_known = i_run
        if mesh_map_known :
            dict_sample['Map_known'] = True
            dict_sample['L_L_i_XYZ_used'] = dict_database['Run_'+str(int(i_known))]['L_L_i_XYZ_used']
        else :
            dict_sample['Map_known'] = False
    else :
        dict_sample['Map_known'] = False

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_sum_mean_etai_c(dict_user, dict_sample):
    '''
    Plot figure illustrating the sum and the mean of etai and c.
    '''
    # compute tracker
    dict_user['L_sum_eta_1'].append(np.sum(dict_sample['eta_1_map']))
    dict_user['L_sum_eta_2'].append(np.sum(dict_sample['eta_2_map']))
    dict_user['L_sum_c'].append(np.sum(dict_sample['c_map']))
    dict_user['L_sum_mass'].append(1/dict_user['V_m']*np.sum(dict_sample['eta_1_map'])+1/dict_user['V_m']*np.sum(dict_sample['eta_2_map'])+np.sum(dict_sample['c_map']))
    dict_user['L_m_eta_1'].append(np.mean(dict_sample['eta_1_map']))
    dict_user['L_m_eta_2'].append(np.mean(dict_sample['eta_2_map']))
    dict_user['L_m_c'].append(np.mean(dict_sample['c_map']))
    dict_user['L_m_mass'].append(1/dict_user['V_m']*np.mean(dict_sample['eta_1_map'])+1/dict_user['V_m']*np.mean(dict_sample['eta_2_map'])+np.mean(dict_sample['c_map']))

    # plot sum eta_i, c
    if 'sum_etai_c' in dict_user['L_figures']:
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))
        ax1.plot(dict_user['L_sum_eta_1'])
        ax1.set_title(r'$\Sigma\eta_1$')
        ax2.plot(dict_user['L_sum_eta_2'])
        ax2.set_title(r'$\Sigma\eta_2$')
        ax3.plot(dict_user['L_sum_c'])
        ax3.set_title(r'$\Sigma C$')
        ax4.plot(dict_user['L_sum_mass'])
        ax4.set_title(r'$\Sigma mass$')
        fig.tight_layout()
        fig.savefig('plot/sum_etai_c.png')
        plt.close(fig)

    # plot mean eta_i, c
    if 'mean_etai_c' in dict_user['L_figures']:
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))
        ax1.plot(dict_user['L_m_eta_1'])
        ax1.set_title(r'Mean $\eta_1$')
        ax2.plot(dict_user['L_m_eta_2'])
        ax2.set_title(r'Mean $\eta_2$')
        ax3.plot(dict_user['L_m_c'])
        ax3.set_title(r'Mean $c$')
        ax4.plot(dict_user['L_m_mass'])
        ax4.set_title(r'Mean mass')
        fig.tight_layout()
        fig.savefig('plot/mean_etai_c.png')
        plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def compute_mass(dict_user, dict_sample):
    '''
    Compute the mass at a certain time.
     
    Mass is sum of etai and c.
    '''
    # sum of masses
    dict_user['sum_eta_1_tempo'] = np.sum(dict_sample['eta_1_map'])
    dict_user['sum_eta_2_tempo'] = np.sum(dict_sample['eta_2_map'])
    dict_user['sum_c_tempo'] = np.sum(dict_sample['c_map'])
    dict_user['sum_mass_tempo'] = 1/dict_user['V_m']*np.sum(dict_sample['eta_1_map'])+1/dict_user['V_m']*np.sum(dict_sample['eta_2_map'])+np.sum(dict_sample['c_map'])
    
#------------------------------------------------------------------------------------------------------------------------------------------ #

def compute_mass_loss(dict_user, dict_sample, tracker_key):
    '''
    Compute the mass loss from the previous compute_mass() call.
     
    Plot in the given tracker.
    Mass is sum of etai and c.
    '''
    # delta masses
    deta1 = np.sum(dict_sample['eta_1_map']) - dict_user['sum_eta_1_tempo']
    deta2 = np.sum(dict_sample['eta_2_map']) - dict_user['sum_eta_2_tempo']
    dc = np.sum(dict_sample['c_map']) - dict_user['sum_c_tempo']
    dm = 1/dict_user['V_m']*np.sum(dict_sample['eta_1_map'])+1/dict_user['V_m']*np.sum(dict_sample['eta_2_map'])+np.sum(dict_sample['c_map']) - dict_user['sum_mass_tempo']
    
    # save
    dict_user[tracker_key+'_eta1'].append(deta1)
    dict_user[tracker_key+'_eta2'].append(deta2)
    dict_user[tracker_key+'_c'].append(dc)
    dict_user[tracker_key+'_m'].append(dm)

    # plot
    if 'mass_loss' in dict_user['L_figures']:
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))
        ax1.plot(dict_user[tracker_key+'_eta1'])
        ax1.set_title(r'$\eta_1$ loss')
        ax2.plot(dict_user[tracker_key+'_eta2'])
        ax2.set_title(r'$\eta_2$ loss')
        ax3.plot(dict_user[tracker_key+'_c'])
        ax3.set_title(r'$c$ loss')
        ax4.plot(dict_user[tracker_key+'_m'])
        ax4.set_title(r'$\eta_1$ + $c$ loss')
        fig.tight_layout()
        fig.savefig('plot/'+tracker_key+'.png')
        plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_performances(dict_user, dict_sample):
    '''
    Plot figure illustrating the time performances of the algorithm.
    '''
    if 'performances' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.plot(dict_user['L_t_dem'], label='DEM')
        ax1.plot(dict_user['L_t_pf'], label='PF')
        ax1.plot(dict_user['L_t_dem_to_pf'], label='DEM to PF')
        ax1.plot(dict_user['L_t_pf_to_dem_1'], label='PF to DEM 1')
        ax1.plot(dict_user['L_t_pf_to_dem_2'], label='PF to DEM 2')
        ax1.legend()
        ax1.set_title('Performances (s)')
        ax1.set_xlabel('Iterations (-)')
        fig.tight_layout()
        fig.savefig('plot/performances.png')
        plt.close(fig)
                
#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_dem(dict_user, dict_sample):
    '''
    Plot figure illustrating the overlap and force transmitted.
    '''
    if 'overlap' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.plot(dict_user['L_overlap'])
        ax1.set_xlabel('iterations (-)')
        ax1.set_ylabel('overlap in DEM (-)')
        fig.tight_layout()
        fig.savefig('plot/dem_overlap.png')
        plt.close(fig)
    if 'normal_force' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.plot(dict_user['L_normal_force'], label='DEM')
        ax1.plot([0, len(dict_user['L_normal_force'])-1],\
                 [dict_user['force_applied'], dict_user['force_applied']], label='target')
        ax1.legend()
        ax1.set_xlabel('iterations (-)')
        ax1.set_ylabel('normal force (-)')
        fig.tight_layout()
        fig.savefig('plot/dem_normal_force.png')
        plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_contact(dict_user, dict_sample):
    '''
    Plot figure illustrating the contact characteristics.
    '''
    if 'contact_box' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.plot(dict_user['L_contact_box_x'], label='x')
        ax1.plot(dict_user['L_contact_box_y'], label='y')
        ax1.plot(dict_user['L_contact_box_z'], label='z')
        ax1.legend()
        ax1.set_xlabel('iterations (-)')
        ax1.set_ylabel('contact box dimensions (-)')
        fig.tight_layout()
        fig.savefig('plot/contact_box.png')
        plt.close(fig)
    if 'contact_volume' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.plot(dict_user['L_contact_volume'])
        ax1.set_xlabel('iterations (-)')
        ax1.set_ylabel('contact volume (-)')
        fig.tight_layout()
        fig.savefig('plot/contact_volume.png')
        plt.close(fig)
    if 'contact_surface' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.plot(dict_user['L_contact_surface'], label='DEM')
        ax1.set_xlabel('iterations (-)')
        ax1.set_ylabel('contact surface (-)')
        fig.tight_layout()
        fig.savefig('plot/contact_surface.png')
        plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_as_pressure(dict_user, dict_sample):
    '''
    Plot figure illustrating the solid activity and pressure at the contact.
    '''
    if 'as' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.plot(dict_user['L_contact_as'])
        ax1.set_xlabel('iterations (-)')
        ax1.set_ylabel('solid activity (-)')
        fig.tight_layout()
        fig.savefig('plot/contact_as.png')
        plt.close(fig)
    if 'pressure' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.plot(dict_user['L_contact_pressure'])
        ax1.set_xlabel('iterations (-)')
        ax1.set_ylabel('solid pressure (-)')
        fig.tight_layout()
        fig.savefig('plot/contact_pressure.png')
        plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_slices(dict_user, dict_sample):
    '''
    Plot slices of the configuration.
    '''
    # determine the index of the slices
    i_x_slice = int(len(dict_sample['x_L'])/2)
    i_y_slice = int(len(dict_sample['y_L'])/2)
    i_z_slice = int(len(dict_sample['z_L'])/2)

    if 'configuration_eta' in dict_user['L_figures'] :
        # plot
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
        # y-z
        im = ax1.imshow(dict_sample['eta_1_map'][i_x_slice,:,:], interpolation = 'nearest')
        fig.colorbar(im, ax=ax1)
        ax1.set_title(r'slice y-z',fontsize = 30)
        # x-z
        im = ax2.imshow(dict_sample['eta_1_map'][:,i_y_slice,:], interpolation = 'nearest')
        fig.colorbar(im, ax=ax2)
        ax2.set_title(r'slice x-z',fontsize = 30)
        # x-y
        im = ax3.imshow(dict_sample['eta_1_map'][:,:,i_z_slice], interpolation = 'nearest')
        fig.colorbar(im, ax=ax3)
        ax3.set_title(r'slice x-y',fontsize = 30)
        # close
        fig.tight_layout()
        fig.savefig('plot/configuration/eta1_'+str(dict_sample['i_DEMPF_ite'])+'.png')
        plt.close(fig)
        
        # plot
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
        # y-z
        im = ax1.imshow(dict_sample['eta_2_map'][i_x_slice,:,:], interpolation = 'nearest')
        fig.colorbar(im, ax=ax1)
        ax1.set_title(r'slice y-z',fontsize = 30)
        # x-z
        im = ax2.imshow(dict_sample['eta_2_map'][:,i_y_slice,:], interpolation = 'nearest')
        fig.colorbar(im, ax=ax2)
        ax2.set_title(r'slice x-z',fontsize = 30)
        # x-y
        im = ax3.imshow(dict_sample['eta_2_map'][:,:,i_z_slice], interpolation = 'nearest')
        fig.colorbar(im, ax=ax3)
        ax3.set_title(r'slice x-y',fontsize = 30)
        # close
        fig.tight_layout()
        fig.savefig('plot/configuration/eta2_'+str(dict_sample['i_DEMPF_ite'])+'.png')
        plt.close(fig)

    if 'configuration_c' in dict_user['L_figures']:
        # plot
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
        # y-z
        im = ax1.imshow(dict_sample['c_map'][i_x_slice,:,:], interpolation = 'nearest')
        fig.colorbar(im, ax=ax1)
        ax1.set_title(r'slice y-z',fontsize = 30)
        # x-z
        im = ax2.imshow(dict_sample['c_map'][:,i_y_slice,:], interpolation = 'nearest')
        fig.colorbar(im, ax=ax2)
        ax2.set_title(r'slice x-z',fontsize = 30)
        # x-y
        im = ax3.imshow(dict_sample['c_map'][:,:,i_z_slice], interpolation = 'nearest')
        fig.colorbar(im, ax=ax3)
        ax3.set_title(r'slice x-y',fontsize = 30)
        # close
        fig.tight_layout()
        fig.savefig('plot/configuration/c_'+str(dict_sample['i_DEMPF_ite'])+'.png')
        plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_displacement(dict_user, dict_sample):
    '''
    Plot figure illustrating the cumulative displacement.
    '''
    # pp data
    L_strain = []
    for i_displacement in range(len(dict_user['L_displacement'])):
        if i_displacement == 0:
            L_displacement_cum = [dict_user['L_displacement'][i_displacement][2]]
        else : 
            L_displacement_cum.append(L_displacement_cum[-1]+dict_user['L_displacement'][i_displacement][2])
        L_strain.append(L_displacement_cum[-1]/(2*dict_user['radius']))
    # plot
    fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
    ax1.plot(L_strain)
    ax1.set_xlabel('iterations (-)')
    ax1.set_ylabel('vertical strain (-)')
    fig.tight_layout()
    fig.savefig('plot/vertical_strain.png')
    plt.close(fig)
