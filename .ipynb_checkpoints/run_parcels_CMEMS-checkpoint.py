# Lexi Jones
# Date Created: 07/15/21
# Last Edited: 12/20/22

# Run an OceanParcels simulation

import time,sys
import numpy as np
import xarray as xr
from glob import glob
from datetime import datetime
from parcels import FieldSet,Variable,JITParticle
from functions_for_parcels import *
from config import *

date_input = sys.argv[1] # user input: particle intialization date

start_year,start_month,start_day = int(str(date_input)[0:4]),int(str(date_input)[4:6]),int(str(date_input)[6:8])
start_date = datetime(start_year,start_month,start_day) # format datetime

### Create Parcels fieldset ###
parcels_input_files = sorted(glob(gos_vel_dir+'dt_global_allsat_phy_l4_*.nc'))
filenames = {'U': parcels_input_files,'V': parcels_input_files}
variables = {'U': 'ugos','V': 'vgos'} #name of the velocity variables in the netCDF file
dimensions = {'U': {'lon':'longitude','lat':'latitude','time':'time'},
              'V': {'lon':'longitude','lat':'latitude','time':'time'}}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
print('Fieldset created.')

### Create particleset ###
class SpinnyParticle(JITParticle):
    u = Variable('u',dtype=np.float64)
    v = Variable('v',dtype=np.float64)
    vort = Variable('vort',dtype=np.float64)
    
pset_dynamic,num_particles = particle_grid2d(fieldset,SpinnyParticle,
                                             [grid_bounds['lat_bound_south'],grid_bounds['lat_bound_north'],grid_bounds['lag_grid_res']],
                                             [grid_bounds['lon_bound_west'],grid_bounds['lon_bound_east'],grid_bounds['lag_grid_res']],
                                             start_date)

### Execute particle simulation ###
print("Running Lagrangian simulation ...")
traj_output_file_path = lag_traj_dir + str(date_input) + '_' + filename_str + '.nc'
simulate_particles2d(pset_dynamic,traj_output_file_path,
                     sim_params['runtime'],sim_params['runtime_unit'],
                     sim_params['timestep'],sim_params['output_freq'],
                     sim_params['backwards'])
print('Trajectory output file: %s'%(traj_output_file_path))

### Calculate Lagrangian average vorticity deviation (LAVD) ###
print("Calculating the LAVD ...")
traj_ds = xr.open_dataset(traj_output_file_path) # open the Lagrangian trajectory dataset that was just produced
vort_premask = traj_ds.variables["vort"]
vort = np.array(vort_premask.where(vort_premask != 0)) #filters out land values
LAVD = calc_LAVD(vort,sim_params['output_freq'],sim_params['runtime'])
LAVD_output_file_path = LAVD_dir + str(date_input) + '_LAVD_' + filename_str + '.npy'
np.save(LAVD_output_file_path,LAVD)
print('LAVD output file: %s'%(LAVD_output_file_path))
