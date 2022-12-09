# Lexi Jones
# Date Created: 07/15/21
# Last Edited: 12/09/22

# Run an OceanParcels simulation with functions from 'functions_for_parcels'. 
# Format: python run_parcels_CMEMS_v3.py lat_start lat_stop lon_start lon_stop spatial_step date_input time_step output_freq runtime runtime_unit backwards
#	lat_start,lat_stop,lon_start,lon_stop: latitude and longitude bounds to intialize particles
#	spatial_step: lat/lon spacing between initial particles
#	date_input: start date with format YYYYMMDD
#	time_step: (minutes) the amount of time that passes in the fieldset before calculating new particle positions & attributes
#	output_freq: (hours) frequency to output information about the particle positions & attributes
#	runtime: how long to run the simulation
#	runtime_units: 'days' or 'hours', corresponds with runtime
#	backwards: 'y' if the particles should be run backwards in time, 'n' if the particles should be run forwards in time

import numpy as np
import time,sys
from glob import glob
from datetime import datetime
from parcels import FieldSet,ParticleSet,Variable,JITParticle,AdvectionRK4,ErrorCode
from functions_for_parcels import DeleteParticle,particle_grid2d,simulate_particles2d

#Collect User Input
#Note: sys.argv[0] is the name of the script
print(sys.argv)

lat_start,lat_stop = float(sys.argv[1]), float(sys.argv[2])
lon_start,lon_stop = float(sys.argv[3]), float(sys.argv[4])
spatial_step = float(sys.argv[5])

date_input = sys.argv[6]
start_year,start_month,start_day  = int(date_input[0:4]),int(date_input[4:6]),int(date_input[6:8])
start_date = datetime(start_year,start_month,start_day) # format datetime

time_step,output_freq = int(sys.argv[7]),int(sys.argv[8])
runtime,runtime_unit = int(sys.argv[9]),sys.argv[10]
backwards = str(sys.argv[11])

########## PARCELS #############

### Create Parcels fieldset ###
input_dir = '/nfs/micklab005/jonesae/gradients4/CMEMS_data/' 
parcels_input_files = sorted(glob(input_dir+'nrt_global_allsat_phy_l4_cropped_*.nc'))
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
    
pset_dynamic,num_particles = particle_grid2d(fieldset,SpinnyParticle,[lat_start,lat_stop,spatial_step],[lon_start,lon_stop,spatial_step],start_date)

output_dir = '/nfs/micklab005/jonesae/gradients4/parcels_trajs/'
output_file_path = '%s%s_%s%s_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc'%(
output_dir,str(start_date)[0:10],runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
print('Output file: %s'%(output_file_path))

### Execute particle simulation ###
start_time = time.time()
simulate_particles2d(pset_dynamic,output_file_path,runtime,runtime_unit,time_step,output_freq,backwards)
print("--- Simulation runtime: %s minutes ---" % (round((time.time() - start_time)/60),5))
