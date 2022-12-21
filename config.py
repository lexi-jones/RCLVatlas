# Configuration file
# The user should define directory paths and adjust parameters here.

import numpy as np

### Directory paths ### 
gos_vel_dir = '/nfs/micklab005/jonesae/NP_2000_to_2020/CMEMS_data/' # Geostrophic velocity directory
lag_traj_dir = '/nfs/micklab005/jonesae/NP_2000_to_2020/RCLVatlas_test/lag_trajs/' # Lagrangian trajectory directory
LAVD_dir = '/nfs/micklab005/jonesae/NP_2000_to_2020/RCLVatlas_test/LAVD/' # LAVD directory 
RCLV_dir = '/nfs/micklab005/jonesae/NP_2000_to_2020/RCLVatlas_test/RCLVs/' #RCLVatlas directory

grid_bounds = {'lon_bound_west':199.0,
               'lon_bound_east':204.0,
               'lat_bound_south':17.0,
               'lat_bound_north':22.0,
               'lag_grid_res':0.03125}

sim_params = {'runtime':8, 
              'runtime_unit':'days',
              'timestep':20, 
              'output_freq':6,
              'backwards':'y'} 

# NOTE: min_dist and min_area are in units of pixels, relative to the LAVD field
RCLV_params = {'min_dist':24,
               'min_area':104,
               'init_contour_step_frac':0.1,
               'convex_def_tol':0.001}

filename_str = '%sdays_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq'%(
    sim_params['runtime'],sim_params['timestep'],grid_bounds['lat_bound_south'],grid_bounds['lat_bound_north'],
    grid_bounds['lon_bound_west'],grid_bounds['lon_bound_east'],grid_bounds['lag_grid_res'],sim_params['output_freq'])

traj_lon_array = np.arange(grid_bounds['lon_bound_west'],grid_bounds['lon_bound_east'],grid_bounds['lag_grid_res'])
traj_lat_array = np.arange(grid_bounds['lat_bound_south'],grid_bounds['lat_bound_north'],grid_bounds['lag_grid_res'])
