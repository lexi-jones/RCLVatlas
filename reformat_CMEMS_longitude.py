# Converts the longitude from CMEMS data to positive longitudes (0 to 360) instead of -180 to 180 
#
# Lexi Jones

import os
import xarray as xr
import numpy as np
from config import *

for data_file in os.listdir(gos_vel_dir):
    ds = xr.open_dataset(gos_vel_dir + data_file)
    ds['longitude'] = ds['longitude'] + 180
    new_file_name = 'dt_global_allsat_phy_l4_' + str(np.asarray(ds['time'])[0])[0:10].replace('-','') + '.nc'
    ds.to_netcdf(gos_vel_dir + new_file_name)
