# Converts the longitude from CMEMS data to positive longitudes (0 to 360) instead of -180 to 180.
# Need to then wrap all of the variables to be accordingly centered on the Pacific Ocean.
#
# Lexi Jones

import os
import xarray as xr
import numpy as np
from config import *

for data_file in os.listdir(gos_vel_dir):
    ds = xr.open_dataset(gos_vel_dir + data_file)

    new_lon = []
    for l in np.array(ds.longitude): 
        if l < 0:
            new_lon.append(l+360)
        else:
            new_lon.append(l)
    lon = np.array(new_lon[720:] + new_lon[0:720])
    ugos = np.concatenate((np.array(ds.ugos)[0,:,720:],np.array(ds.ugos)[0,:,:720]),axis=1)[np.newaxis,:,:]
    vgos = np.concatenate((np.array(ds.vgos)[0,:,720:],np.array(ds.vgos)[0,:,:720]),axis=1)[np.newaxis,:,:]
    ugosa = np.concatenate((np.array(ds.ugosa)[0,:,720:],np.array(ds.ugosa)[0,:,:720]),axis=1)[np.newaxis,:,:]
    vgosa = np.concatenate((np.array(ds.vgosa)[0,:,720:],np.array(ds.vgosa)[0,:,:720]),axis=1)[np.newaxis,:,:]
    sla = np.concatenate((np.array(ds.sla)[0,:,720:],np.array(ds.sla)[0,:,:720]),axis=1)[np.newaxis,:,:]

    # Edit the longitude attributes
    ds.attrs['geospatial_lon_min'] = min(lon)
    ds.attrs['geospatial_lon_max'] = max(lon)
    ds['longitude'].attrs['valid_min'] = min(lon)
    ds['longitude'].attrs['valid_max'] = max(lon)

    # Write new netCDF file
    new_ds = xr.Dataset({'ugos':(['t','y','x'], ugos),
                     'vgos':(['t','y','x'], vgos),
                    'ugosa':(['t','y','x'], ugosa),
                        'vgosa':(['t','y','x'], vgosa),
                        'sla':(['t','y','x'], sla)},
			 coords = {'time':(['t'], ds.time),
				   'latitude':(['y'], ds.latitude),
				   'longitude':(['x'], lon)})

    new_ds['ugos'].attrs = ds['ugos'].attrs
    new_ds['vgos'].attrs = ds['vgos'].attrs
    new_ds['ugosa'].attrs = ds['ugosa'].attrs
    new_ds['vgosa'].attrs = ds['vgosa'].attrs
    new_ds['sla'].attrs = ds['sla'].attrs
    new_ds['time'].attrs = ds['time'].attrs
    new_ds['latitude'].attrs = ds['latitude'].attrs
    new_ds['longitude'].attrs = ds['longitude'].attrs
    new_ds.attrs = ds.attrs

    new_file_name = 'dt_global_allsat_phy_l4_' + str(np.asarray(ds['time'])[0])[0:10].replace('-','') + '.nc'
    new_ds.to_netcdf(gos_vel_dir + new_file_name)
