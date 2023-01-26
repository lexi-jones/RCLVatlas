# RCLVatlas Pipeline
# 1. Identify RCLVS using floater package
#      - Checks multiple CDs (0.01,0.02,0.03)
#      - Sets dispersal threshold with CI (>=-0.5)
#      - 85% of particles must have the same sign of the vorticity
# 2. Track the RCLVs through time and ID them
# 3. QC: Check if any RCLVs skipped a date (or 2)
#      - Interpolate over the missing timestep using particle trajectories
# 4. Interpolate eddy bounds for first 3 timesteps
#      - Tracks the RCLV birth 
# 5. QC: Address instances of overlapping contours
# 6. Give RCLVs an age

# Lexi Jones
# Last Edited: 12/21/22

import os,itertools,sys
import xarray as xr
import numpy as np
from floater import rclv
from skimage.feature import peak_local_max
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon 
from subfunctions_for_RCLV_atlas import *

sys.path.append('../')
from config import *

####################################### 1. Identify RCLVS #######################################

def LAVD_to_RCLV(initial_date):
    """
    Extract the RCLVs from the LAVD using the floater package (https://github.com/ocean-transport/floater). 
    To qualify, a feature must have a convex deficiency <= 0.03, area >= min_area, coherency index >= -0.5, and
    85% of the particles need to have the same sign of the vorticity on day 1 and the last day of the sim. 
    
    Input
        initial_date: Lagrangian particle initialization date; determines which LAVD file to open
    Output
        plms: Indices of the peak local max of the LAVD, i.e. center of an RCLV
        cons: Boundaries (in grid index format) for the RCLV contour
        areas: Areas (km^2) of the RCLVs
        cds: Convex deficiencies of the RCLVs
        pols: Polarity (cyclonic or anticyclonic) of the vortex
    
    """
    # Load data
    traj = xr.open_dataset('%s%s_%s.nc'%(lag_traj_dir,initial_date,filename_str)) #Load Lagrangian trajectories
    LAVD = np.load('%s%s_LAVD_%s.npy'%(LAVD_dir,initial_date,filename_str)) #Load LAVD data
    LAVD = np.ma.masked_where(np.isnan(LAVD),LAVD) #Land mask required for the peak_local_max function to work
    LAVD_reshape = np.transpose(np.reshape(LAVD,(len(traj_lon_array),len(traj_lat_array))))
    
    # Iterate through all of the local maxima to find structures that meet our RCLV criteria 
    plms,cons,areas,cds,pols = [],[],[],[],[] # arrays to store the output
    plm = peak_local_max(LAVD_reshape,min_distance=RCLV_params['min_dist']) # grid indices of the local maxima in the LAVD field
    for ji in plm:
        args = {'lon0':traj_lon_array[ji[1]],
                'lat0':traj_lat_array[ji[0]],
                'dlon':np.abs(traj_lon_array[1]-traj_lon_array[0]),
                'dlat':np.abs(traj_lat_array[1]-traj_lat_array[0])}
        
        for target_CD in [0.03,0.02,0.01]: #target CDs to test, will iterate through each of these unless if a successful contour is identified
            con,area,cd = rclv.convex_contour_around_maximum(LAVD_reshape,ji,
                                                             init_contour_step_frac=RCLV_params['init_contour_step_frac'],
                                                             convex_def_tol=RCLV_params['convex_def_tol'],
                                                             convex_def=target_CD)
            
            try: # Code errors if the local maximum is too close to the land where there is missing values
                if (area >= RCLV_params['min_area']) and (cd <= 0.03): # check area and convex deficiency meet thresholds; actual cd identified will be different from target CD
                    lon_inds = [round(i) for i in con[:,1]]
                    lat_inds = [round(i) for i in con[:,0]]

                    # Reformat points of the polygon with parentheses
                    poly_pts = [(traj_lon_array[lon_inds][pt],traj_lat_array[lat_inds][pt]) for pt in np.arange(0,len(lon_inds))]

                    # Get the CI of the particles after simulation run time to make sure it meets threshold set
                    x_mask,y_mask = find_polygon_pts(poly_pts,traj_lon_array,traj_lat_array)                    
                    eddy_day0_lons,eddy_day0_lats,_ = extract_particles_after_time(traj,x_mask,y_mask,traj_lat_array,sim_params,0) #day 0 data
                    eddy_dayx_lons,eddy_dayx_lats,eddy_dayx_vorts = extract_particles_after_time(traj,x_mask,y_mask,traj_lat_array,sim_params,sim_params['runtime']) # last day
                    contour_CI = CI(eddy_day0_lons,eddy_day0_lats,eddy_dayx_lons,eddy_dayx_lats) #calculate the coherency index (measure of dispersal)

                    if (contour_CI >= -0.5): # Coherency index criteria, minimal dispersal required
                        # Check that the vorticty is consistent within the vortex over the course of the simulation (filters out some weird saddle point features,etc)
                        # Timestep to get initial vorticity needs to be day 1 since at day 0 all particles are initialized with vort=0
                        eddy_day1_lons,eddy_day1_lats,eddy_day1_vorts = extract_particles_after_time(traj,x_mask,y_mask,traj_lat_array,sim_params,1) 
                        day1_vort_signs = [np.sign(v) for v in eddy_day1_vorts]
                        day1_anti_count,day1_cyc_count = day1_vort_signs.count(-1.0),day1_vort_signs.count(1.0)  
                        day1_vort_freq = np.max((day1_anti_count,day1_cyc_count))/(day1_anti_count+day1_cyc_count) #percentage of particles with the same vorticity

                        if (day1_vort_freq >= 0.85): # 85% or more have the same vort sign on day 1
                            dayx_vort_signs = [np.sign(v) for v in eddy_dayx_vorts]
                            dayx_anti_count,dayx_cyc_count = dayx_vort_signs.count(-1.0),dayx_vort_signs.count(1.0)  
                            dayx_vort_freq = np.max((dayx_anti_count,dayx_cyc_count))/(dayx_anti_count+dayx_cyc_count)

                            if day1_anti_count > day1_cyc_count: #get the polarity of the vortex on day 1
                                day1_vort_pol = 'anti'
                            else:
                                day1_vort_pol = 'cyc'
                            if dayx_anti_count > dayx_cyc_count: #get the polarity of the vortex on last day
                                dayx_vort_pol = 'anti'
                            else:
                                dayx_vort_pol = 'cyc'

                            # If all criteria are passed, add values to output arrays & move on to next local maxima
                            if (dayx_vort_freq >= 0.85) and (day1_vort_pol == dayx_vort_pol): # Final criteria: vort is consistent for 85% of particles on the last day
                                region_area = rclv.polygon_area(rclv.project_vertices(con, **args)) #Convert area units from pixels to m^2
                                plms.append(ji)
                                cons.append(con)
                                areas.append(region_area/(10**6)) #convert from m^2 -> km^2
                                cds.append(cd)
                                pols.append(day1_vort_pol)
                                break # do not need to look at other CDs if the criteria to this point has been met 
            except TypeError: # too close to land to define an RCLV
                break
    return plms,cons,areas,cds,pols

def set_up_RCLV_atlas(date_list):
    """
    Generate the structure for an atlas based on RCLVs identified each date in the date_list
    
    Input
        date_list: dates to search for RCLVs; need to have already produced an LAVD field for this date
    """
    print('Identifying RCLVs...')
    
    # Set up RCLV atlas
    RCLV_data = [['Date', 'RCLV ID', 'Orientation', 'Age (days)', 'Area (km^2)', 'Center Lon', 'Center Lat', 'CD', 'flag', 'Boundary Coords']]
    for date in date_list:
        print(date)

        plms,cons,areas,cds,pols = LAVD_to_RCLV(date)
        for i in np.arange(0,len(cons)):
            temp = [date,np.nan,pols[i],np.nan,areas[i],traj_lon_array[plms[i][1]],traj_lat_array[plms[i][0]],cds[i],0] 
            
            # Convert from grid indices to lat/lon coordinates for the output
            lon_bounds = traj_lon_array[[round(j) for j in cons[i][:, 1]]]
            lat_bounds = traj_lat_array[[round(j) for j in cons[i][:, 0]]]
            for j in np.arange(0,len(lon_bounds)):
                temp.extend([lon_bounds[j],lat_bounds[j]])
            RCLV_data.append(temp) # Add RCLV data to main array
    return RCLV_data
        
####################################### 2. Track the RCLVs through time and ID them #######################################

def track_and_ID_RCLVs(RCLV_data,date_list):
    """
    Track RCLVs through time based on if the center particle in located inside an identified RCLV in the following timestep.
    """
    print('Tracking eddies through time & IDing...')
    
    # QC step: Bug in the peak_local_max function sometimes outputs the same plm twice..
    # This loop ensures we don't get the same contour twice in the dataset
    remove_contour_inds = []
    for date_to_check in date_list:
        # Look for RCLVs on the same date with identical bounds
        RCLV_inds = np.where(np.array([row[0] for row in RCLV_data]) == date_to_check)[0]
        for a, b in itertools.combinations(RCLV_inds, 2): #tries all combinations from the list
            if RCLV_data[a][9:] == RCLV_data[b][9:]:
                remove_contour_inds.append(a)
    print('Number of identical contours to remove: %s'%(len(remove_contour_inds)))
    RCLV_data = np.delete(RCLV_data,remove_contour_inds)
    
    RCLV_id = 1
    for d in np.arange(0,len(date_list)): # Iterate through the dates in order to ID eddies that have not been yet
        date = date_list[d]
        print(date)

        # Open the Lagrangian trajectory data for the current date
        traj = xr.open_dataset('%s%s_%s.nc'%(lag_traj_dir,date,filename_str))
        particle_lon,particle_lat = traj.variables["lon"],traj.variables["lat"]
        
        day_inds = np.where((np.array([row[0] for row in RCLV_data]) == date))[0] # Get indeces of the RCLV data from the current date
        for r in np.arange(0,len(day_inds)): # iterate through RCLVs on the current date         
            row = RCLV_data[day_inds[r]]              
            if np.isnan(row[1]): #only check rows that don't already have an ID
                row[1] = RCLV_id # Assign the RCLV a new ID
                x_ind = int(np.where(traj_lon_array==row[5])[0]) # center point lon ind
                y_ind = int(np.where(traj_lat_array==row[6])[0]) # center point lat ind
                particle_num = x_ind*len(traj_lat_array) + y_ind #identifier for the center particle in the trajectory file

                # Check if the particle is in an RCLV on the days before the search date, but next in the 'date_list' (because date_list is in backwards order)
                if (d < (len(date_list) - 1)): # Can potentially go up until the end of the date list
                    eddy_stop,counter = 0,1
                    while ((eddy_stop == 0) & (counter < (len(date_list)-d))): #stops when eddy is not found on the next date OR after going through all available dates
                        new_date = date_list[d+counter] # this is the date of the data we will look through in RCLV_data (back ~8 days from date)

                        # Track center particle back in time to the previous initialization date
                        if counter == 1: # get trajectory from the open dataset for 1 timestep back
                            days_btwn = num_days_between(new_date,date) # num of days between date & new_date (typically 8 except for when in between years)
                            # Get lat,lon coords of center particle back 1 timestep
                            backxdays_center_lon = float(particle_lon[particle_num,days_btwn*int(24/sim_params['output_freq'])])
                            backxdays_center_lat = float(particle_lat[particle_num,days_btwn*int(24/sim_params['output_freq'])])

                        else: # open traj dataset for previous date (sequentially forward one timestep) to track particle to current date (counter will be > 1)
                            prev_date = date_list[d+counter-1] # forward one time step from new date
                            prev_traj = xr.open_dataset('%s%s_%s.nc'%(lag_traj_dir,prev_date,filename_str))
                            prev_particle_lon,prev_particle_lat = prev_traj.variables["lon"],prev_traj.variables["lat"]

                            # The new 'center position' will not line up perfectly with the grid of the next initialization, so need to use find_nearest()
                            days_btwn = num_days_between(new_date,prev_date)
                            
                            x_ind = int(np.where(traj_lon_array==backxdays_center_lon)[0]) # center point lon ind
                            y_ind = int(np.where(traj_lat_array==backxdays_center_lat)[0])
                            particle_num = x_ind*len(traj_lat_array) + y_ind #identifier for the center particle in prev_traj_file

                            # Replace old variables with the locations tracked back 8 more days
                            backxdays_center_lon = float(prev_particle_lon[particle_num,days_btwn*int(24/sim_params['output_freq'])])
                            backxdays_center_lat = float(prev_particle_lat[particle_num,days_btwn*int(24/sim_params['output_freq'])])

                        # Check if the center particle back 1 timestep is inside another contour; if so, assign the contour the RCLV ID
                        new_center_point = Point(backxdays_center_lon,backxdays_center_lat)
                        new_day_inds = np.where((np.array([irow[0] for irow in RCLV_data]) == new_date))[0]
                        for nr in np.arange(0,len(new_day_inds)): #iterate through all of the RCLVs at the search date   
                            nrow = RCLV_data[new_day_inds[nr]]

                            if (np.isnan(nrow[1])) and (nrow[2] == row[2]): #only want to check RCLVs that don't have an ID & have the same orientation
                                x_bnds, y_bnds = nrow[9:][0::2], nrow[9:][1::2]
                                polygon = Polygon([(x_bnds[pt],y_bnds[pt]) for pt in np.arange(0,len(x_bnds))])
                                if (polygon.contains(new_center_point)): # check if point is inside RCLV boundaries on the previous date
                                    nrow[1] = RCLV_id #change nan to ID
                                    backxdays_center_lon,backxdays_center_lat = nrow[5],nrow[6] # Assign LAVD local max of this contour as the new particle to track 
                                    counter += 1
                                    break #stop look through RCLVs for this date and go to the next (back up to the top of the while loop)

                            if nr == (len(new_day_inds) - 1): #on the last RCLV and still haven't found one that contains the center point                            
                                eddy_stop = 1 
                RCLV_id += 1  
                
    return RCLV_data

####################################### 3. QC: Check if any RCLVs skipped a date (or 2) #######################################

def draw_eddy_boundaries_back_timesteps(log_file,date_list,eddy_id,current_date,traj,poly_lons,poly_lats,center_particle_num,n):
    """
    Input
        log_file: Notes about any issues are outputted here
        eddy_id: assigned from RCLV tracking
        current_date: intialization date to draw boundaries back from 
        traj: trajectory dataset
        poly_lons,poly_lats: coordinates of the polygon
        center_particle_num: particle ID number in the trajectory dataset
        n: number of times to go back a timestep (i.e., to draw boundaries back ~8 days 4 times, n = 4)
            NOTE: there is a limit on doing this; can't do more than length of the simulation
    Output
        successful_dates: dates that a boundary was successfully drawn
        traj_back_timesteps_boundaries: boundary coordinates
        eddy_center_lons,eddy_center_lats: eddy center coordinates 
        areas_km2: area in kilometers of the feature
    """
    # Project the contour on the grid to find all particles inside the eddy
    poly_pts = [(poly_lons[pt],poly_lats[pt]) for pt in np.arange(0,len(poly_lons))] # Reformat lat/lon boundary points
    x_mask,y_mask = find_polygon_pts(poly_pts,traj_lon_array,traj_lat_array)

    # Get the dates of the timesteps
    traj_back_timesteps_dates = []
    current_date_ind = date_list.index(current_date)    
    for step in np.arange(1,n+1):
        if (current_date_ind+step) <= (len(date_list)- 1): # do not want to project back from the earliest time
            traj_back_timesteps_dates.append(date_list[current_date_ind+step])
        
    # Get the coordinates of the particles
    particle_lon,particle_lat = traj.variables["lon"],traj.variables["lat"] # particle lat/lons from the trajectory file
    
    eddy_center_lons,eddy_center_lats,traj_back_timesteps_boundaries,areas_km2 = [],[],[],[]
    successful_dates = []
    for t in traj_back_timesteps_dates: #iterate through each timestep
        days_btwn = num_days_between(t,current_date)
        time_ind = int(24/sim_params['output_freq'])*days_btwn
        
        # All particles in the contour locations at timestep t
        eddy_timestepx_lons,eddy_timestepx_lats,_ = extract_particles_after_time(traj,x_mask,y_mask,traj_lat_array,sim_params,days_btwn)
        pts = np.array([[eddy_timestepx_lons[i],eddy_timestepx_lats[i]] for i in np.arange(0,len(eddy_timestepx_lons))])
        new_edges = stitch_boundaries(alpha_shape(pts,0.1)) # alpha=0.1 to get a contour resolution similar to floater output
        lon_bounds,lat_bounds = [],[]
        for i,j in new_edges[0]:
            lon_bounds.append(round(pts[[i,j],0][0],5))
            lon_bounds.append(round(pts[[i,j],0][1],5))
            lat_bounds.append(round(pts[[i,j],1][0],5))
            lat_bounds.append(round(pts[[i,j],1][1],5))
        
        # Stop drawing contours if the trajectories go outside of the grid domain
        if (np.max(lat_bounds)>traj_lat_array[-1]) or (np.min(lat_bounds)<traj_lat_array[0]) or (np.max(lon_bounds)>traj_lon_array[-1]) or (np.min(lon_bounds)<traj_lon_array[0]):
            log_file.write("RCLV %s went out of bounds on %s\n"%(eddy_id,t))
            break
        
        successful_dates.append(t)
        
        # Center particle location at timestep t
        eddy_center_lons.append(float(particle_lon[center_particle_num,time_ind]))
        eddy_center_lats.append(float(particle_lat[center_particle_num,time_ind]))
        
        areas_km2.append(calc_area_of_stitched_bounds(lon_bounds,lat_bounds,traj_lon_array,traj_lat_array)) #calculate area of new contour
        
        # Reformat boundaries to be appended to the RCLV dataset
        bounds = []
        for j in np.arange(0,len(lon_bounds)):
            bounds.extend([lon_bounds[j],lat_bounds[j]])
        traj_back_timesteps_boundaries.append(bounds)
            
    return successful_dates,traj_back_timesteps_boundaries,eddy_center_lons,eddy_center_lats,areas_km2

def interpolate_skipped_contours(RCLV_data,log_file,date_list):
    """
    Look for a skipped contour over one or two timesteps. This will appear in the dataset as two RCLVs with different IDs, but they are actually the same
    water mass. By allowing skips, we correct for accidentally counting the second instance water mass as a new eddy. 
    """
    
    print('Checking if any RCLVs skipped a date...')
    RCLV_data_no_header = np.array(RCLV_data)[1:]
    all_ids = np.unique([r[1] for r in RCLV_data[1:]])
    
    # Sort eddy ids into dictionaries that contain their first and last date they appear in the dataset
    first_date_dict,last_date_dict = {},{}
    for d in date_list:
        first_date_dict[d],last_date_dict[d] = [],[]
        
    for eddy_id in all_ids:
        eddy_data = RCLV_data_no_header[np.where(np.array([int(row[1]) for row in RCLV_data_no_header]) == int(eddy_id))]
        eddy_dates = [int(row[0]) for row in eddy_data]
        first_date,last_date = min(eddy_dates),max(eddy_dates)
        first_date_dict[str(first_date)].append(eddy_id)
        last_date_dict[str(last_date)].append(eddy_id)
        
    # Checking for 1 or 2 timesteps skipped between RCLV contours
    for d in np.arange(0,len(date_list[:-3])): 
        current_date = date_list[d]
        print(current_date)
        
        traj = xr.open_dataset('%s%s_%s.nc'%(lag_traj_dir,current_date,filename_str)) 
        particle_lon,particle_lat = traj.variables["lon"],traj.variables["lat"]
        
        for n in [2,3]: #check 2 timesteps back, then 3 timesteps back; contours 1 step back have been given ID already in the tracking function
            back_ntimesteps_date,back_nminus1timesteps_date = date_list[d+n],date_list[d+n-1]
            days_btwn_ntimesteps = num_days_between(back_ntimesteps_date,current_date)
            days_btwn_nminus1timesteps = num_days_between(back_nminus1timesteps_date,current_date) # This is used later in code only if n = 3

            # Subset the data by the dates of interest
            current_date_eddy_data = RCLV_data_no_header[np.where(np.array([int(row[0]) for row in RCLV_data_no_header]) == int(current_date))] 
            back_ntimesteps_eddy_data = RCLV_data_no_header[np.where(np.array([int(row[0]) for row in RCLV_data_no_header]) == int(back_ntimesteps_date))]
            
            # Iterate through RCLVs that start on the current_date in the dictionary
            for eddy_id in first_date_dict[current_date]: 
                eddy_data = current_date_eddy_data[np.where(np.array([int(row[1]) for row in current_date_eddy_data]) == int(eddy_id))][0]
                orientation,center_lon,center_lat = eddy_data[2],eddy_data[5],eddy_data[6]
                x_ind = int(np.where(traj_lon_array==float(center_lon))[0])
                y_ind = int(np.where(traj_lat_array==float(center_lat))[0])
                particle_num = x_ind*len(traj_lat_array) + y_ind #identifier for the center particle

                # Get location of center particle back 2 time steps (~16 days) or 3 time steps (~24 days)
                back_ntimesteps_center_lon = float(particle_lon[particle_num,int(24/sim_params['output_freq'])*days_btwn_ntimesteps])
                back_ntimesteps_center_lat = float(particle_lat[particle_num,int(24/sim_params['output_freq'])*days_btwn_ntimesteps])
                new_center_point = Point(back_ntimesteps_center_lon,back_ntimesteps_center_lat)

                # Iterate through the RCLVs that end n timesteps back to check if it is the same water mass
                for back_ntimesteps_eddy_id in last_date_dict[back_ntimesteps_date]: 
                    if back_ntimesteps_eddy_id != eddy_id: # ignore contours already belonging to this eddy
                        back_ntimesteps_eddy_data_to_check = back_ntimesteps_eddy_data[np.where(np.array([int(row[1]) for row in back_ntimesteps_eddy_data]) == int(back_ntimesteps_eddy_id))][0] 
                        if back_ntimesteps_eddy_data_to_check[2] == orientation: #speeds up computation by only checking eddies with the same orientation
                            bnds = back_ntimesteps_eddy_data_to_check[9:]
                            x_bnds, y_bnds = bnds[0::2],bnds[1::2]
                            polygon = Polygon([(x_bnds[pt],y_bnds[pt]) for pt in np.arange(0,len(x_bnds))])

                            # Check if the center point is inside eddy boundaries
                            if (polygon.contains(new_center_point)): # RCLVs are the same water mass; change the IDs to be the same
                                inds_to_change_id = np.where((np.array([int(irow[1]) for irow in RCLV_data[1:]]) == int(back_ntimesteps_eddy_id)))[0] #look for id in main dataset, ignoring header
                                inds_to_change_id = [i+1 for i in inds_to_change_id] #fix index for main dataset WITH header
                                for ind in inds_to_change_id:
                                    RCLV_data[ind][1] = eddy_id #change ID

                                # Remove the "false" ID from dictionaries
                                last_date_dict[back_ntimesteps_date].remove(back_ntimesteps_eddy_id)
                                for key,values in first_date_dict.items():
                                    if back_ntimesteps_eddy_id in values:
                                        values.remove(back_ntimesteps_eddy_id)
                                        first_date_dict[key].append(eddy_id)
                                        
                                # Draw a boundary around the particles for the date that was skipped in between these two eddy contours & add it to the dataset
                                traj_nminus1_timesteps_back_dates,traj_nminus1_timesteps_back_boundaries,traj_nminus1_timesteps_back_center_lons,traj_nminus1_timesteps_back_center_lats,areas_km2 = draw_eddy_boundaries_back_timesteps(log_file,date_list,eddy_id,current_date,traj,eddy_data[9::][0::2],eddy_data[9::][1::2],particle_num,n-1)

                                # Add new data to RCLV_data
                                for k in np.arange(0,len(traj_nminus1_timesteps_back_dates)):
                                        #Format: ['Date', 'RCLV ID', 'Orientation', 'Age (days)', 'Area (km^2)', 'Center Lon', 'Center Lat', 'CD', 'flag', 'Boundary Coords']
                                        #           0        1          2               3             4             5             6            7     8              9
                                    temp = [traj_nminus1_timesteps_back_dates[k],eddy_id,orientation,np.nan,areas_km2[k],
                                            traj_nminus1_timesteps_back_center_lons[k],traj_nminus1_timesteps_back_center_lats[k],np.nan,1] #skip interp flag=1
                                    temp.extend(traj_nminus1_timesteps_back_boundaries[k])
                                    RCLV_data.append(temp)                                    

                                break #stop looking for an eddy match, go back up looking for a match for the next eddy (2 for loops back)
                            
    return RCLV_data
    

####################################### 4. Interpolate eddy bounds for first 3 timesteps #######################################

def interpolate_first_3timesteps(RCLV_data,log_file,date_list):
    """
    Since the RCLVs are identified every 8 days from 32 day simulations, the first appearance of an RCLV will not have a history before it is
    32 days old. To extract the contours from these early days of the eddy birth, I track the particles back in time and draw a contour around
    the particles. 
    """
    
    print('Interpolating eddy bounds for first 3 timesteps...')    
    RCLV_data_no_header = np.array(RCLV_data)[1:]
    all_ids = np.unique([r[1] for r in RCLV_data[1:]])
    for eddy_id in all_ids:

        # Get date & data from the RCLV's earliest appearance
        eddy_data = RCLV_data_no_header[np.where(np.array([int(row[1]) for row in RCLV_data_no_header]) == int(eddy_id))]
        first_date = min([int(row[0]) for row in eddy_data])
        print(first_date) 
        
        first_date_data = eddy_data[np.where(np.array([int(row[0]) for row in eddy_data]) == first_date)][0]
        first_date_lons,first_date_lats = first_date_data[9::2],first_date_data[10::2] #boundary coords
        x_ind = int(np.where(traj_lon_array==first_date_data[5])[0])
        y_ind = int(np.where(traj_lat_array==first_date_data[6])[0])
        particle_num = x_ind*len(traj_lat_array) + y_ind #identifier for the center particle

        # Track the contour back (up to) 3 timesteps
        first_date = str(first_date)
        if first_date != date_list[-1]: #only track backwards if not the first day of the dataset
            if first_date == date_list[-2]: # track backwards 1 timestep
                n = 1
            elif first_date == date_list[-3]: # track backwards 2 timesteps
                n = 2
            else: # track backwards 3 timesteps
                n = 3
            traj = xr.open_dataset('%s%s_%s.nc'%(lag_traj_dir,first_date,filename_str))
            traj_back_timesteps_dates,traj_back_timesteps_boundaries,traj_back_timesteps_center_lon,traj_back_timesteps_center_lat,areas_km2 = draw_eddy_boundaries_back_timesteps(log_file,date_list,eddy_id,first_date,traj,first_date_lons,first_date_lats,particle_num,n)                
            for j in np.arange(0,len(traj_back_timesteps_boundaries)):
                #['Date', 'RCLV ID', 'Orientation', 'Age (days)', 'Area (km^2)', 'Center Lon', 'Center Lat', 'CD', 'flag', 'Boundary Coords']
                #   0        1          2               3             4             5             6            7         8              9
                temp = [traj_back_timesteps_dates[j],eddy_id,first_date_data[2],np.nan,areas_km2[j],
                            traj_back_timesteps_center_lon[j],traj_back_timesteps_center_lat[j],np.nan,2] # flag=2 for first 3 timestep interps
                temp.extend(traj_back_timesteps_boundaries[j])
                RCLV_data.append(temp) # Add data to main RCLV array
    
    return RCLV_data
   
####################################### 5. QC: Address instances of overlapping contours #######################################

def overlapping_RCLV_QC(RCLV_data,log_file,date_list):
    """
    In some rare cases there will be overlapping contours of RCLVs from one of the interpolated boundaries steps. Some example 
    phenomena of what could cause this include eddy splitting, eddy merging, an "escaping" center particle, etc. 
    """
    print('Find instances of overlapping RCLV contours...')

    count3,count4,count5 = 0,0,0
    contain_IDs = {}
    for d in date_list:
        print(d)
            
        RCLV_inds = np.where(np.array([row[0] for row in RCLV_data]) == d)[0]# get all of the data for the RCLVs on the current date (d)     
        for a, b in itertools.combinations(RCLV_inds, 2): #tries all combinations from the list
            ID1,ID2 = RCLV_data[a][1],RCLV_data[b][1]
            bnds1,bnds2 = RCLV_data[a][9:],RCLV_data[b][9:]
            contour_pts1 = Polygon([(float(bnds1[0::2][i]),float(bnds1[1::2][i])) for i in np.arange(0,len(bnds1[0::2]))])
            contour_pts2 = Polygon([(float(bnds2[0::2][i]),float(bnds2[1::2][i])) for i in np.arange(0,len(bnds2[0::2]))])        

            ############ First check if one contour fully contains the other ##############
            if contour_pts1.contains(contour_pts2):
                count5 += 1
                log_file.write('RCLV %s contains %s on %s\n'%(ID1,ID2,d))
                key_string = '%s,%s'%(ID1,ID2)
                if key_string not in contain_IDs: #Dictionary key order matters here: first ID contains the second
                    contain_IDs[key_string] = [d]
                else:
                    contain_IDs[key_string].append(d)
            elif contour_pts2.contains(contour_pts1): 
                count5 += 1
                log_file.write('RCLV %s contains %s on %s\n'%(ID2,ID1,d))
                key_string = '%s,%s'%(ID2,ID1)
                if key_string not in contain_IDs:
                    contain_IDs[key_string] = [d]
                else:
                    contain_IDs[key_string].append(d)

            ############ Next check if there is an intersection ##############
            elif contour_pts1.intersects(contour_pts2): # order does not matter here            
                # Intersection is one polygon
                intersection = contour_pts1.intersection(contour_pts2) # get the shape of the intersection
                if intersection.type == 'Polygon':
                    try:
                        xx, yy = intersection.exterior.coords.xy
                    except:
                        print(intersection.type)
                        print(intersection)
                    lon_bounds,lat_bounds = xx.tolist(),yy.tolist()
                    try:
                        area_intersect = calc_area_of_stitched_bounds(lon_bounds,lat_bounds,traj_lon_array,traj_lat_array)
                    except: #area too small to compute 
                        area_intersect = 0
                # Intersection is several polygons
                elif intersection.type == 'MultiPolygon': 
                    for poly in list(intersection.geoms): # iterate through all of the polygons from the intersection (can also be 1)
                        xx, yy = poly.exterior.coords.xy
                        lon_bounds,lat_bounds = xx.tolist(),yy.tolist()
                        area_intersect = 0
                        try:
                            area_intersect += calc_area_of_stitched_bounds(lon_bounds,lat_bounds,traj_lon_array,traj_lat_array)
                        except: #area too small to compute 
                            pass

                # If the intersection is small (<5% of the area of both polygons), then remove the intersecting polygon from both RCLVs 
                if (area_intersect/RCLV_data[a][4] <= 0.05) and (area_intersect/RCLV_data[b][4] <= 0.05):
                    count3 += 1
                    log_file.write('Removing intersection (<=5 percent of area) between RCLV %s and %s on %s \n'%(ID1,ID2,d))
                    temp_data_dict = {a:[contour_pts1],b:[contour_pts2]}
                    for ind in (a,b):

                        # Adjust the RCLV data boundaries
                        if ind==a:
                            primary_contour,secondary_contour = contour_pts1,contour_pts2
                        elif ind==b:
                            primary_contour,secondary_contour = contour_pts2,contour_pts1

                            
                        intersection = primary_contour.difference(secondary_contour)
                        if intersection.type != 'GeometryCollection': 
                            new_poly_x,new_poly_y = intersection.exterior.coords.xy
                        else: # rare case
                            for poly in list(intersection.geoms):
                                if poly.type == 'LineString':
                                    pass
                                else:
                                    new_poly_x,new_poly_y = poly.exterior.coords.xy
                                
                        reformat_bounds = []
                        for l in np.arange(0,len(new_poly_x)):
                            reformat_bounds.append(new_poly_x[l])
                            reformat_bounds.append(new_poly_y[l])
                        RCLV_data[ind][9:] = reformat_bounds
                        RCLV_data[ind][8] = str(int(RCLV_data[ind][8])) + str(3) #<- 3 flag will mean adjusted boundaries

                # If the intersection is large, flag it and do nothing 
                else:
                    count4 += 1
                    log_file.write('Intersection (>5 percent of area) between RCLV %s and %s on %s \n'%(ID1,ID2,d))
                    for ind in (a,b):
                        RCLV_data[ind][8] = str(int(RCLV_data[ind][8])) + str(4) #<- 4 flag will mean there is an RCLV overlap / split
                        
    #### Combine RCLVs that are fully contained in another into one RCLV ID 
    remove_contour_inds = []
    for key,values in contain_IDs.items():
        ID_to_keep,ID_to_remove = key.split(',')
        ID_to_remove_inds = np.where([row[1] == int(ID_to_remove) for row in RCLV_data])[0]
        for ind in ID_to_remove_inds:
            if RCLV_data[ind][0] in values: #check if this is an overlap date
                remove_contour_inds.append(ind) # remove this data point later to not mess up the indexing
            else: # change ID & flag 
                RCLV_data[ind][1] = ID_to_keep # change the ID 
                RCLV_data[ind][8] = str(int(RCLV_data[ind][8])) + str(5) #<- 5 flag will mean there was an ID change
                
    RCLV_data = np.delete(RCLV_data,remove_contour_inds)
    log_file.write('%s small intersections removed\n'%(count3))
    log_file.write('%s large intersections flagged\n'%(count4))
    log_file.write('%s overlapping RCLV intances collapsed\n'%(count5))
    return RCLV_data
    
            
####################################### 6. Give RCLVs an age #######################################

def age_RCLVs(RCLV_data):
    print('Giving RCLVs an age & saving the final atlas...')

    # Get the unique IDs
    all_ids = np.unique([r[1] for r in RCLV_data[1:]])
    
    # Age each instance of each ID
    for i in all_ids:
        eddy_inds = np.where(np.array([str(row[1]) for row in RCLV_data]) == str(i))[0]        
        eddy_dates = [int(RCLV_data[e][0]) for e in eddy_inds]
        eddy_dates_sorted, eddy_inds_sorted = (list(t) for t in zip(*sorted(zip(eddy_dates, eddy_inds))))
        
        eddy_flags = np.array([int(RCLV_data[e][8]) for e in eddy_inds_sorted])
        ref_ind = np.where((eddy_flags==0)|(eddy_flags==3)|(eddy_flags==4)|(eddy_flags==5))[0][0] #first identified RCLV (32 days old)
        ref_date = str(eddy_dates_sorted[ref_ind])
        
        # Iterate through each instance of this eddy ID
        for j in np.arange(0,len(eddy_inds_sorted)):
            if j == ref_ind:
                RCLV_data[eddy_inds_sorted[j]][3] = sim_params['runtime'] #find the first non-interpolated contour, this is a 32-day old instance
            elif j < ref_ind: #subtract dates that come before reference
                RCLV_data[eddy_inds_sorted[j]][3] = sim_params['runtime'] - num_days_between(ref_date,str(eddy_dates_sorted[j])) 
            else: #add dates that come after reference
                RCLV_data[eddy_inds_sorted[j]][3] = sim_params['runtime'] + num_days_between(ref_date,str(eddy_dates_sorted[j])) 
    return RCLV_data
