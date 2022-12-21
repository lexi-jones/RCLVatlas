# Functions for the RCLV atlas
#
# Lexi Jones
# Last edited: 12/13/22

import csv,datetime,math
import xarray as xr
import numpy as np
from floater import rclv
from matplotlib.path import Path
from scipy.spatial import Delaunay

####################################### set_up_RCLV_atlas() functions #######################################

def find_polygon_pts(poly_pts,traj_lon_array,traj_lat_array):
    """
    Extract the points from a grid that are inside the polygon.
    
    Input
        poly_pts: 
        grid_points: grid set up with dimensions of the longitude & latitude arrays
    Output
        x_mask,y_mask: indeces from the longitude & latitude arrays that are inside of the polygon
    """
    poly = Path(poly_pts) # make a polygon
    x,y = np.meshgrid(traj_lon_array,traj_lat_array) #convert lat/lon arrays to a grid
    x,y = x.flatten(), y.flatten() #flatten the grid points
    grid_points = np.vstack((x,y)).T # vertical stack the grid points
    grid = poly.contains_points(grid_points) #find what grid points are inside the polygon
    mask = grid.reshape(len(traj_lat_array), len(traj_lon_array)) # now you have a mask with points inside a polygon
    x_mask,y_mask = np.where(mask == True)[1], np.where(mask == True)[0]    
    return x_mask,y_mask

def extract_particles_after_time(traj,x_mask,y_mask,traj_lat_array,sim_params,days):
    """
    Get the lat/lons of particles from a ploygon after some number of days along the Lagrangian trajectory.
    
    Input
        traj: trajectory file 
        x_mask,y_mask: indeces from the longitude & latitude arrays that are inside of the polygon
        days: number of days from the initialization time to retreive particle locations (back trajectories will be back in time)
    Output
        eddy_xdays_lons,eddy_xday_lats: lon/lat coordinates of the particles of interest on day x
    
    """
    particle_lon,particle_lat,particle_vort = traj.variables["lon"],traj.variables["lat"],traj.variables["vort"] #read in particle location lat, lons, vorts
    particle_nums = x_mask*len(traj_lat_array) + y_mask #formula to get the integer ID of the particles inside a polygon
    eddy_xday_lons = [float(particle_lon[p,int((24/sim_params['output_freq'])*days)]) for p in particle_nums]   
    eddy_xday_lats = [float(particle_lat[p,int((24/sim_params['output_freq'])*days)]) for p in particle_nums]
    eddy_xday_vorts = [float(particle_vort[p,int((24/sim_params['output_freq'])*days)]) for p in particle_nums]
    return eddy_xday_lons,eddy_xday_lats,eddy_xday_vorts

def distance_from_lat_lon(lat1,lon1,lat2,lon2):
    """
    Returns the distance between two geographic coordinate points. Accepts negative (-180 to 180) or positive coordinate systems (0 to 360). 

    Input
        lat1,lon1: geographic coordinates for point 1
        lat2,lon2: geographic coordinates for point 2
    Output
        dist: distance between the two point (units of kilometers)
    """
    R = 6371 # Radius of the earth in km
    delta_lat,delta_lon = math.radians(lat2-lat1),math.radians(lon2-lon1)
    lat1_radians,lat2_radians = math.radians(lat1),math.radians(lat2)
    a = math.sin(delta_lat/2) * math.sin(delta_lat/2) + math.cos(lat1_radians) * math.cos(lat2_radians) * ((math.sin(delta_lon/2))**2)
    dist = R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) # Distance in km
    return dist

def positional_variance(particle_lons,particle_lats):
    """
    Returns the mean variance in distance of the particles from the mean 
    
    Input
        particle_lons,particle_lats: geographic coordinates
    """
    mean_lon,mean_lat = np.mean(particle_lons),np.mean(particle_lats)
    
    dist = []
    for p in np.arange(0,len(particle_lons)):
        lon,lat = particle_lons[p],particle_lats[p]
        dist.append(distance_from_lat_lon(mean_lat,mean_lon,lat,lon)**2)
    return np.mean(dist)

def CI(particle_lons_t0,particle_lats_t0,particle_lons_tx,particle_lats_tx):
    """
    Calculates the coherency index after x days
    
    Input
        particle_lons_t0,particle_lats_t0: geographic coordinates of Lagrangian particles at t=0
        particle_lons_tx,particle_lats_tx: geographic coordinates of Lagrangian particles at t=x
    Ouput
        CI: coherency index of particles after x days; a unitless measure of dispersion
    """
    variance_t0 = positional_variance(particle_lons_t0,particle_lats_t0)
    variance_tx = positional_variance(particle_lons_tx,particle_lats_tx)
    return ((variance_t0 - variance_tx)/variance_t0)

def save_RCLV_CSV(RCLV_data,output_file):
    """
    Save RCLV data to a CSV file.
    
    Input
        RCLV_data: data matrix
        filename: filename for the RCLV data file
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(RCLV_data)
        
####################################### track_and_ID_RCLVs()functions #######################################

def find_nearest(array, value):
    """
    Find nearest item in the array to the input value
    
    Input
        array: array to find the value in
        value: value of interest
    Output:
        idx: index of the value in the array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def num_days_between(earlier_date,later_date):
    """
    Convert string dates to datetime dates to calculate the number of days between timesteps
    
    Input
        earlier_date: sequentially first date
        later_date: sequentially second date
    """
    d0 = datetime.date(int(earlier_date[0:4]),int(earlier_date[4:6]),int(earlier_date[6:8]))
    d1 = datetime.date(int(later_date[0:4]),int(later_date[4:6]),int(later_date[6:8]))
    delta = d1 - d0
    return abs(delta.days)

def read_RCLV_CSV_untracked(output_file,header_flag):
    """
    Read in the CSV from the output of the `set_up_RCLV_atlas()` function
    
    Input
        filename: filename of the RCLV dataset to open
        header_flag: 1 to read in header; 0 to not read in header
    Output
        RCLV_data: 2D array of format ['Date', 'RCLV ID', 'Orientation', 'Age (days)', 'Area (km^2)', 'Center Lon', 'Center Lat', 'CD', 'interp_flag', 'Boundary Coords']
    """
    RCLV_data = []
    with open(output_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        i = 0
        for row in csv_reader:
            if i == 0:
                if header_flag == 0: # no header
                    pass
                else: #include header
                    RCLV_data.append(row)
            else:
                #reformat so that everything except for the date and orientation is a float
                RCLV_data.append([row[0],float(row[1]),row[2]] + [float(i) for i in row[3:]])
            i += 1
    return RCLV_data

####################################### interpolate_skipped_contours() functions #######################################

def alpha_shape(points, alpha, only_outer=True):
    """
    ADAPTED FROM https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points
    
    Compute the alpha shape (concave hull) of a set of points.
    
    Input
        points: np.array of shape (n,2) points.
        alpha: alpha value (a > 0), the larger the alpha value, the more smooth the boundary is
        only_outer: boolean value to specify if we keep only the outer border or also inner edges.
    Output
        set of (i,j) pairs representing edges of the alpha-shape. (i,j) are the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        """
        Add an edge between the i-th and j-th points, if not in the list already
        """
        if (i, j) in edges or (j, i) in edges: # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer: # if both neighboring triangles are in shape, it's not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    # Loop over triangles: ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa,pb,pc = points[ia],points[ib],points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges

def find_edges_with(i, edge_set):
    """ 
    ADAPTED FROM https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points
    """
    i_first = [j for (x,j) in edge_set if x==i]
    i_second = [j for (j,x) in edge_set if x==i]
    return i_first,i_second

def stitch_boundaries(edges):
    """
    ADAPTED FROM https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points
    
    Stitches together the edges in the appropriate order.
    """
    edge_set = edges.copy()
    boundary_lst = []
    while len(edge_set) > 0:
        edge0 = edge_set.pop()
        boundary = [edge0]
        last_edge = edge0
        while len(edge_set) > 0:
            i,j = last_edge
            j_first, j_second = find_edges_with(j, edge_set)
            if j_first:
                edge_set.remove((j, j_first[0]))
                edge_with_j = (j, j_first[0])
                boundary.append(edge_with_j)
                last_edge = edge_with_j
            elif j_second:
                edge_set.remove((j_second[0], j))
                edge_with_j = (j, j_second[0])  # flip edge rep
                boundary.append(edge_with_j)
                last_edge = edge_with_j
            if edge0[0] == last_edge[1]:
                break
        boundary_lst.append(boundary)
    return boundary_lst

def calc_area_of_stitched_bounds(lon_bounds,lat_bounds,traj_lon_array,traj_lat_array):
    """
    Use floater package to calculate the area of the interpolated contour. This is keeping 
    consistent with the area calculated for the non-interpolated contours.
    
    Input
        lon_bounds,lat_bounds: lat/lon coords of RCLV boundary
        traj_lon_array,traj_lat_array: array of the grid coordinates
    Output
        area of contour (km^2)
    """
    interp_con = []
    for vert in np.arange(0,len(lat_bounds)): #skip every other because the vertices repeat after stitching together
        lon_index = find_nearest(traj_lon_array,lon_bounds[vert])
        lat_index = find_nearest(traj_lat_array,lat_bounds[vert])
        interp_con.append([lat_index,lon_index])
        
    args = {'lon0':np.max(lon_bounds)-np.min(lon_bounds),
            'lat0':np.max(lat_bounds)-np.min(lat_bounds),
            'dlon':np.abs(traj_lon_array[1]-traj_lon_array[0]),
            'dlat':np.abs(traj_lat_array[1]-traj_lat_array[0])}

    region_area = rclv.polygon_area(rclv.project_vertices(np.array(interp_con),**args)) 
    return region_area/(10**6) # convert from m^2 to km^2

def read_RCLV_CSV_tracked(output_file):
    """
    This function is specific to RCLV atlas files that have IDs & orientation already attributed
    """
    RCLV_data = []
    with open(output_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        i = 0
        for row in csv_reader:
            if i == 0: # header
                row_reformat = row
            else: # reformat so that everything except for the date, ID, and orientation is a float
                row_reformat = [row[0],int(row[1]),row[2]] + [float(i) for i in row[3:]] 
            RCLV_data.append(row_reformat)
            i += 1
    return RCLV_data