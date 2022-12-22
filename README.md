# RCLVatlas

Lexi Jones (MIT)

Email: jonesae@mit.edu


## About

This code builds off of the core packages (OceanParcels & floater) to create an atlas of Rotationally Coherent Lagrangain Vortices (RCLVs) from the Lagrangian Averaged Vorticity Deviation (LAVD - Haller et al 2016 - doi:10.1017/jfm.2016.151). The pipeline includes functions to initialize particles in a grid formation, and advect them through satellite geostrophic velocity fields from CMEMS. Moreover, we built a custom kernel to calculate the vorticity along a Lagrangian particle trajectory, which is necessary for calculating the LAVD.

The LAVD is calculated by integrating the vorticity along a particles trajectory, and subtracting the domain average vorticity. Each particle will then have a single LAVD value assigned to it. By plotting the LAVD at the particle initialization location, RCLVs are easily identifiable as circular local maxima in the LAVD field. The floater package is used to identify closed contours around these local maxima in order to define the boundary of a coherent structure. 

We build on the floater package by iterating through viable convex deficiency parameters, requiring minimal dispersal, and ensuring that the sign of the vorticity is consistent and the beginning and end of the features lifetime. Finally, this code can be used to track RCLVs through time and space by assigning water masses a unique ID. The tracking methodology allows one to create an RCLV atlas that can be easily compared with Eulerian eddy atlases. 

## Core Packages

1. OceanParcels v2.2 (Delandmeter and van Sebille 2019 - https://doi.org/10.5194/gmd-12-3571-2019)
	- Package documentation: https://oceanparcels.org/index.html
	- Used to run Lagrangian trajectories

2. floater (Tarshish et al. 2018 - https://doi.org/10.1016/j.ocemod.2018.07.001)
	- Package documentation: https://floater.readthedocs.io/en/latest/
	- Used to identify RCLVs from the LAVD field

## Installation

1. Install conda following the instructions here: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
	- Select Python 3

2. From the terminal (Linux/macOS) or Anaconda prompt (Windows), create a new environment with the dependencies needed to use this package:
	- `conda create -n myenvname python`
	- `conda install -c conda-forge parcels jupyterlab numpy xarray scipy matplotlib shapely`
	- `conda install scikit-image`

3. Clone the floater repository
	- `git clone https://github.com/ocean-transport/floater.git`

3. Clone the repository
	- `git clone https://github.com/lexi-jones/RCLVatlas.git`

## Usage

1. Activate the conda environment
    `source activate myenvname`

2. Download CMEMS daily geostrophic velocity data (ftp recommended)
    - https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_MY_008_047/description
    - Code is setup to run on geostrophic velocity data with each day stored in a seperate file of format: 'dt_global_allsat_phy_l4_YYYYMMDD.nc'

3. Define local directory paths and adjust parameters with the `config.py` file

4. Reformat the CMEMS longitude array and file name with `reformat_CMEMS_longitude.py`
    - Longitude coordinates need to be 0 -> 360 (default is -180 to 180 when downloaded from CMEMS)

4. Run Lagrangian particle simulation & calculate LAVD with command `python run_parcels_CMEMS.py YYYYMMDD`
    - The geostrophic velocity data must be downloaded for each day needed for the timeframe of your particle simulation
    
5. Practice using the RCLV tracking tools with the Jupyter notebook `example_usage.ipynb`; `produce_RCLV_atlas.py` can be used
   for producing larger datasets 
