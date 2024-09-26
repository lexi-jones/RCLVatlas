# RCLVatlas

Lexi Jones-Kellett (MIT)

Email: jonesae@mit.edu


## About

This is the supplemental code for ESSD Manucript Jones-Kellett & Follows 2024 (https://doi.org/10.5194/essd-16-1475-2024).


This code builds off of the core packages described below (OceanParcels & floater) to create an atlas of Rotationally Coherent Lagrangian Vortices (RCLVs) from the Lagrangian Averaged Vorticity Deviation (LAVD; Haller et al 2016 - doi:10.1017/jfm.2016.151). The pipeline includes functions to initialize particles in a grid formation and advect them through satellite geostrophic velocity fields distributed from CMEMS (https://doi.org/10.48670/moi-00148). We built a custom kernel to calculate the relative vorticity along a Lagrangian particle trajectory, which is necessary for calculating the LAVD.

The LAVD is calculated by integrating the vorticity along a particle trajectory, and subtracting the domain average vorticity. Each particle will then have a single LAVD value assigned to it. By plotting the LAVD at the particle initialization location, RCLVs are identified from closed contours surrounding local maxima in the LAVD field. The floater package is used to identify the appropriate contours in the LAVD field. 

We build on the floater package by iterating through viable convex deficiency parameters, requiring minimal dispersal, and ensuring that the sign of the vorticity is consistent at the beginning and end of the feature lifetime. Finally, this code can be used to track RCLVs through time and space by assigning water masses a unique ID. The tracking methodology allows one to create an RCLV atlas that can be easily compared with Eulerian eddy atlases. 

## Core Packages

1. OceanParcels v2.2 (Delandmeter and van Sebille 2019 - https://doi.org/10.5194/gmd-12-3571-2019)
	- Package documentation: https://oceanparcels.org/index.html
	- Used to run Lagrangian trajectories
	- Note that Version 3 of Parcels outputs zarr files instead of netCDF files

2. floater (Tarshish et al. 2018 - https://doi.org/10.1016/j.ocemod.2018.07.001)
	- Package documentation: https://floater.readthedocs.io/en/latest/
	- Used to identify RCLVs from the LAVD field

## Installation

1. Install conda following the instructions here: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
	- Select Python 3

2. Set up environment with the YML file (preferred method with correct versioning):
	- `conda env create -f py3_parcels_v2.yml`

~~~OR~~~
2. From the terminal (Linux/macOS) or Anaconda prompt (Windows), create a new environment with the dependencies needed to use this package:
	- `conda create -n myenvname python`
	- `conda install -c conda-forge parcels jupyterlab numpy xarray scipy matplotlib shapely`
	- `conda install scikit-image`
   NOTE: Need parcels v2.2 to output netCDF files
~~~~~~~~

3. Clone the floater repository
	- `git clone https://github.com/ocean-transport/floater.git`

4. Clone the repository
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

## Citation

If you use this software, please cite the following:

Jones-Kellett, A. (2023). RCLVatlas (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.7702978

Jones-Kellett, AE & Follows, MJ. (2024). A Lagrangian coherent eddy atlas for biogeochemical applications in the North Pacific Subtropical Gyre. Earth Syst. Sci. Data. 16, 1475–1501. https://doi.org/10.5194/essd-16-1475-2024


