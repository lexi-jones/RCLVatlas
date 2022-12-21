# RCLVatlas

Lexi Jones (MIT)

Email: jonesae@mit.edu

Date created: 12/08/22

Last edited: 12/20/22


## About


## Key Packages

1. OceanParcels (Delandmeter and van Sebille 2019 - https://doi.org/10.5194/gmd-12-3571-2019)
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

2. Download CMEMS daily geostrophic velocity data
    - https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_MY_008_047/description
    - Code is setup to run on geostrophic velocity data with each day stored in a seperate file of format: 'dt_global_allsat_phy_l4_YYYYMMDD.nc'

3. Define local directory paths and adjust parameters with the `config.py` file

4. Run Lagrangian particle simulation with command `python run_parcels_CMEMS.py YYYYMMDD`
    - The geostrophic velocity data must be downloaded for each day needed for the timeframe of your particle simulation
    
5. Practice using the RCLV tracking tools with the Jupyter notebook `example_usage.ipynb`




