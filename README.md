# RCLVatlas

Lexi Jones (MIT)

Email: jonesae@mit.edu

Date created: 12/08/22

Last edited: 12/09/22


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

4. Practice using the package tools with the Jupyter notebook `example_usage.ipynb`

#########################

