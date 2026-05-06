# gcp-chem-sim-private
# Oslofjord simulation with FjordSim

## Data preparation

In the preprocess folder there are Python Jupyter notebooks to make:
- A grid file;
- A forcing file for boundary conditions, rivers, etc.;
- An atmospheric forcing file can be prepared using <https://github.com/limash/atm-forcing.git>.

## Installation

Python:
1. Install conda <https://conda-forge.org/download/>.
2. Create a conda environment `conda create --name oslofjord python=3.12`
3. Activate the envoronment `conda activate oslofjord`
4. Navigate to the oslofjord directory
5. Install the dependencies `pip install -e .` (note - some dependencies will not install 
      successfully on an Intel mac, so before running this step, install them with:
      `mamba install -c conda-forge numpy pandas matplotlib seaborn jupyterlab xarray dask netcdf4 bottleneck numba llvmlite scipy h5netcdf zarr cftime cartopy pyproj shapely`
      )
6. For some cases xesmf is required, see <https://xesmf.readthedocs.io/en/stable/installation.html>;
   install with `conda install -c conda-forge xesmf`.
7. Set your base directory as an environment variable in conda `conda env config vars set PROJECT_ROOT=/Path/to/this/repo`

Julia:

1. Install julia <https://julialang.org/downloads/>.
2. Clone the repository: `git clone https://github.com/NIVANorge/oslofjord-sim.git` 
3. Run Julia REPL from the directory with `Project.toml` and activate the environment: `julia --project`.
4. Enter the Pkg REPL by pressing `]` from Julia REPL.
5. Type `instantiate` to 'resolve' a `Manifest.toml` from a `Project.toml` to install and precompile dependency packages.
(you may need to `add https://github.com/NIVANorge/FjordSim.jl.git` from Pkg REPL to install the latest FjordSim).

One of the options to run a FjordsSim simulation requires 2 files:

- A bathymetry netcdf file.
It should contain a 2d array variable "h" with depths (they should be negative),
Two 1d arrays with "lat" and "lon" corresponding to depths in variable "h",
a 1d array "z_faces" with the desired layer depths (also negative values).

- A forcing netcdf file.
This file contains the information about the forcing fields.
To 'force' any variable, one need to define two 4d arrays called, for example, "T" and "T_lambda".
"T" is an oceananigans name for temperature, lets use it further as example;
it is possible to provide forcing for any variable defined in an oceananigans simulation.
Spatial dimensions should have the shape of the corresponding [stagerred grid](https://clima.github.io/OceananigansDocumentation/stable/fields/#Staggered-grids-and-field-locations).
The forth dimension is time in seconds, in Python one can use a datetime format.
"T_lambda" defines a type of forcing to be used in a simulation.
"T" value should correspond to "T_lambda".
If 0 < "T_lambda" < 1, [relaxation](https://clima.github.io/OceananigansDocumentation/stable/model_setup/forcing_functions/#Relaxation) is used.
If "T_lambda" > 1, horizontal flux in "T" should be provided.
If "T_lambda" < -1, vertical flux in "T" should be provided.
there are 2 options:
1. Download the prepared in advance files (`bathymetry_105to232.nc, forcing_105to232.nc, JRA55 files or NORA3.nc`) from [here](https://www.dropbox.com/scl/fo/gc3yc155b5eohi7998wgh/AGN2Yt3HyQ0LlZGImpcca6o?rlkey=x6okc3uxe2avud6sbxgd00l14&st=093llyqp&dl=0) to run a simulation.
2. Use scripts in the preprocess folder to download and prepare the bathymetry and forcing files.
In this case you can add rivers, other sinks and sources, change other forcing for any variable.

## Scenarios

Simulations are configured through **named scenarios** defined in `scenarios.json`. Each scenario maps to a set of input files (bathymetry grid, forcing, atmosphere, and optionally a hotstart file):

```json
{
  "example": {
    "grid_file": "bathymetry_105to232.nc",
    "forcing_file": "forcing_105to232.nc",
    "atmospheric_forcing_dir": "JRA55",
    "hotstart_file": "snapshots_ocean.nc"
  }
}
```

### Adding a new scenario

1. Ensure the files are in the /data/input
2. Add an entry to `scenarios.json` referencing those filenames
3. Commit and push

## Running a simulation

### Locally

e.g.
```bash
export PROJECT_ROOT=/path/to/repo
julia --project simulation.jl --scenario example --arch cpu --stop_days 1
```
