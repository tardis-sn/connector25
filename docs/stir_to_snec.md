# STIR to SNEC conversion

## Inputs to STIR
`STIR` takes as input a stellar model (density, temperature, radius, radial 
velocity, ye, 1/Abar). `STIR` requires an input file `flash.par` which 
configures, e.g., grid options such as radial extent and adaptive mesh 
refinement levels and criteria. This is where the run base name (`basenm`), 
model file, and output directory are set. To run `STIR`, navigate to 
`obj_ccsn1dMLTnoHyb/run` and then `./stir -par_file inputs/stir.par` 
(or, with MPI, e.g.m `mpirun -n 8 ./stir -par_file inputs/stir.par`).

## Outputs from STIR
`STIR` outputs hydrodynamic and radiation radial profiles as a function of time
in the `output` directory. These are HDF5 files. Quantities of interest might 
include density, temperature, pressure, energies, ye, and velocity.
`STIR` outputs integrated time series data in a `.dat` file prefixed with the 
run `basenm` -- these are not in the `outputs` directory. The columns are 
labeled and include quantities such as central density, diagnostic explosion 
energy, mean shock radius, and mass accretion rate as a function of time.
There is also a `.log` file which summarizes options for the current run.

## Inputs to SNEC
