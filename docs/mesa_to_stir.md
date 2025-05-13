# MESA to STIR conversion

## Inputs to MESA

MESA input consists of a `work` directory which includes project `inlists` (user configuration files), a `run_star_extras/` folder for custom Fortran functions that hook into the MESA stellar evolution solve, and make/run/restart bash scripts.
To carry out a MESA simulation locally on your laptop, compile with `./mk` and run with `./rn`.

## Outputs from MESA

MESA outputs include history and profile data (will appear in the `LOGS` directory in your `work` folder), photos for restarting the simulation, and images (optional). We include output files for our fiducial mesa methods. You should be able to reproduce this output by running the inlists as they appear in the inputs section. 

## Inputs to STIR

`STIR` takes as input a stellar model (density, temperature, radius, radial 
velocity, ye, 1/Abar). `STIR` requires an input file `flash.par` which 
configures, e.g., grid options such as radial extent and adaptive mesh 
refinement levels and criteria. This is where the run base name (`basenm`), 
model file, and output directory are set. To run `STIR`, navigate to 
`obj_ccsn1dMLTnoHyb/run` and then `./stir -par_file inputs/stir.par` 
(or, with MPI, e.g.m `mpirun -n 8 ./stir -par_file inputs/stir.par`).

## Conversion from MESA to STIR

To convert the MESA output to STIR input, we take a final MESA profile 
at the onset of core-collapse (roughly 1000 km/s maximum infall) and 
translate this into a form suitable for `STIR`.
There is a script.... (#TODO)

## Outputs from STIR
`STIR` outputs hydrodynamic and radiation radial profiles as a function of time
in the `output` directory. These are HDF5 files. Quantities of interest might 
include density, temperature, pressure, energies, ye, and velocity.
`STIR` outputs integrated time series data in a `.dat` file prefixed with the 
run `basenm` -- these are not in the `outputs` directory. The columns are 
labeled and include quantities such as central density, diagnostic explosion 
energy, mean shock radius, and mass accretion rate as a function of time.
There is also a `.log` file which summarizes options for the current run.
