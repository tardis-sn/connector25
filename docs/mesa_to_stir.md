# MESA to STIR conversion

## Inputs to MESA

MESA input consists of a `work` directory which includes project `inlists` (user configuration files), a `run_star_extras/` folder for custom Fortran functions that hook into the MESA stellar evolution solve, and make/run/restart bash scripts.
To carry out a MESA simulation locally on your laptop, compile with `./mk` and run with `./rn`.

## Outputs from MESA

MESA outputs include history and profile data (will appear in the `LOGS` directory in your `work` folder), photos for restarting the simulation, and images (optional). We include output files for our fiducial mesa methods. You should be able to reproduce this output by running the inlists as they appear in the inputs section. 

## Inputs to STIR

## Conversion from MESA to STIR

To convert the MESA output to STIR input, ...

## Outputs from STIR

