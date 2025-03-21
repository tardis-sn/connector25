# MESA to STIR conversion

## Inputs to MESA

MESA input consists of a `work` directory which includes project `inlists` (user configuration files), a `run_star_extras/` folder for custom Fortran functions that hook into the MESA stellar evolution solve, and make/run/restart bash scripts.
To carry out a MESA simulation locally on your laptop, compile with `./mk` and run with `./rn`.

## Outputs from MESA

MESA outputs history and profile data in the `LOGS` directory in your `work` folder, photos for restarting the simulation, and images (optional).

## Inputs to STIR

## Conversion from MESA to STIR

To convert the MESA output to STIR input, ...

## Outputs from STIR

