# STIR to STELLA conversion

## Inputs to STIR

## Outputs from STIR

## Inputs to STELLA

A simplified version of the STELLA code is packaged in MESA in `$MESA_DIR/stella/`.
STELLA is a 1D multi-group (frequency-dependent) Radiation-
Hydrodynamics code which is useful for modeling bolometric
lightcurves and other observables for core-collapse supernovae.

STELLA takes a `mesa.hyd` (contains hydrodynamic data) and `mesa.abd` (contains abundance data) file as input.

To run STELLA, go to your STELLA folder `cd $MESA_DIR/stella` and make sure the input files (`mesa.hyd` and `mesa.abd`) are placed in `stella/modmake/`. From the stella directory, run just like in MESA, via the command `./rn`

## Conversion from STIR to STELLA

To convert the STIR output to STELLA input, ...

## Outputs from STELLA

Output is contained in `stella/res/`. This includes a `mesa.lbol` lightcurve file, among others.