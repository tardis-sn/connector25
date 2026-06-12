# SNEC + fallback

MESA-explosion can provide velocity profiles which include fall back
matter with negative velocity.
[SNEC](https://stellarcollapse.org/index.php/SNEC.html) default inner
boundary condition is `v(1)=0` which can cause discrepancies in the
inner velocity profile and artifacts in the dynamics.

For this reason, we employ a modified version of `SNEC` including the
possibility of fallback (without backreaction on the outflows),
`SNEC+fb`. This is available from
[gitHub](https://github.com/mathren/SNEC/tree/tardis_fallback), and
described in Renzo et al., in prep.

The fallback is implemented by implementing a no-gradients inner
boundary and taking `v(1)` to be the minimum between `0` and `v(2)` if
the latter is negative.

Hereafter, we interchangeably refer to `SNEC` or `SNEC+fb`.

# How to use SNEC within the pipeline

`SNEC` takes a stellar model as input. The input profile needs to have
the columns `zone`, `logRho`, `logT`, `velocity`, `ye`, `radius`,
`mass`, the mass fractions of all isotopes included, `mixing_type`,
and optionally `omega`.

It also needs to have the isotopes tracked by the nuclear reaction
network in the stellar model. These will be bundled into categories,
and depending on the species included adjustment may be necessary.

## Generating the input models

The scripts to generate input models are in `./scripts`

The script `mesa_to_GR1D.py` takes a `MESA` profile and converts it to
a `*.short` format, which is the input format for both
[GR1D](https://github.com/evanoconnor/GR1D) and `SNEC`. The path of
the input profile and output file (including filenames) are command
line arguments for this script. You will need the output file path for
the `SNEC` configuration below.

The script `MESA_isotopes.py` will generate from the input
`MESA` profile the file for the composition input of `SNEC`, by
grouping isotopes together in 15 specific groups. Once again, the
paths of the input and output files (including filenames) are command
line arguments. The script should work with `approx21.net` and
`mesa_206.net` but may need adjustments for other nuclear reaction
networks in `MESA`.

*N.B.:* `ni56` needs to be logged in a separate (Z,A) pair in
`isolist`, and `SNEC` will pick it up. All other isotopes are assigned
to the previous group of elements.

For convenience, the script `short_and_iso.py` is a wrapper that will
run both the scripts described above, where the only argument is the
path and filename of the input `MESA` profile.


### Generating the `SNEC` grid

`SNEC` needs to have a file telling it the grid to be used. This
typically is *not* the same as the `MESA` input profile grid. This can be
done with the script `./script/grid_setup.py` which is a modified
version of the script that is part of `SNEC`. It takes as input only the
number of mesh points in the `MESA` profile and it produces a file
`GridPattern.dat`.

*N.B.:* The original `SNEC` **needs** this file to be in
`SNEC-1.01/tables/GridPattern.dat`, `SNEC+fb` allows for a custom
location (see below). The file name is hardcoded in both cases.

*N.B.:* The number of mesh points in a MESA profile can be found in
the header under `num_zones`.

## The `SNEC` parameter file

An example `parameters` file for `SNEC` is in `./input`. Please refer
to the [SNEC documentation
(pdf)](https://stellarcollapse.org/codes/snec_notes-1.00.pdf) for a
complete description. We highlight the most important modifications
below.

### File location

For each input `MESA` profile, one needs to hardcode in `parameters` the
number of mesh points in the input profile in the variable `imax` (same
number as the input of `grid_setup.py` described above). 

`SNEC+fb` allows for customized location of `GridPattern.dat`
providing a `grid_pattern_name`. Empty string will resume default
behavior.

```
# These must exist
outdir              = "./output/SNEC-tardis/"  
profile_name 		= "/path/to/profile.data.short"
comp_profile_name	= "/path/to/profile.data.iso.dat"
grid_pattern_name   = "/path/to/GridPattern.dat"
```

## Allow for fallback with `SNEC+fb`

Assuming you have installed `SNEC+fb`, use the inflow inner boundary
condition.

```
innerBC                = "inflow"
```

(Empty string resumes default boundary conditions, not implemented
without radiation).

### Explosion energy

This was provided by `STIR`, and used in `MESA-explosion`. `SNEC`
should not input more energy. Use

```
initial_data = "Thermal_Bomb"
bomb_mode	 = 2       # do not subtract binding energy from bomb
final_energy = 0.0d40  # do not inject energy
bomb_tstart	 = 0.0d0
```

### Mixing

This has been done by MESA-explosion using a 1D diffusive
approximation for RTI. Turn off `boxcar_smoothing`.

### Nickel heating 

If `ni56` is present in the model, it is logged as its own isotope
class (no other isotopes are included). Then SNEC if
`parameters` is correctly set, `SNEC` can pick it up and calculate
it's contribution to the energetics. For example:

```
	Ni_by_hand = 0  # do not manually add nickel -- it is read from MESA profile instead
	Ni_switch  = 1  # 0 = no nickel heating, 1 = mass within Ni_boundary_mass (will be smoothed by boxcar)
	Ni_mass	   = 0.0d0
```


## Running SNEC

If everything is in the right place, type `./snec` and see the
explosion being simulated!
