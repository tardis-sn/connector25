
# How to use SNEC within the pipeline

[SNEC](https://stellarcollapse.org/index.php/SNEC.html) takes a stellar model as input. The input profile needs to have
the columns `zone`, `logRho`, `logT`, `velocity`, `ye`, `radius`, `mass`, the mass
fractions of all isotopes included, `mixing_type`, and optionally `omega`.


## Generating the input models

The scripts to generate input models are in `./scripts`

The script `mesa_to_GR1D.py` takes a `MESA` profile and converts it to a
`*.short` format, which is the input format for both [GR1D](https://github.com/evanoconnor/GR1D) and `SNEC`. The
path of the input profile and output file (including filenames) are
command line arguments for this script.

The script `MESA_isotopes.py` script will generate from the input `MESA`
profile the file for the composition input of `SNEC`. Once again the
paths of the input and output files (including filenames) are command
line arguments. For untested nuclear networks, the script may need
adjustments, but it should work with many networks out-of-the-box.

For convenience, the script `short_and_iso.py` is a wrapper that will
run both the scripts described above, where the only argument is the
path and filename of the input `MESA` profile.


## Generating the `SNEC` grid

`SNEC` needs to have a file telling it the grid to be used. This
typically is *not* the same as the `MESA` input profile grid. This can be
done with the script `./script/grid_setup.py` which is a modified
version of the script that is part of `SNEC`. It takes as input only the
number of mesh points in the `MESA` profile and it produces a file that
**needs** to be in `SNEC-1.01/tables/GridPattern.dat` (with this specific
filename).


## The `SNEC` parameter file

An example `parameters` file for `SNEC` is in `./input`. Please refer to the
[SNEC documentation (pdf)](https://stellarcollapse.org/codes/snec_notes-1.00.pdf) for a complete description.

For each input `MESA` profile, one needs to hardcode in this file the
number of mesh points in the input profile in the variable `imax` (same
number as the input of `grid_setup.py` described above). Possibly opacity
floors, and thermal bomb or piston properties need to be adjusted too.
