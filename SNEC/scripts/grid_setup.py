"""

 This file generates a pattern used for the gridding of the SNEC models

 Total number of grid points is equal to imax.
 Current pattern consists of two geometric progressions, connected at the
 point itran and ensuring finer resolution of the grid in the inner (explosion)
 region and the surface region of the models.

 ratio1 gives the ratio between the size of the first and itran-th cell
 ratio2 gives the ratio between the size of the last  and itran-th cell

 Modified by Mathieu Renzo 20250321 to make imax an input
"""

from __future__ import division
import math
import sys

def make_SNEC_grid(output_file: str):
    # read from command line
    imax = int(sys.argv[1])  # number of mesh points wanted, needs to match input model
    delta = []
    grid_pattern = []
    itran = 100
    ratio1 = 0.1
    ratio2 = 0.001

    f1 = ratio1 ** (1 / (itran - 2))
    f2 = ratio2 ** (1 / (imax - itran - 1))

    delta_tran = (
        (1 - f2)
        * (1 - f1)
        / ((1 - f2) * (1 - f1 ** (itran - 1)) + (1 - f1) * (1 - f2 ** (imax - itran)))
    )

    for l in range(0, itran - 1):
        delta.append(delta_tran * f1 ** (itran - 2 - l))

    for l in range(0, imax - itran):
        delta.append(delta_tran * f2 ** (l))

    grid_pattern.append(0.0)
    for l in range(0, imax - 1):
        grid_pattern.append(grid_pattern[l] + delta[l])

    print("grid pattern consists of", len(grid_pattern), "points")

    with open(output_file, "w") as outfile:
        for l in range(0, imax):
            outfile.write(str(grid_pattern[l]) + "\n")

    print("done!")



if __name__ == "__main__":
    make_SNEC_grid("SNEC-1.01/table/sGridPattern.dat")
