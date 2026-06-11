#!/usr/bin/env python

import sys
import re
import numpy as np

msun = 1.99e33
rsun = 6.96e10

ELEM_Z = {
    'h': 1, 'he': 2, 'li': 3, 'be': 4, 'b': 5,
    'c': 6, 'n': 7, 'o': 8, 'f': 9, 'ne': 10,
    'na': 11, 'mg': 12, 'al': 13, 'si': 14, 'p': 15,
    's': 16, 'cl': 17, 'ar': 18, 'k': 19, 'ca': 20,
    'sc': 21, 'ti': 22, 'v': 23, 'cr': 24, 'mn': 25,
    'fe': 26, 'co': 27, 'ni': 28, 'cu': 29, 'zn': 30,
}

GROUP_Z = [2,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]
GROUP_A = [4, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56]
N_GROUPS = len(GROUP_Z)  # 13

# Inclusive upper-Z boundary for groups 0-11; group 12 is ni56-only
_UPPER_Z = [4, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 10**6]

NON_COMP_COLS = {'mass', 'radius', 'mixing_type', 'h1', 'neut'}
_ISOTOPE_RE = re.compile(r'^([a-z]+)\d+$')


def isotope_group(name):
    """Return group index 0-12, or None if not a bundleable isotope."""
    if name == 'ni56':
        return N_GROUPS - 1          # reserved last group
    m = re.match(r'^([a-z]+)(\d+)$', name)
    if not m:
        return None
    z = ELEM_Z.get(m.group(1))
    if z is None:
        return None
    for g, upper in enumerate(_UPPER_Z):
        if z <= upper:
            return g
    return None


def main():
    path, pathout = sys.argv[1], sys.argv[2]

    with open(path) as f:
        lines = f.readlines()

    hdr = [lines[i].split() for i in range(3)]
    zones = int(hdr[2][hdr[1].index("num_zones")])

    col_names = lines[5].split()
    arr = np.loadtxt(path, skiprows=6)

    imass   = col_names.index("mass")
    iradius = col_names.index("radius")
    ineut   = col_names.index("neut")
    iprot   = col_names.index("h1")

    # Accumulate mass fractions into groups
    group_fracs = np.zeros((N_GROUPS, zones))
    unassigned = []

    for ci, name in enumerate(col_names):
        if name in NON_COMP_COLS:
            continue
        m = _ISOTOPE_RE.match(name)
        if not m or m.group(1) not in ELEM_Z:
            continue
        g = isotope_group(name)
        if g is not None:
            group_fracs[g] += arr[:, ci]
        else:
            unassigned.append(name)

    assert not unassigned, f"Unassigned isotopes: {unassigned}"

    # Normalize group fracs relative to total of all composition columns
    total = arr[:, iprot] + arr[:, ineut] + group_fracs.sum(axis=0)
    group_fracs /= total
    neut_norm = arr[:, ineut] / total
    prot_norm = arr[:, iprot] / total

    with open(pathout, 'w') as out:
        out.write(f"{zones}\t{N_GROUPS + 2}\n")
        out.write(" ".join(f"{a:.1f}d0" for a in [1.0, 1.0] + GROUP_A) + " \n")
        out.write(" ".join(f"{z:.1f}d0" for z in [0.0, 1.0] + GROUP_Z) + " \n")

        # MESA stores surface→core; SNEC wants core→surface
        for i in range(zones, 0, -1):
            j = i - 1
            line  = "%15.6E" % (arr[j, imass] * msun)
            line += "%15.6E" % (arr[j, iradius] * rsun)
            line += "%15.6E" % neut_norm[j]
            line += "%15.6E" % prot_norm[j]
            for g in range(N_GROUPS):
                line += "%15.6E" % group_fracs[g, j]
            out.write(line + "\n")


if __name__ == "__main__":
    main()
