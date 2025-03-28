#!/usr/bin/env python
import sys
import numpy as np

# hardcode constants as defined in MESA
global mu_sun, G, msun, rsun
mu_sun = 1.3271244e26
G = 6.67430e-8
msun = mu_sun / G
rsun = 6.957e10

# read profile*.data
infilename = sys.argv[1]
outfilename = sys.argv[2]
modfile = open(infilename, "r")
profile = modfile.readlines()
modfile.close()

# first three lines are going to be header
header = [[]] * 3
header[0] = profile[0].split()
header[1] = profile[1].split()
header[2] = profile[2].split()

# profile[3] is a blank line
# stuff starts at profile[4]
# get column names
nzones = len(profile[6:])
colnames = profile[5].split()

izone = colnames.index("zone")
irho = colnames.index("logRho")
itemp = colnames.index("logT")
ivel = colnames.index("velocity")
iye = colnames.index("ye")
try:
    iomega = colnames.index("omega")
except:
    pass
irad = colnames.index("radius")
imass = colnames.index("mass")

# data starts at profile[6]
arr = np.zeros((8, nzones))
for i in range(nzones):
    sline = profile[nzones + 5 - i].split()
    # mass
    # print i, sline[imass]
    arr[0, i] = float(sline[imass]) * msun
    # radius
    arr[1, i] = float(sline[irad]) * rsun
    # temp
    arr[2, i] = 10.0 ** float(sline[itemp])
    # rho
    arr[3, i] = 10.0 ** float(sline[irho])
    # vel
    arr[4, i] = float(sline[ivel])
    # Y_e
    arr[5, i] = float(sline[iye])
    # omega
    try:
        arr[6, i] = float(sline[iomega])
    except:
        arr[6, i] = 0.0

# radius of cell i is at its outer edge,
# so are omega, velocity, and mass
# density, temperature, Y_e are at cell centers

# GR1D assumes all quantities except r at the
# cell center, so we must move velocity
# and omega to cell center

# set up radius at cell centers
crad = np.zeros(nzones)
crad[0] = arr[1, 0] / 2.0
for i in range(1, nzones):
    crad[i] = (arr[1, i - 1] + arr[1, i]) / 2.0

# inner boundary: extrapolate
# velocity and omega
arr[4, 0] = (arr[4, 1] - arr[4, 0]) / (arr[1, 1] - arr[1, 0]) * (crad[0] - arr[1, 0]) + arr[4, 0]
arr[6, 0] = (arr[6, 1] - arr[6, 0]) / (arr[1, 1] - arr[1, 0]) * (crad[0] - arr[1, 0]) + arr[6, 0]
for i in range(1, nzones):
    arr[4, i] = (arr[4, i - 1] + arr[4, i]) / 2.0
    arr[6, i] = (arr[6, i - 1] + arr[6, i]) / 2.0

outfile = open(outfilename, "w")
buff = "%d\n" % nzones
outfile.writelines(buff)

for i in range(nzones):
    outline = "%6d " % (i + 1)
    for j in range(7):
        outline += " %18.9E " % (arr[j, i])
    outline += "\n"
    outfile.writelines(outline)


outfile.close()
print("done")
