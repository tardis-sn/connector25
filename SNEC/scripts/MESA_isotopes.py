#!/usr/bin/env python

import sys
import numpy as np

path = sys.argv[1]  # first argument is the path of the MESA profile we want to convert
pathout = sys.argv[2] # second argument is the path where to save the output (must include the extension)

msun = 1.99e33
rsun = 6.96e10

modfile = open(path, "r")
profile = modfile.readlines()
modfile.close()
# first three lines are going to be header
header = [[]] * 3
header[0] = profile[0].split()
header[1] = profile[1].split()
header[2] = profile[2].split()

# get the number of zones (from the header)
izones = header[1].index("num_zones")
zones = int(
    header[2][izones]
)  # header[2] is a list, the second pair of brackets indicates which element of the list
# get the total mass
itotmass = header[1].index("star_mass")
totalMass = float(header[2][itotmass]) * msun

# get the columns names
columnNames = profile[5].split()
# get some more indexes
imixing = columnNames.index("mixing_type")
imass = columnNames.index("mass")  # mass is in M_\odot units
iradius = columnNames.index("radius")  # outer radius of the cell in R_\odot units
iprot = columnNames.index("h1")  # index for column containing protons
ineut = columnNames.index("neut")  # index for column containing free neutrons

arr = np.loadtxt(path, skiprows=6)

for i in range(0, len(arr[:, imixing]), 1):  # start reading the the data
    if arr[i, imixing] == -1:
        arr[i, imixing] = 0  # crystallized
    elif arr[i, imixing] == 0:
        arr[i, imixing] = 0  # no mixing
    elif arr[i, imixing] == 1:
        arr[i, imixing] = 0  # convection
    elif arr[i, imixing] == 2:
        arr[i, imixing] = -1  # overshooting
    elif arr[i, imixing] == 3:
        arr[i, imixing] = 0.5  # semiconv
    elif arr[i, imixing] == 4:
        arr[i, imixing] = 0  # thermohaline
    elif arr[i, imixing] == 5:
        arr[i, imixing] = 0  # rotation
    elif arr[i, imixing] == 6:
        arr[i, imixing] = 0  # RTI
    elif arr[i, imixing] == 7:
        arr[i, imixing] = 0  # minimum mixing
    elif arr[i, imixing] == 8:
        arr[i, imixing] = 0  # anonymous
    elif arr[i, imixing] == 9:
        arr[i, imixing] = 0  # decaying convection (TDC)
    elif arr[i, imixing] == 10:
        arr[i, imixing] = 0  # phase separation mixing
    else:
        arr[i, imixing] = 0

isolist = []
isolist.append(["he3", "he4", "h2", "h3", "li6", "li7", "li8", "be7", "be9", "b8"])
isolist.append(["c11", "c12", "c13", "c14", "n13", "n14", "n15", "b10", "b11"])
isolist.append(["o15", "o16", "o17", "o18", "f17", "f18", "f19"])
isolist.append(["ne20", "ne21", "ne22", "ne23", "na21", "na22", "na23", "na24"])
isolist.append(["mg23", "mg24", "mg25", "mg26", "mg27", "al25", "al26", "al27", "al28"])
isolist.append(
    [
        "si27",
        "si28",
        "si29",
        "si30",
        "si31",
        "si32",
        "p29",
        "p30",
        "p31",
        "p32",
        "p33",
        "p34",
    ]
)
isolist.append(
    [
        "s31",
        "s32",
        "s33",
        "s34",
        "s35",
        "s36",
        "s37",
        "cl33",
        "cl34",
        "cl35",
        "cl36",
        "cl37",
        "cl38",
    ]
)
isolist.append(
    [
        "ar36",
        "ar37",
        "ar38",
        "ar39",
        "ar40",
        "ar41",
        "k37",
        "k38",
        "k39",
        "k40",
        "k41",
        "k42",
    ]
)
isolist.append(
    [
        "ca40",
        "ca41",
        "ca42",
        "ca43",
        "ca44",
        "ca45",
        "ca46",
        "ca47",
        "ca48",
        "ca49",
        "sc41",
        "sc42",
        "sc43",
        "sc44",
        "sc45",
        "sc46",
        "sc47",
        "sc48",
        "sc49",
        "sc50",
    ]
)
isolist.append(
    [
        "ti44",
        "ti45",
        "ti46",
        "ti47",
        "ti48",
        "ti49",
        "ti50",
        "ti51",
        "v45",
        "v46",
        "v47",
        "v48",
        "v49",
        "v50",
        "v51",
        "v52",
    ]
)
isolist.append(
    [
        "cr48",
        "cr49",
        "cr50",
        "cr51",
        "cr52",
        "cr53",
        "cr54",
        "cr55",
        "cr56",
        "mn51",
        "mn52",
        "mn53",
        "mn54",
        "mn55",
        "mn56",
        "mn57",
    ]
)
isolist.append(
    [
        "fe52",
        "fe53",
        "fe54",
        "fe55",
        "fe56",
        "fe57",
        "fe58",
        "fe59",
        "fe60",
        "fe61",
        "co56",
        "co57",
        "co58",
        "co59",
        "co60",
        "co61",
        "co62",
    ]
)
isolist.append(
    [
        "ni56",
        "ni57",
        "ni58",
        "ni59",
        "ni60",
        "ni61",
        "ni62",
        "ni63",
        "ni64",
        "ni65",
        "cu57",
        "cu58",
        "cu59",
        "cu60",
        "cu61",
        "cu62",
        "cu63",
        "cu64",
        "cu65",
        "cu66",
        "zn60",
        "zn61",
        "zn62",
        "zn63",
        "zn64",
        "zn65",
        "zn66",
        "zn67",
        "zn68",
        "zn69",
        "ga62",
        "ga63",
        "ga64",
        "ga65",
        "ga66",
        "ga67",
        "ga68",
        "ga69",
        "ga70",
        "ge64",
        "ge65",
        "ge66",
        "ge67",
        "ge68",
        "ge69",
        "ge70",
        "ge71",
    ]
)
# added in the 11th group cr56, which is traced by MESA with approx26.net

iso_num = len(isolist)

Z = [2, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]
A = [4, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56]

massfracsShort = (
    []
)  # will be a list of the isotope mass fraction profiles once the isotopes are grouped

"""
massfracs_p = arr[:,columnNames.index('h1')]
ye = arr[:,columnNames.index('ye')]
norm = (ye - massfracs_p)*2.0
massfracs_n = 1.0 - massfracs_p - norm
massnorm = norm/(totalMass)
print massnorm
"""
for i in range(
    len(isolist)
):  # here we add mass fractions of the isotope groupings together
    group = isolist[i]
    temp = []
    groupMassfracs = []
    for isotope in group:
        if (
            isotope in columnNames
        ):  # skip the isotopes MESA doesn't trace (as far as MESA is concerned, they are just not there)
            ind = int(columnNames.index(isotope))
            groupMassfracs.append(arr[:, ind])  # at the end of the loop groupMassfracs
            # is a matrix: every line is an isotope, every column a zone
        groupArr = np.array(groupMassfracs)  # cast it into numbers
        TgroupArr = (
            groupArr.T
        )  # .T means transpose: we get a matrix in which each column is for an isotope and each row for a zone
    # we have a matrix for each group
    zonesum = []
    for j in range(0, zones):
        zonesum.append(sum(TgroupArr[j, :]))
        temp.append(
            zonesum[j]
        )  # no need to normalize with MESA (judging from the output values orders of magnitude)
    massfracsShort.append(temp)

# OUTPUT
outfile = open(pathout, "w")

outfile.write(
    str(zones) + "\t" + str(iso_num + 2) + "\n"
)  # this is the first line of the output file, with numbers of zones and isotopes

outfile.write(
    "1.0d0 1.0d0 4.0d0 12.0d0 16.0d0 20.0d0 24.0d0 28.0d0 32.0d0 36.0d0 40.0d0 44.0d0 48.0d0 52.0d0 56.0d0 \n"
)  # values of A for each isotope, starting with n, p, He-4, C-12...

outfile.write(
    "0.0d0 1.0d0 2.0d0 6.0d0 8.0d0 10.0d0 12.0d0 14.0d0 16.0d0 18.0d0 20.0d0 22.0d0 24.0d0 26.0d0 28.0d0 \n"
)  # values of Z for each isotope

# MESA stores the output in reverse order with respect to KEPLER
for i in range(zones, 0, -1):
    writeNewLine = "%15.6E" % (arr[i - 1, imass] * msun)  # mass in grams
    writeNewLine += "%15.6E" % (arr[i - 1, iradius] * rsun)  # radius in cm
    writeNewLine += "%15.6E" % (
        arr[i - 1, ineut]
    )  # abundance of neutrons as MESA gives it
    writeNewLine += "%15.6E" % (
        max(1e-40, arr[i - 1, iprot]) # floor H to small value for stripped stars
    )  # abundance of h1 (no protons from photodisintegration)

    for newIso in massfracsShort:
        writeNewLine += "%15.6E" % (
            newIso[i - 1]
        )  # sum of the abundances of each element in the group
        # elements that I don't have in the output are just skipped
    outfile.write(writeNewLine + "\n")


outfile.close()
