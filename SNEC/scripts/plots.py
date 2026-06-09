import numpy as np
import sys
from SNEC_output_parser import SNEC_output_parser
import matplotlib.pyplot as plt


def plot_v_r(mass_file, velocity_file):
    data = SNEC_output_parser(velocity_file)
    times = np.array(list(data.keys())) # in sec

    # to map on radius
    data_r = SNEC_output_parser(mass_file)

    fig, ax = plt.subplots()
    for i, t in enumerate(times):
        if i < 917: continue
        # ax.plot(data[t][:, 0]/2e33, data[t][:, 1]/1e8)
        ax.plot(data_r[t][:, 0], data[t][:, 1]/1e8)
    # ax.set_xlabel(r"mass coordinate [$M_{\odot}$]")
    ax.set_ylabel(r"velocity [km/s]")
    ax.set_xlabel(r"r  [cm]")
    # ax.set_yscale('log')
    plt.show()


def plot_v_m(velocity_file):
    data = SNEC_output_parser(velocity_file)
    times = np.array(list(data.keys())) # in sec
    fig, ax = plt.subplots()
    for i, t in enumerate(times):
        if i < 917: continue
        ax.plot(data[t][:, 0]/2e33, data[t][:, 1]/1e8)
    ax.set_xlabel(r"mass coordinate [$M_{\odot}$]")
    ax.set_ylabel(r"velocity [km/s]")
    plt.show()


def plot_LC(obs_lum, ax=None, **kwargs):
    # print(folder)
    # obs_lum = folder+'Data/lum_observed.dat'
    src = np.genfromtxt(obs_lum)
    t = src[:, 0]
    L = src[:,1]
    if not ax:
        fig, ax = plt.subplots()
    # ax.scatter(t.to(u.d), np.log10(L.value), **kwargs)
    ax.scatter(t[L>0]/86400, np.log10(L[L>0]), **kwargs)
    plt.show()
    return L, t



if __name__ == "__main__":
    root = "../output/large_net/"
    # plot_v_r(root+"mass.xg", root+"vel.xg")
    # plot_v_m(root+"vel.xg")
    L, t = plot_LC(root+"lum_photo.dat")
