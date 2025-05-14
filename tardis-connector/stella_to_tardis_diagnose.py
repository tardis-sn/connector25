import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import astropy.units as u

from tardis.io.model import read_stella_model
from tardis.util.base import is_valid_nuclide_or_elem


def plot_profile_data(
    stella_model_path, x_col, y_col_s, tau_upper_limit=False, tardis_config_folder=None
):
    """Plot the profiles of the stella model and resulant tardis mapping

    Parameters
    ----------
    stella_model_path : str
        path to the folder where holds the stella models
    tardis_config_folder : str
        path to the folder where holds the tardis configuration files
    x_col : str, x axis in the plot
    y_col_s : list, y axis in the plot
    """
    if tardis_config_folder is not None:
        tardis_csvy_files = sorted(glob.glob(f"{tardis_config_folder}/*.csvy"))
        # stripp the days from the csvy files
        days_str = [
            csvy_file.split("/")[-1].split("_")[1] for csvy_file in tardis_csvy_files
        ]
    else:
        profile_files = sorted(glob.glob(f"{stella_model_path}/res/mesa.day*"))
        days_str = [
            file.split("/")[-1].split("_")[0].split("day")[1] for file in profile_files
        ]

    # generate a color set for each day
    color_set = cm.cool(np.linspace(0, 1, len(days_str)))

    fig, axes = plt.subplots(
        len(y_col_s), 1, figsize=(8, 2.5 * len(y_col_s)), sharex=True
    )
    fig.subplots_adjust(hspace=0)

    for i, day_str in enumerate(days_str):
        stella_model = read_stella_model(
            f"{stella_model_path}/res/mesa.day{day_str}_post_Lbol_max.data"
        )
        for j, y_col in enumerate(y_col_s):
            axes[j].plot(
                stella_model.data[x_col],
                stella_model.data[y_col],
                color=color_set[i],
                alpha=0.3,
                ls="--",
            )
        if tau_upper_limit is not False:
            if ("tau" in y_col_s) and i == 0:
                axes[y_col_s.index("tau")].axhline(
                    y=tau_upper_limit, color="k", linestyle="--", lw=0.5
                )
            df_model_stella = stella_model.data[
                stella_model.data["tau"] <= tau_upper_limit
            ]
            for j, y_col in enumerate(y_col_s):
                axes[j].plot(
                    df_model_stella[x_col],
                    df_model_stella[y_col],
                    color=color_set[i],
                    alpha=0.7,
                    label=f"Day {float(day_str)}",
                )
        else:
            df_model_stella = stella_model.data

    axes[0].legend(loc="upper right", fontsize=10, bbox_to_anchor=(1.2, 1))
    axes[-1].set_xlabel(x_col, fontsize=16)
    for i, ax in enumerate(axes):
        ax.set_yscale("log")
        ax.set_ylabel(y_col_s[i], fontsize=16)
        ax.tick_params(axis="both", which="major", labelsize=14)
        ax.grid(alpha=0.3)

    return fig


def plot_bolometric_LC(
    STELLA_model_folder,
    L_NUC_RATIO_UPPER_LIMIT=None,
    tardis_config_folder=None,
    ax=None,
):
    lbol_file = f"{STELLA_model_folder}/res/mesa.lbol_lnuc.txt"
    df_bol = pd.read_csv(
        lbol_file,
        sep=r"\s+",
        skiprows=1,
        header=None,
        names=["time", "logL_bol", "logL_nuc"],
    )
    df_bol_with_uBVri = pd.read_csv(
        f"{STELLA_model_folder}/res/mesa.lbol",
        sep=r"\s+",
    )
    if ax is None:
        fig, axes = plt.subplots(2, 1, figsize=(8, 12))
        ax = axes[0]

    if L_NUC_RATIO_UPPER_LIMIT is not None:
        df_bol["L_nuc_ratio"] = 10 ** (df_bol["logL_nuc"] - df_bol["logL_bol"])
        max_photospheric_day = df_bol.time.values[
            df_bol["L_nuc_ratio"].values <= L_NUC_RATIO_UPPER_LIMIT
        ][-1]
        df_bol = df_bol[df_bol["time"] <= max_photospheric_day + 10]
        ax.axvline(max_photospheric_day, color="k", ls="--", label="Photospheric limit")

    # mark the time have stella turned tardis configs
    if tardis_config_folder is not None:
        tardis_csvy_files = glob.glob(f"{tardis_config_folder}/*.csvy")
        days_tardis = [
            float(csvy_file.split("/")[-1].split("_")[1])
            for csvy_file in tardis_csvy_files
        ]
        for day in days_tardis:
            day_index = np.argmin(np.abs(df_bol["time"].values - day))
            ax.scatter(
                day, df_bol["logL_bol"].values[day_index], marker="o", color="tab:blue"
            )

    # plot the bolometric luminosity
    ax.plot(df_bol["time"], df_bol["logL_bol"], label="logL_bol")
    ax.plot(df_bol["time"], df_bol["logL_nuc"], label="logL_nuc")
    ax.plot(df_bol_with_uBVri["time"], df_bol_with_uBVri["L_ubvri"], label="logL_UBVRI")

    axes[1].plot(df_bol["time"], df_bol["L_nuc_ratio"], label="L_nuc/L_bol")
    axes[1].plot(
        df_bol_with_uBVri["time"],
        10 ** (df_bol_with_uBVri["L_ubvri"] - df_bol_with_uBVri["L_bol"]),
        label="L_UBVRI/L_bol",
    )
    axes[1].legend(fontsize=16)
    axes[1].grid(alpha=0.5)
    axes[1].set_ylim(-1, 1.2)

    ax.legend(fontsize=16)
    ax.grid(alpha=0.5)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_xlabel("Days since peak", fontsize=18)
    ax.set_ylabel("$log_{10}$ $L_{bol}$  (erg/s)", fontsize=18)
    ax.set_xlim(-2, df_bol["time"].values[-1])

    return ax
