import logging
import glob
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import interpolate

from tardis.io.model import read_stella_model

from tardis_utils import write_tardis_csvy, write_tardis_config
import astropy.units as u

logger = logging.getLogger(__name__)

##########################################
###### Key options that affect the parser
INTERPOLATE_MASS_FRACTIONS = True  # boolean, if True then interpolate the mass fractions from MESA profile onto STELLA mass grid
SKIP_NONHOMOLOGOUS_MODELS = (
    True  # boolean, if True then skip non-homologous models and not save them
)
MAX_NONHOMOLOGOUS_SHELLS = 5  # int, if active if SKIP_NONHOMOLOGOUS_MODELS is True, the maximum number of non-homologous shells to skip
TAU_UPPER_LIMIT = 1e3  # False or float, filter out the shells that has tau larger than this value
SHRINK_SHELL_NUMBER = False  # False or int, if int then truncate the shell number to this value by skipping even number of shells in between
L_NUC_RATIO_UPPER_LIMIT = (
    0.8  # default 0.8, criteria to determine if the photosphere holds, means L_nuc/L_bol <= 0.8
)
##########################################

STELLA_model_folder = sys.argv[1]
parse_stella_models_to_tardis_configs(STELLA_model_folder, "TARDIS_example_configs")


def parse_stella_models_to_tardis_configs(
    stella_folder_path,
    tardis_example_config_folder_path,
    tardis_config_output_folder_path=None,
    interpolate_mass_fractions=INTERPOLATE_MASS_FRACTIONS,
    skip_nonhomologous_models=SKIP_NONHOMOLOGOUS_MODELS,
    max_nonhomologous_shells=MAX_NONHOMOLOGOUS_SHELLS,
    tau_upper_limit=TAU_UPPER_LIMIT,
    shrink_shell_number=SHRINK_SHELL_NUMBER,
    l_nuc_ratio_upper_limit=L_NUC_RATIO_UPPER_LIMIT,
):
    """
    Purpose:
    Parse the Stella model outputs (time series) to TARDIS config files with csvy format
    -------------------
    stella_folder_path: str, folder path to the stella model, need the res folder and other essential files
    tardis_example_config_folder_path: str, folder path to the tardis example config files
    tardis_config_output_folder_path: default None, will create a tardis_configs folder in stella model, path to save the tardis config
    interpolate_mass_fractions: boolean, if True then interpolate the mass fractions from MESA profile onto STELLA mass grid
    skip_nonhomologous_models: boolean, if True then skip non-homologous models and not save them
    max_nonhomologous_shells: int, if active if skip_nonhomologous_models is True, the maximum number of non-homologous shells to skip
    tau_upper_limit: False or float, filter out the shells that has tau larger than this value
    shrink_shell_number: False or int, if int then truncate the shell number to this value by skipping even number of shells in between
    l_nuc_ratio_upper_limit: float, criteria to determine if the photosphere holds, a fraction between 0 - 1, L_nuc/L_bol <= 0.8 is the default
    -------------------

    Returns
    -------------------
    Convert the stella model to a tardis config with csvy format, saving at
    """

    ###  check all the required files exsits in the folder
    # check the stella model files
    stella_output_files = sorted(glob.glob(f"{stella_folder_path}/res/mesa.day*"))
    if len(stella_output_files) == 0:
        logger.error(f"No stella model files found in the folder: {stella_folder_path}/res")
        return None

    # check the stella mesa.lbol_lnuc.txt file
    L_bol_file = f"{stella_folder_path}/res/mesa.lbol_lnuc.txt"
    if not Path(L_bol_file).exists():
        logger.error(f"No mesa.lbol_lnuc.txt file found in the folder: {stella_folder_path}/res")
        return None

    # check if the MESA profile exisits
    if interpolate_mass_fractions:
        mesa_profile_file = f"{stella_folder_path}/profile1.data"
        if not Path(mesa_profile_file).exists():
            logger.error(
                f"No mesa profile1.data file found in the folder: {stella_folder_path}/res"
            )
            return None

    # check if the tardis examples files exsits
    tardis_sample_config_path = (
        f"{tardis_example_config_folder_path}/tardis_example_config_SESN.yml"
    )
    tardis_sample_csvy_path = f"{tardis_example_config_folder_path}/tardis_example_csvy.csvy"
    if not Path(tardis_sample_config_path).exists():
        logger.error(
            f"No tardis example config file found in the folder: {tardis_example_config_folder_path}"
        )
        return None
    if not Path(tardis_sample_csvy_path).exists():
        logger.error(
            f"No tardis example csvy file found in the folder: {tardis_example_config_folder_path}"
        )
        return None

    # make the output folder if it doesn't exsit yet
    if tardis_config_output_folder_path is None:
        tardis_config_output_folder_path = f"{stella_folder_path}/tardis_configs"
    Path(tardis_config_output_folder_path).mkdir(parents=True, exist_ok=True)

    #######
    # extract the maximum day within photospheric assumption using l_nuc_ratio_upper_limit
    df_bol = pd.read_csv(
        L_bol_file,
        sep=r"\s+",
        skiprows=1,
        header=None,
        names=["time", "logL_bol", "logL_nuc"],
    )
    df_bol["L_nuc_ratio"] = 10 ** (df_bol["logL_nuc"] - df_bol["logL_bol"])
    max_photopheric_day = np.where(df_bol["L_nuc_ratio"].values <= l_nuc_ratio_upper_limit)[0][-1]

    # extract the days that have available stella profiles
    days_stella_profiles = np.array(
        [float(file.split("/")[-1].split("_")[0].split("day")[1]) for file in stella_output_files]
    )
    days_stella_profiles_photospheric = days_stella_profiles[
        days_stella_profiles <= max_photopheric_day
    ]

    ###### Read the stella model and convert to tardis config
    # loop through each profile check if they are monotonically increasing
    for i, day in enumerate(days_stella_profiles_photospheric):
        # read the stella model
        stella_model = read_stella_model(stella_output_files[i])
        df_stella_data = stella_model.data

        # filter out the optical thick shells
        if tau_upper_limit is not False:
            df_stella_data = df_stella_data[df_stella_data["tau"] <= tau_upper_limit].reset_index(
                drop=True
            )

        # check if the model is homologous
        if skip_nonhomologous_models is not False:
            non_homologous_shell = np.where(np.diff(df_stella_data["cell_center_v"]) < 0)[0]
            if non_homologous_shell.shape[0] > max_nonhomologous_shells:
                continue
            else:
                # filter out the non homologous shells
                df_stella_data = df_stella_data.drop(non_homologous_shell, axis=0).reset_index(
                    drop=True
                )

        # check if the user want to shrink the shell number
        if shrink_shell_number is not False:
            skip_shells = int(df_stella_data.shape[0] / shrink_shell_number)
            df_stella_data = df_stella_data.iloc[::skip_shells].reset_index(drop=True)

        # filter out the columns that are not needed for TARDIS
        matter_columns = ["cell_center_v", "avg_density", "radiation_temperature"]
        composition_columns_stella = [
            col for col in df_stella_data.columns if col[0].isalpha() and col[-1].isdigit()
        ]

        # interpolate the mass fractions based on MESA profile instead of using STELLA composition
        if interpolate_mass_fractions:
            df_stella_for_tardis = df_stella_data[matter_columns].copy()

            # read the mesa profile composition
            df_profile = pd.read_csv(mesa_profile_file, sep=r"\s+", skiprows=5)

            # get the isotopes from the composition columns
            composition_columns_profile = [
                col
                for col in df_profile.columns
                if col[0].isalpha()
                and col[-1].isdigit()
                and "_" not in col
                and col not in ["gamma1", "pnhe4"]
            ]

            # interpolate the mass fractions from mesa profile onto stella mass grid
            for isotope in composition_columns_profile:
                # first get a function that represent the MESA profile mass fraction
                mesa_mass_grid = (
                    df_profile["mass"].astype(np.float64).values[::-1]
                )  # MESA going inwards
                mesa_mass_fraction = df_profile[isotope].astype(np.float64).values[::-1]
                f_itp = interpolate.interp1d(
                    mesa_mass_grid,
                    mesa_mass_fraction,
                    kind="linear",
                    bounds_error=False,
                    fill_value=(mesa_mass_grid[0], mesa_mass_grid[-1]),
                )

                # interpolate onto the stella mass grid
                stella_mass_grid = df_stella_data["cell_center_m"].astype(
                    np.float64
                ).values * u.g.to(u.Msun)  # Stella has unit of g but MESA has units of Msun
                itp_mass_fraction = f_itp(stella_mass_grid)
                df_stella_for_tardis.loc[:, isotope[0].capitalize() + isotope[1:]] = (
                    itp_mass_fraction
                )

            # Stella seems to have na23 which is not in the appro21 net in MESA
            stella_unique_isotopes = [
                isotope
                for isotope in composition_columns_stella
                if isotope not in composition_columns_profile
            ]
            for isotope in stella_unique_isotopes:
                df_stella_for_tardis.loc[:, isotope[0].capitalize() + isotope[1:]] = df_stella_data[
                    isotope
                ]
        else:
            df_stella_for_tardis = df_stella_data[matter_columns + composition_columns_stella]
            df_stella_for_tardis = df_stella_for_tardis.rename(
                columns={col: col[0].capitalize() + col[1:] for col in composition_columns_stella}
            )

        # update the column names
        df_stella_for_tardis = df_stella_for_tardis.rename(
            columns={
                "cell_center_v": "velocity",
                "avg_density": "density",
                "radiation_temperature": "t_rad",
            }
        )

        # get the bolometric luminosity at the chosen day
        L_bol_at_chosen_day = (
            10 ** df_bol["logL_bol"].values[df_bol["time"].sub(day).abs().idxmin()]
        )

        # extract the time of the data relative to SBO
        day_since_SBO = df_bol["time"].min() + stella_model.metadata["t_max"].value

        # roughly estimate the photosphere index using tau = 2/3 for initial T_inner
        ph_idx = df_stella_data.index[df_stella_data["tau"].sub(1).abs().idxmin()]
        T_inner_guess = df_stella_data.loc[ph_idx, "radiation_temperature"]

        # write the tardis config file
        modify_parameters = {
            "supernova": {
                "luminosity_requested": f"{L_bol_at_chosen_day} erg/s",
                "time_explosion": f"{day_since_SBO:.4f} day",
            },
            "plasma": {"initial_t_inner": f"{T_inner_guess} K"},
        }
        new_config_path = f"{tardis_config_output_folder_path}/Day_{day}_mesa_stella_tardis.yml"
        write_tardis_config(
            tardis_sample_config_path,
            modify_parameters,
            new_config_path,
            csvy_model_path=f"Day_{day}_mesa_stella_model.csvy",
        )

        # write the tardis csvy file
        new_csvy_path = f"{tardis_config_output_folder_path}/Day_{day}_mesa_stella_model.csvy"
        modify_csvy_headers = {
            "name": "mesa_stella_model.csvy",
            "model_density_time_0": f"{day_since_SBO:.4f} day",
            "model_isotope_time_0": "0.0 s",
            "description": "mesa stella model converted to csvy format for tardis simulation",
            "v_inner_boundary": f"{df_stella_for_tardis['velocity'].min():.5e} cm/s",
        }
        write_tardis_csvy(
            tardis_sample_csvy_path, modify_csvy_headers, df_stella_for_tardis, new_csvy_path
        )
