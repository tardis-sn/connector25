import glob
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import interpolate

from tardis_utils import write_tardis_csvy, write_tardis_config
from tardis.util.base import (
    atomic_number2element_symbol,
    is_valid_nuclide_or_elem,
)

logger = logging.getLogger(__name__)

XG_FILE_NAMES = ["vel", "rho", "temp", "tau"]
DAT_FILE_NAMES = ["lum_observed", "T_eff", "vel_photo", "lum_photo", "index_photo"]


def xg_to_dict(fname):
    """
    parse the .xg file from SNEC output into a dictionary.
    credit: Brandon Barker and Chelsea Harris

    Parameters
    ----------
    fname : str

    Returns
    -------
    snec_xg_data : dictionary
    """
    snec_xg_data = {}
    with open(fname) as rf:
        for line in rf:
            cols = line.split()
            # Beginning of time data - make key for this time
            if "Time" in line:
                time = float(cols[-1])
                snec_xg_data[time] = []
            # In time data -- build x,y arrays
            elif len(cols) == 2:
                snec_xg_data[time].append(np.fromstring(line, sep=" "))
            # End of time data (blank line) -- make list into array
            else:
                snec_xg_data[time] = np.array(snec_xg_data[time])

    return snec_xg_data


def parse_snec_to_tardis(
    snec_folder_path,
    tardis_example_config_folder_path,
    tardis_config_output_path=None,
    time_in_days=None,
    tau_upper_limit=1e3,
    tau_lower_limit=1e-10,
    comps_profile_file_path=None,
    comp_use_boxcared=False,  # TODO: if true, need to rerun snec to have full isotopes
    snec_data_folder_path_comp_boxcar=None,
    use_vel_diff=False,
    num_keep_shells=45,
):
    """
    Purpose:
    ---------
    Parse the SNEC output to the TARDIS input config and csvy files.

    TARDIS Configuration Mapping from SNEC:
        supernova:
            luminosity_requested -- from lum_observed.dat at selected time (first col is time in seconds, second col is luminosity in erg/s)
            time_explosion -- the time key from any xg file (which has time in seconds as keys, the xg files are profile snapshots)

        plasma:
            initial_t_inner -- from T_eff.dat (Effective T of the photosphere) at selected time

        csvy_model: see description in subsection below

    CSVY model information:
        model_density_time_0 -- same as time_explosion
        model_isotope_time_0 -- the initial time of the composition profile used in SNEC
        v_inner_boundary -- vel_photo.dat at selected time (first col is time in seconds, second col is photospheric velocity in cm/s)
        velocity -- vel.xg (second column at selected time, first col is mass grid in unit of g)
        density  -- rho.xg (second column at selected time, first col is mass grid in unit of g)
        t_rad (optional)  -- temp.xg (second column at selected time, first col is mass grid in unit of g)
        composition -- profile/xx_comps.snec (need to match mass coordinate with velocity) at selected time

    ----------

    Parameters
    ----------
        snec_folder_path: str
            The path to the folder that contains the SNEC Data and profiles
        tardis_example_config_folder_path: str
            The path to the folder that contains the TARDIS example config files
        tardis_config_output_path: None or str
            If str -- The path to the folder that contains the input from TARDIS
            If None -- The output tardis files will not be produced
        time_in_days: list/array or None, default None
            If given, only the selected time (unit in days) will be used to generate the TARDIS input files, if None will use all available time steps
        tau_upper_limit: if not False, then setting a float value to the upper limit of the tau value based on the snec model profile to filter out the shells
        tau_lower_limit: if not False, then setting a float value to the lower limit of the tau value based on the snec model profile to filter out the shells
        comps_profile_file_path: str
            The path to the composition profile file from SNEC if specified
        comp_use_boxcared: bool, default True
            Use the boxcar smoothed composition profile from SNEC instead of the input file
        snec_data_folder_path_comp_boxcar: str or None, default None
            The path to the folder that contains the boxcar smoothed composition profile, if None will use the snec_folder_path
        use_vel_diff: bool
            If True, use the velocity adjacent difference to check homologous expansion,
            otherwise assume homologous for time step except the first one
        num_keep_shells: int
            The ROUGH number of shells to keep in the TARDIS model
    ----------
    output: (saved in tardis_config_output_path)
        time_series_log.csv
        tardis_config: yml files at each time
        tardis_model: csvy at each time

    """
    # get the folder name
    snec_folder_name = str(snec_folder_path).split("/")[-1]
    snec_data_folder_path = snec_folder_path / "output"

    # check if the tardis examples files exsits
    tardis_sample_config_path = (
        tardis_example_config_folder_path / "tardis_template_config_SESN.yml"
    )
    tardis_sample_csvy_path = tardis_example_config_folder_path / "tardis_example_csvy.csvy"

    # read the snec data
    if comps_profile_file_path is None:
        try:
            comp_files = glob.glob(f"{snec_folder_path}/input/*.iso.dat")
            snec_comps_profile_file_path = comp_files[0]
            if len(comp_files) > 1:
                Warning("More than one composition profile file found, please specify one.")
        except:
            ValueError("No composition profile file found.")
    else:
        snec_comps_profile_file_path = comps_profile_file_path

    dict_SNEC_output = snec_data_to_dict(snec_data_folder_path)

    # read the composition profile
    if comp_use_boxcared == True:
        if snec_data_folder_path_comp_boxcar is None:
            snec_data_folder_path_comp_boxcar = snec_data_folder_path
            print("Using the same folder for boxcar smoothed composition profile.")
        df_snec_comps = snec_boxcar_comps_profile_to_dataframe(
            snec_comps_profile_file_path,
            snec_data_folder_path_comp_boxcar,
            dict_SNEC_output["mass"],
        )
    else:
        df_comps = snec_comps_profile_to_dataframe(snec_comps_profile_file_path)
        # interpolate the composition profile to the mass grid of the SNEC output
        df_snec_comps = interpolate_composition_profile(df_comps, dict_SNEC_output)

    # generate the time mask for the selected time steps
    selected_time_mask = generate_time_mask(
        dict_SNEC_output["time"],
        dict_SNEC_output["vel"],
        dict_SNEC_output["vel_photo_profile"],
        use_vel_diff=use_vel_diff,
    )

    # filter the data to the selected time steps
    for param in ["time"] + XG_FILE_NAMES + [item + "_itp" for item in DAT_FILE_NAMES]:
        dict_SNEC_output[param] = dict_SNEC_output[param][selected_time_mask]

    if tardis_config_output_path is not None:
        # create the output folder if not exsit
        if not os.path.exists(tardis_config_output_path):
            os.makedirs(tardis_config_output_path)

        if time_in_days is None:
            # write tardis config and csvy file for each time step
            for time_index, _ in enumerate(dict_SNEC_output["time"]):
                new_csvy_path = (
                    f"{tardis_config_output_path}/{snec_folder_name}_tardis_csvy_{time_index}.csvy"
                )
                new_config_path = (
                    f"{tardis_config_output_path}/{snec_folder_name}_tardis_config_{time_index}.yml"
                )
                save_tardis_config_and_csvy(
                    dict_SNEC_output,
                    time_index,
                    df_snec_comps,
                    tardis_sample_csvy_path,
                    tardis_sample_config_path,
                    new_csvy_path,
                    new_config_path,
                    tau_upper_limit=tau_upper_limit,
                    tau_lower_limit=tau_lower_limit,
                    num_keep_shells=num_keep_shells,
                )
        else:
            # write tardis config and csvy file for the selected time steps only (unit in days)
            for time in time_in_days:
                time = float(time)
                time_index = np.argmin(np.abs(dict_SNEC_output["time"] - time * 24 * 3600))
                if np.abs(dict_SNEC_output["time"][time_index] / 24 / 3600 - time) > 1:
                    Warning(f"Time {time} day is not found in the SNEC output within +/-1d range.")
                    continue
                new_csvy_path = (
                    f"{tardis_config_output_path}/{snec_folder_name}_tardis_csvy_{time}_day.csvy"
                )
                new_config_path = (
                    f"{tardis_config_output_path}/{snec_folder_name}_tardis_config_{time}_day.yml"
                )
                save_tardis_config_and_csvy(
                    dict_SNEC_output,
                    time_index,
                    df_snec_comps,
                    tardis_sample_csvy_path,
                    tardis_sample_config_path,
                    new_csvy_path,
                    new_config_path,
                    tau_upper_limit=tau_upper_limit,
                    tau_lower_limit=tau_lower_limit,
                    num_keep_shells=num_keep_shells,
                )
                print("Saved TARDIS config and csvy file for time:", time)

    return dict_SNEC_output, df_snec_comps


def save_tardis_config_and_csvy(
    dict_SNEC_output,
    time_index,
    df_snec_comps,
    tardis_sample_csvy_path,
    tardis_sample_config_path,
    new_csvy_path,
    new_config_path,
    tau_upper_limit=False,
    tau_lower_limit=False,
    num_keep_shells=60,
):
    # get the time in day
    time_in_day = dict_SNEC_output["time"][time_index] / (60 * 60 * 24)

    # attached velocity, density, t_rad, dilution_factor, and composition profile data into a df
    df_profiles_non_comp = pd.DataFrame(
        {
            "velocity": dict_SNEC_output["vel"][time_index],
            "density": dict_SNEC_output["rho"][time_index],
            "t_rad": dict_SNEC_output["temp"][time_index],
            "tau": dict_SNEC_output["tau"][time_index],
        },
    )
    df_profiles = df_profiles_non_comp.join(df_snec_comps.reset_index(drop=True))

    # limit the shells numebers -- computational cost - This need to filtered before anything else ensure the index is correct
    if tau_upper_limit is not False:
        df_profiles = df_profiles.loc[df_profiles["tau"] <= tau_upper_limit]
    if tau_lower_limit is not False:
        df_profiles = df_profiles.loc[df_profiles["tau"] >= tau_lower_limit]
    if "tau" in df_profiles.columns:
        df_profiles = df_profiles.drop(columns=["tau"])

    # filter out the zero density and velocity shells on the outer region
    df_profiles = df_profiles.loc[(df_profiles.density > 0) & (df_profiles.velocity > 0)]

    # # filter out the inner shells that's purposely replaced with pure He
    # df_profiles = df_profiles.loc[df_profiles.He4 < 1]

    # filter out the outer shells that has t_radiative too low for TARDIS -> but this cause trouble though so replace low T shells with 500K instead
    df_profiles = df_profiles.loc[df_profiles.t_rad > 0]

    # discard the shells that has velocity backwards (which technically is in non-homologous expansion, but we can cut those out if there are only a few of them)
    while (df_profiles["velocity"].diff() <= 0).any():
        df_profiles = delete_non_increasing_neighbour(df_profiles, "velocity")

    # discard the element columns that has all zero values
    df_profiles = df_profiles.loc[:, (df_profiles != 0).any(axis=0)]

    if num_keep_shells is not None:
        shell_gap_length = max([int(np.floor(df_profiles.shape[0] / num_keep_shells)), 1])
        df_csv = df_profiles.groupby(np.arange(len(df_profiles)) // shell_gap_length).mean()
    else:
        df_csv = df_profiles

    # write the tardis csvy file
    modify_csvy_headers = {
        "name": new_csvy_path.split("/")[-1],
        "model_density_time_0": f"{time_in_day:.3f} day",
        "model_isotope_time_0": "0.0 s",
        # "v_inner_boundary": f"{dict_SNEC_output['vel_photo_itp'][time_index]:.6e} cm/s", # as of April2025, the v_inner workflow csvy model HAVE to the use the first shell to start to avoid dim error
        "v_inner_boundary": f"{df_profiles['velocity'].min():.6e} cm/s",
    }
    write_tardis_csvy(tardis_sample_csvy_path, modify_csvy_headers, df_csv, new_csvy_path)

    # write the tardis config file
    modify_parameters = {
        "supernova": {
            "luminosity_requested": f"{dict_SNEC_output['lum_observed_itp'][time_index]} erg/s",
            "time_explosion": f"{time_in_day:.3f} day",
        },
        "plasma": {"initial_t_inner": f"{dict_SNEC_output['T_eff_itp'][time_index]} K"},
    }
    write_tardis_config(
        tardis_sample_config_path,
        modify_parameters,
        new_config_path,
        csvy_model_path=new_csvy_path.split("/")[-1],
    )


def delete_non_increasing_neighbour(df, col_name):
    subdf = df[(df[col_name].diff() > 0)]
    if df.iloc[0][col_name] < subdf.iloc[0][col_name]:
        subdf = pd.concat([df.iloc[:1], subdf])
    return subdf


def generate_time_mask(time, velocity_arrays, vel_photo_profile, use_vel_diff=False):
    """
    Purpose:
    ---------
    Generate the time mask for the selected time steps.

    ----------

    Parameters
    ----------
        time: numpy array
            The time array
        velocity_arrays: numpy array
            The velocity 2D arrays
        vel_photo_profile: dictionary
            The dictionary that contains the photospheric velocity profile (time, vel_photo)
        use_vel_diff: bool
            If True, use the velocity adjacent difference to check homologous expansion,
            otherwise assume homologous for time step except the first one
    ----------
    output:
        selected_time_mask: numpy array
            The time mask for the selected time steps
    """
    ## select the time window that the ejecta is in homologous expansion and in photospheric phase
    if use_vel_diff == True:
        # calculate the velocity adjacent difference for use of checking homologous expansion
        diff_vel = np.diff(velocity_arrays, axis=1)
        homologous_mask = np.all(diff_vel >= 0, axis=1)
    else:
        homologous_mask = time > 0

    # taking the time when photospheric velocity gets to 0 as the time upper limit
    photospheric_time_limit = vel_photo_profile["time"][vel_photo_profile["vel_photo"] > 0][-1]
    photospheric_mask = time <= photospheric_time_limit

    selected_time_mask = np.where(homologous_mask & photospheric_mask)

    if selected_time_mask[0].size == 0:
        ValueError("No time step selected for homologous + photospheric phase.")

    return selected_time_mask


def snec_data_to_dict(snec_data_folder_path):
    """
    Purpose:
    ---------
    Parse the output from SNEC (a time series) to a single dictionary.


    ----------

    Parameters
    ----------
        snec_data_folder_path: str
            The path to the folder that contains the Data output from SNEC

    ----------
    output:
        A dictionary contains selected information of the SNEC output
    """
    dict_SNEC_output = {"vel": [], "rho": [], "temp": [], "tau": []}

    # read in the time steps that larger than 0 seconds (skipping the first time step)
    for i, param in enumerate(XG_FILE_NAMES):
        param_data = xg_to_dict(f"{snec_data_folder_path}/{param}.xg")
        if param == "vel":
            dict_SNEC_output["time"] = np.array(list(param_data.keys()))[
                1:
            ]  # the first time step is time 0
            dict_SNEC_output["mass"] = param_data[0].T[0]
        else:
            # check if the simulation time matches
            assert np.array_equal(dict_SNEC_output["time"], np.array(list(param_data.keys()))[1:])

        for time, data in param_data.items():
            if time > 0:
                # check if the mass grid matches
                assert np.array_equal(dict_SNEC_output["mass"], data.T[0])
                dict_SNEC_output[param].append(data.T[1])

        dict_SNEC_output[param] = np.array(dict_SNEC_output[param])

    for i, param in enumerate(DAT_FILE_NAMES):
        param_data = np.loadtxt(f"{snec_data_folder_path}/{param}.dat")
        # check if the simulation time matches
        dict_SNEC_output[param + "_profile"] = {
            "time": param_data.T[0],
            param: param_data.T[1],
        }
        # interpolate the data to the time grid
        f_itp = interpolate.interp1d(
            param_data.T[0],
            param_data.T[1],
            kind="linear",
            fill_value=(
                param_data.T[1][0],
                param_data.T[1][-1],
            ),
            bounds_error=False,
        )
        dict_SNEC_output[param + "_itp"] = f_itp(dict_SNEC_output["time"])
        if param == "index_photo":
            dict_SNEC_output[param + "_itp"] = dict_SNEC_output[param + "_itp"].astype(int)

    # # get the minimum above-zero photospheric velocity
    # if cut_inner_region:
    #     first_shell_above_photosphere_all_time = np.min(
    #         dict_SNEC_output["index_photo_profile"]["index_photo"][
    #             dict_SNEC_output["vel_photo_profile"]["vel_photo"] > 0
    #         ]
    #     )

    #     # cut the inner regions (defined as the shells below the lowest above-zero photospheric velocity)
    #     cut_index = max(
    #         [(int(first_shell_above_photosphere_all_time) - 1), 1]
    #     )  # avoid the first shell which has velocity 0
    #     dict_SNEC_output["mass"] = dict_SNEC_output["mass"][cut_index:]
    #     for i, param in enumerate(xg_params):
    #         dict_SNEC_output[param] = dict_SNEC_output[param][:, cut_index:]

    return dict_SNEC_output


def snec_comps_profile_to_dataframe(snec_comps_profile_file_path):
    """
    Purpose:
    ---------
    Parse the SNEC composition profile to a dataframe.
    """
    # Read the file and extract the second and third lines into arrays
    with open(snec_comps_profile_file_path) as file:
        lines = file.readlines()
        mass_numbers = np.array(lines[1].split()).astype(float).astype(int)
        atomic_numbers = np.array(lines[2].split()).astype(float).astype(int)

    # Convert the atomic numbers to element symbols
    element_symbols = [
        atomic_number2element_symbol(atomic_number) + str(mass_number)
        for atomic_number, mass_number in zip(atomic_numbers[1:], mass_numbers[1:])
    ]
    # check if the nuiclide is valid
    for element_symbol in element_symbols:
        if is_valid_nuclide_or_elem(element_symbol) == False:
            Warning(f"{element_symbol} is not valid nuiclide in tardis database.")

    # Read the remaining lines into a DataFrame
    df_abundance = pd.read_csv(
        snec_comps_profile_file_path,
        skiprows=3,
        sep=r"\s+",
        header=None,
    )
    df_abundance.columns = ["mass", "radius", "neutron"] + element_symbols

    return df_abundance


def snec_boxcar_comps_profile_to_dataframe(
    snec_comps_profile_file_path, snec_data_folder_path, snec_mass_grid
):
    """
    Purpose:
    ---------
    Parse the SNEC composition profile that are artificially smoothed by boxcar to a dataframe.
    """
    # Read the file and extract the second and third lines into arrays
    with open(snec_comps_profile_file_path) as file:
        lines = file.readlines()
        mass_numbers = np.array(lines[1].split()).astype(int)
        atomic_numbers = np.array(lines[2].split()).astype(int)

    # Convert the atomic numbers to element symbols
    element_symbols = [
        atomic_number2element_symbol(atomic_number) + str(mass_number)
        for atomic_number, mass_number in zip(atomic_numbers[1:], mass_numbers[1:])
    ]

    # get composition profile by isotopes
    df_abundance = pd.DataFrame(index=snec_mass_grid, columns=element_symbols)
    for i, element_symbol in enumerate(element_symbols):
        # check if the nuiclide is valid
        if is_valid_nuclide_or_elem(element_symbol) == False:
            Warning(f"{element_symbol} is not valid nuiclide in tardis database.")

        # Read the smoothed composition profile from the snec data folder
        iso_id = i + 2  # +2 due to 1.fortran start with 1, 2. first col is neutrinos
        file = f"{snec_data_folder_path}/iso_id_{iso_id}_init_frac.dat"
        df_init_frac = pd.read_csv(
            file,
            sep=r"\s+",
            header=None,
        )
        df_abundance[element_symbol] = df_init_frac[1].values

    return df_abundance


def interpolate_composition_profile(df_comps, dict_SNEC_output):
    """
    Purpose:
    ---------
    Interpolate the composition profile to the mass grid of the SNEC output.

    ----------

    Parameters
    ----------
        df_comps: dataframe
            The dataframe that contains the composition profile
        dict_SNEC_output: dictionary
            The dictionary that contains the SNEC output

    ----------
    output:
        A dataframe that contains the interpolated composition profile
    """
    # Get the mass grid from the data
    data_mass_grid = dict_SNEC_output["mass"]

    # Get the composition profile mass grid
    profile_mass_grid = df_comps["mass"].values

    # Get the composition profile data
    profile_data = df_comps.drop(columns=["mass", "radius", "neutron"])

    # Interpolate the composition profile to the data mass grid
    df_interpolated_comps = pd.DataFrame(index=data_mass_grid, columns=profile_data.columns)
    for column in profile_data.columns:
        f_itp = interpolate.interp1d(
            profile_mass_grid,
            profile_data[column],
            kind="linear",
            fill_value=(
                profile_data[column].values[0],
                profile_data[column].values[-1],
            ),
            bounds_error=False,
        )
        df_interpolated_comps[column] = f_itp(data_mass_grid)

    return df_interpolated_comps
