from scipy.interpolate import RegularGridInterpolator as rgi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
import h5py
import yt

# Load configuration
stream = open('config.yaml', 'r')
configs = yaml.safe_load(stream)
progenitor_directory = configs['progenitor_directory']
progenitor_suffix = configs['progenitor_suffix']
stir_profiles_directory = configs['stir_profiles_directory']
stir_profiles_suffix = configs['stir_profiles_suffix']
output_directory = configs['output_directory']
output_suffix = configs['output_suffix']
cell_edge_velocity = configs['cell_edge_velocity']
eos_file_path = configs['eos_file_path']
velocity_col_name = "v" if cell_edge_velocity else "u"

# Constants in CGS units
M_sun = 1.989E33 # Mass of the sun in grams
R_sun = 6.959E10 # Radius of the sun in centimeters
sigma_b = 5.669E-5 # Stefan-Boltzmann constant
G = 6.67430E-8 # Gravitational constant

default_plotted_profiles = ['enclosed_mass', 'r', 'density', 'temp', velocity_col_name, 'total_specific_energy', 'pressure']
log_plots = ["enclosed_mass", "r", "density", "pressure", "temp"]

    
def convert(model_name, stir_alpha, stir_portion = 0.8, plotted_profiles = ["DEFAULT"]):
    """
    Loads in a MESA progenitor and STIR profiles, making various modifications to make 
    them compatible, then creates a MESA-readable .mod file.

    Parameters:
        model_name (str) :   
            The name of the model (what comes before .mod or .data in the MESA progenitors).
        stir_alpha (float) : 
            The alpha value used for the STIR simulations.
        stir_portion (float) :
            What fraction of the STIR domain to include. A value of 0.8 will use progenitor data for the last 20% of the STIR domain.
        plotted_profiles (numpy array (str)) : 
            The names of each profile/variable you want to see plotted. Use ["COMPOSITION"] to plot only composition. Default of ["DEFAULT"] will plot enclosed mass, radius, density, temperature, velocity, total specific energy, and pressure.
            The profiles available for plotting are: enclosed_mass, density, temp, r, L, dq, v, mlt_vc, ener, pressure, and any nuclear network composition
    """

    prog = load_progenitor(model_name)
    stir = load_stir_profiles(f"{model_name}_a{stir_alpha}", prog["nuclear_network"])
    data = combine_data(stir, prog, stir_portion)
    write_mesa_model(data, prog, f"{model_name}_a{stir_alpha}")
    
    # Plot the desired profiles
    if "DEFAULT" in plotted_profiles: plotted_profiles = default_plotted_profiles
    if "COMPOSITION" in plotted_profiles: plotted_profiles = prog["nuclear_network"]
    for profile in plotted_profiles:
        plot_profile(data, profile)
    

def combine_data(stir, prog, stir_portion):
    """
    Combines data from the STIR domain with the progenitor data outside that domain.
    """

    # Determines the end of the STIR domain and start of the progenitor domain
    stir_domain = stir["r"].values[stir["r"].values <= np.max(stir["r"].values) * stir_portion]
    data = { "stir_domain_end": np.argmax(stir_domain) }
    prog_domain = len(prog['profiles'].loc[prog['profiles']['enclosed_mass'].values > np.max(stir['enclosed_mass'].values[:data["stir_domain_end"]])])

    # Combines STIR and progenitor data, simply placing the progenitor data at the end of the STIR domain
    data["profiles"] = pd.DataFrame(index = pd.RangeIndex(data["stir_domain_end"] + prog_domain), columns = stir.columns, dtype=float)
    for col in stir.columns:
        
        # If the column doesn't exist in progenitor data, notify us and skip it
        if not(col in prog['profiles']):
            print(f"Column {col} missing from progenitor data.")
            continue
        
        # Fill in data for each domain
        data["profiles"][col].values[:data["stir_domain_end"]] = stir[col].values[:data["stir_domain_end"]]
        data["profiles"][col].values[data["stir_domain_end"]:] = prog['profiles'][col].values[-prog_domain:]

    # Determine the total specific energy in each zone
    data["profiles"] = data["profiles"].assign(total_specific_energy = data["profiles"]['ener'].values + data["profiles"]['gpot'].values)

    # Set the PNS mass as the enclosed mass within which all cells have a negative total specific energy
    data["pns_masscut_index"] = np.min(np.where(data["profiles"]['total_specific_energy'].values >= 0)) - 1
    data["pns_masscut"] = data["profiles"]['enclosed_mass'].values[data["pns_masscut_index"]]
    data["pns_radius"] = data["profiles"]['r'].values[data["pns_masscut_index"]]

    # Find the mass of the star, and of the star outside the PNS
    data["total_mass"] = np.sum(data["profiles"]['density'] * data["profiles"]['cell_volume'])
    data["xmstar"] = data["total_mass"] - data["pns_masscut"] * M_sun
    data["profiles"]["dq"] = data["profiles"]['density'] * data["profiles"]['cell_volume'] / data["xmstar"]

    # Excise all zones within the PNS radius, since MESA does not want them
    data["pre_masscut_profiles"] = data["profiles"].copy()
    data["profiles"] = data["profiles"].drop(np.arange(data["pns_masscut_index"] + 1))
    
    # Calculate the remaining energy of the star without the PNS
    data["total_energy"] = data["profiles"]["total_specific_energy"].values * data["profiles"]['density'].values * data["profiles"]['cell_volume'].values

    # Prepares the data for MESA output by adding missing columns
    data["profiles"] = data["profiles"].assign(lnR = np.log(data["profiles"]['r'].values), lnd = np.log(data["profiles"]['density'].values), lnT = np.log(data["profiles"]['temp'].values))
    data["profiles"] = data["profiles"].assign(mlt_vc = np.zeros(data["profiles"].shape[0]))
    
    # MESA needs the surface luminosity which is the same as in progenitor, but all other values should just be a copy of that value
    data["profiles"] = data["profiles"].assign(L = np.ones(data["profiles"].shape[0]) * prog['profiles']["L"].values[-1])

    return data


def load_progenitor(model_name):
    '''Loads a MESA progenitor.'''

    prog = {}
    profile_path = f"{progenitor_directory}/{model_name}{progenitor_suffix}.data"
    model_path = f"{progenitor_directory}/{model_name}{progenitor_suffix}.mod"

    # Reads from the .mod file to get necessary header and footer information
    with open(model_path, "r") as file:

        # Copies the first few lines of the header to ensure bit flags and units are preserved
        lines = file.readlines()
        prog["header_start"] = lines[0] + lines[1] + lines[2]

        # Find the table header and from it grab the nuclear network
        for line in lines:
            if "                lnd" in line:
                prog["table_header"] = line
                break
        prog["nuclear_network"] = prog["table_header"].split()[7:]

        # Load each header and footer variable individually
        for i in np.concat((range(4, 20), range(-3, -1))):
            line_data = lines[i].split()
            value = line_data[-1]
            if "D+" in value or "D-" in value: prog[line_data[0]] = float(value.replace("D", "e")) 
            elif value.isdigit(): prog[line_data[0]] = int(value)
            else: prog[line_data[0]] = value

        # These are on the same line so they have to be handled separately
        prog["cumulative_error/total_energy"] = float(lines[20].split()[1].replace("D", "e"))
        prog["log_rel_run_E_err"] = float(lines[20].split()[3].replace("D", "e"))

    # Load profiles of each variable
    with open(profile_path, "r") as file:
        lines = file.readlines()

        # Find the columns that contain the necessary data
        needed_profiles = np.concat((['mass', 'logRho', 'temperature', 'radius_cm', 'luminosity', 
                                    'logdq', 'velocity', 'conv_vel', 'energy', 'pressure'], prog["nuclear_network"]))
        input_column_names = lines[5].split()
        column_indices = [input_column_names.index(col) for col in needed_profiles if col in input_column_names]

        # Create a 2D array containing the numerical data
        structured_data = [list(map(float, np.array(line.split())[column_indices])) for line in lines[6:]]

        # If any columns are missing from the progenitor, fill them with a very small number. 
        # This is a failsafe, but shouldn't be necessary since the composition we need is pulled from the .mod file.
        for i, col in enumerate(needed_profiles):
            if col not in input_column_names:
                print(f"Column {col} missing from progenitor data. Filling with very small number for all cells.")
                for j in range(len(structured_data)): 
                    structured_data[j].insert(i, 1e-99)

        # Write the numerical data into a pandas dataframe, renaming columns for compatibility and ease of use
        output_column_names = np.concat((['enclosed_mass', 'density', 'temp', 'r', 'L', 'dq', velocity_col_name, 
                                        'mlt_vc', 'ener', 'pressure'], prog["nuclear_network"]))
        prog["profiles"] = pd.DataFrame(structured_data, columns=output_column_names)
        
        # Invert order of zones since loaded MESA data will have first zone as outer radius
        prog["profiles"] = prog["profiles"].iloc[::-1].reset_index(drop=True)

        # Convert data to the correct scale for compatibility with STIR output
        prog["profiles"]["density"] = (10 ** prog["profiles"]['density'])
        prog["profiles"]["dq"] = 10 ** prog["profiles"]['dq']

        # Calculates the progenitor volume and gravitational potential
        prog_volume = (4/3) * np.pi * np.diff(np.concatenate(([0], prog["profiles"]['r'].values**3)))
        prog["profiles"] = prog["profiles"].assign(cell_volume = prog_volume) 
        prog["profiles"] = prog["profiles"].assign(gpot = -G * prog["profiles"]['enclosed_mass'].values / prog["profiles"]['r'].values)

    return prog


def load_stir_profiles(model_name, nuclear_network):
    '''Reads in data from a STIR checkpoint (or plot) file, making modifications and returning it as a dataframe.'''
    
    # Load the STIR checking/plot data into a dataframe
    profile_path = f"{stir_profiles_directory}/{model_name}{stir_profiles_suffix}"
    stir_ds = yt.load(profile_path)
    stir_data = stir_ds.all_data()
    needed_profiles = [("gas", "density"), ("flash", "temp"), ("gas", "r"), ("flash", "velx"), ("gas", "pressure"),# ("flash", "ye  "),
                ("flash", "cell_volume"), ("flash", "ener"), ("flash", "gpot")]
    stir = stir_data.to_dataframe(needed_profiles)
    
    # Try to pull composition from STIR output. If anything is missing, fill it with 1e-99
    for nuc in nuclear_network: 
        if ('flash', f'{nuc}') in stir_ds.field_list: 
            stir[nuc] = stir_data[nuc].v
        else: 
            print(f"{nuc} missing from the STIR profiles. Filling with 1e-99.")
            stir[nuc] = np.ones(stir.shape[0]) * 1e-99
    
    # Renormalize the composition mass fractions since they need to sum to exactly 1 with double precision
    for i in range(stir.shape[0]):
        mass_fraction_sum = stir.loc[i, nuclear_network].sum()
        stir.loc[i, nuclear_network] /= mass_fraction_sum

    # To better match MESA outputs, shift radius to cell edge and rename velocity
    stir['r'] = shift_to_cell_edge(stir['r'].values) 
    stir.rename(columns={"velx": velocity_col_name}, inplace=True)

    # If using cell edge velocity, shift the STIR velocities to the cell edge as well
    if cell_edge_velocity: 
        stir[velocity_col_name] = shift_to_cell_edge(stir[velocity_col_name].values)
    
    # Calculate the enclosed mass for every zone
    enclosed_mass = np.cumsum(stir['cell_volume'].values * stir['density'].values) / M_sun
    stir = stir.assign(enclosed_mass = enclosed_mass)

    # Load the equation of state used by STIR
    EOS = h5py.File(eos_file_path,'r')
    mif_logenergy = rgi((EOS['ye'], EOS['logtemp'], EOS['logrho']), EOS['logenergy'][:,:,:], bounds_error=False)

    # Use the STIR EOS to calculate the total specific energy
    lye = stir_data['ye  ']
    llogtemp = np.log10(stir_data['temp'] * 8.61733326e-11)
    llogrho = np.log10(stir_data['dens'])
    energy = 10.0 ** mif_logenergy(np.array([lye, llogtemp, llogrho]).T)
    llogtemp = llogtemp * 0.0 - 2.0
    energy0 = 10.0 ** mif_logenergy(np.array([lye, llogtemp, llogrho]).T)
    stir["ener"] = (0.5 * stir_data['velx'] ** 2 + (energy - energy0) * yt.units.erg / yt.units.g).v

    return stir


def write_mesa_model(data, prog, model_name):
    '''Writes the star's data into MESA input files.''' 
    
    # Easy way to convert values into the correct format for MESA model files
    def format_float(num):
        return f"{num:.16e}".replace('e', 'D')
    
    def format_int(num):
        return ' ' * (25 - int(np.floor(np.log10(num)))) + str(num)
    
    # Header Info
    avg_core_density = data["pns_masscut"] * M_sun / (4/3 * np.pi * data["pns_radius"]**3) 
    file_header = prog["header_start"] + f"""
                  version_number   {prog["version_number"]}
                          M/Msun      {format_float(data["total_mass"] / M_sun)}
                    model_number      {format_int(prog["model_number"])}
                        star_age      {format_float(prog["star_age"])}
                       initial_z      {format_float(prog["initial_z"])}
                        n_shells      {format_int(data["profiles"].shape[0])}
                        net_name   {prog["net_name"]}
                         species      {format_int(prog["species"])}
                          xmstar      {format_float(data["xmstar"])}  ! above core (g).  core mass: Msun, grams:      {format_float(data["pns_masscut"])}    {format_float(data["pns_masscut"] * M_sun)}
                        R_center      {format_float(data["pns_radius"])}  ! radius of core (cm).  R/Rsun, avg core density (g/cm^3):      {format_float(data["pns_radius"] / R_sun)}    {format_float(avg_core_density)}
                            Teff      {format_float(prog["Teff"])}
                  power_nuc_burn      {format_float(prog["power_nuc_burn"])}
                    power_h_burn      {format_float(prog["power_h_burn"])}
                   power_he_burn      {format_float(prog["power_he_burn"])}
                    power_z_burn      {format_float(prog["power_z_burn"])}
                     power_photo      {format_float(prog["power_photo"])}
                    total_energy      {format_float(np.sum(data["total_energy"]))}
         cumulative_energy_error      {format_float(prog["cumulative_energy_error"])}
   cumulative_error/total_energy      {format_float(prog["cumulative_error/total_energy"])}  log_rel_run_E_err      {format_float(prog["log_rel_run_E_err"])}
                     num_retries                               0

""" + prog["table_header"]

    # Data that will be written in table form in the stir_output.mod file
    # Also reverses the order of rows in the table so that the first cell is the outer radius and last cell is the center
    mesa_columns = np.concat((['lnd', 'lnT', 'lnR', 'L', 'dq', velocity_col_name, 'mlt_vc'], prog["nuclear_network"]))
    mesa_input = data["profiles"][mesa_columns].iloc[::-1].reset_index(drop=True)

    # Add one line for each cell, consisting of all it's properties
    new_lines = []
    for line_index in range(data["profiles"].shape[0]):

        # Writes the cell/line index 
        spaces = 4 - int(np.floor(np.log10(line_index + 1)))
        new_line = ' ' * spaces + str(line_index + 1)

        # Writes each of the properties
        for column_name in mesa_input.columns:
            spaces = 5 if mesa_input.at[line_index, column_name] >= 0 else 4
            new_line += spaces * ' ' + format_float(mesa_input.at[line_index, column_name])

        new_lines.append(new_line)

    # Footer containing info about the previous model
    file_footer = f"""
    
        previous model

               previous n_shells      {format_int(prog['profiles'].shape[0])}
           previous mass (grams)      {format_float(prog["M/Msun"] * M_sun)}
              timestep (seconds)      {format_float(prog["timestep"])} 
               dt_next (seconds)      {format_float(prog["dt_next"])}

"""

    # Write all of the above to the stir_output.mod file
    output_path = f"{output_directory}/{model_name}{output_suffix}.mod"
    with open(f'{output_path}', 'w') as file:
        file.writelines(file_header + '\n'.join(new_lines) + file_footer)
        print(f"Successfully created/updated '{output_path}'")


def plot_profile(data, profile, zoom_width = 80):
    """
    Plots both the full star and a zoomed in region around the point at which STIR and the progenitor are stitched together.
    
    Parameters:
        data (dict) :
            The combined STIR and MESA data containing paramaters useful drawing domains.
        profile (numpy array) :
            The profile to be plotted.
        zoom_width (int) :
            The full width of the zoomed in region around the stitch point.
    """

    # Most plots should use enclosed mass as the x-axis, except the enclosed mass plot
    if profile == "enclosed_mass":
        xaxis = np.arange(data["pre_masscut_profiles"].shape[0])
        xlabel = "zone"
    else:
        #xaxis = np.log10(data["pre_masscut_profiles"]["r"].values)
        xaxis = data["pre_masscut_profiles"]["enclosed_mass"].values
        #xlabel = "log(radius)"
        xlabel = "enclosed_mass (M_sun)"

    ylog = profile in log_plots
    ylabel = f"log({profile})" if ylog else profile
    yaxis = data["pre_masscut_profiles"][profile].values
    if ylog: yaxis = np.log10(yaxis)

    # Plot the entire curve of stitched data
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    axs[0].plot(xaxis, yaxis)
    axs[0].set_xlabel(xlabel)
    axs[0].set_ylabel(ylabel)
    axs[0].set_title(f"Full Profile")
    axs[0].axvspan(xaxis[0], xaxis[data["pns_masscut_index"]], color="red", alpha=0.2, ec=None, label = "PNS")
    axs[0].axvspan(xaxis[data["pns_masscut_index"]], xaxis[data["stir_domain_end"]], color="blue", alpha=0.2, ec=None, label = "STIR")
    axs[0].axvspan(xaxis[data["stir_domain_end"]], xaxis[-1], color="green", alpha=0.2, ec=None, label = "MESA")
    axs[0].legend()

    # Plot again but zoomed in on the stitch region
    zoom_left = max(0, data["stir_domain_end"] - zoom_width//2)
    zoom_right = min(data["pre_masscut_profiles"].shape[0], data["stir_domain_end"] + zoom_width//2)
    axs[1].plot(xaxis[zoom_left : zoom_right], yaxis[zoom_left : zoom_right])
    axs[1].axvspan(xaxis[zoom_left], xaxis[data["stir_domain_end"]], color="blue", alpha=0.2, ec=None)
    axs[1].axvspan(xaxis[data["stir_domain_end"]], xaxis[zoom_right], color="green", alpha=0.2, ec=None)
    axs[1].set_title(f"Interface Region")
    
    plt.show()


def shift_to_cell_edge(profile):
    """
    Takes values which are cell centered and shifts them to the edge of the cell.

    Parameters:
        profile (numpy array) :
            e.g. mass density

    Returns:
        shifted_data (numpy array) :
    """

    shifted_data = np.zeros_like(profile)
    shifted_data[:-1] = profile[:-1] + 0.5 * (profile[1:] - profile[:-1])
    shifted_data[-1] = profile[-1] + 0.5 * (profile[-1] - profile[-2])
    return shifted_data