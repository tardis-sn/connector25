# general imports
import logging

# TARDIS imports for the gamma-ray code
from tardis.io.atom_data import AtomData
from tardis.model import SimulationState
from tardis.io.configuration import config_reader
from tardis.energy_input.main_gamma_ray_loop import (run_gamma_ray_loop, 
                                                     get_effective_time_array)
from tardis.energy_input.gamma_ray_channel import (create_isotope_dicts, 
            create_inventories_dict, calculate_total_decays, 
            create_isotope_decay_df, time_evolve_cumulative_decay)


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class TARDISHE:
    def __init__(self, atom_data, configuration, config_type="csvy"):

        if config_type == "csvy":
            self.simulation_state = SimulationState.from_csvy(
                configuration, atom_data
            )
        else:
            self.simulation_state = SimulationState.from_config(
                configuration, atom_data
            )

        self.gamma_ray_lines = atom_data.decay_radiation_data

        self.shell_masses = (
            self.simulation_state.volume * self.simulation_state.density
        )

        self.isotopic_mass_fraction = (
            self.simulation_state.composition.isotopic_mass_fraction
        )

    def calculate_total_decays(self, time_start, time_end):
        """
        Calculate the total number of decays for each isotope in the simulation
        state between t_start and t_end for each shell.

        Parameters
        ----------
        time_start : float
            The start time in days.
        time_end : float
            The end time in days.

        Returns
        -------
        total_decays : dataframe
            A dataframe containing the total number of decays for each isotope
            in the simulation state between t_start and t_end for each shell.
        """
        isotopes = create_isotope_dicts(
            self.isotopic_mass_fraction, self.shell_masses
        )
        inventories = create_inventories_dict(isotopes)

        total_decays = calculate_total_decays(
            inventories, time_end - time_start
        )

        return total_decays
    
    def decay_isotopes_expanded(self, total_decays):
        """
        Expand the total decays dataframe to include the decay energy for each gamma-X ray
        transition in each shell between t_start and t_end.

        Parameters
        ----------
        total_decays : dataframe
            A dataframe containing the total number of decays for each isotope
            in the simulation state between t_start and t_end for each shell.
        
        Returns
        -------
        decay_isotopes : dataframe
            A dataframe containing the total number of decays for each isotope
            in the simulation state between t_start and t_end for each shell
            expanded to include the decay energy for each gamma-X ray transition.
        """
        return create_isotope_decay_df(total_decays, self.gamma_ray_lines)
    
    def time_evolve_cumulative_decay_expanded(self, times):
        """
        Time evolve the cumulative decay for each isotope in the simulation state
        between t_start and t_end for each shell.

        Parameters
        ----------
        times : array
            An array of times in days.
        
        Returns
        -------
        cumulative_decay : dataframe
            A dataframe containing the cumulative decay for each isotope, each transition 
            each shell for all time steps between t_start and t_end.
        """
        return time_evolve_cumulative_decay(self.isotopic_mass_fraction, self.shell_masses,
                                            self.gamma_ray_lines, times)
    
    def get_times(self, time_start, time_end, time_space, time_steps):
        """
        Create an array of times between t_start and t_end.

        Parameters
        ----------
        time_start : float
            The start time in days.
        time_end : float
            The end time in days.
        number_of_steps : int
            The number of time steps to take between t_start and t_end.
        
        Returns
        -------
        times : array
            An array of times in days.
        """
        return get_effective_time_array(time_start, time_end, time_space, time_steps)
    
    def run(self, time_start, time_end, number_of_packets,
            time_steps, time_space, seed, fp, spectrum_bins, grey_opacity=-1):
        """
        Run the gamma-ray transport simulation.
        
        Parameters
        ----------
        time_start : float
            The start time in days.
        time_end : float
            The end time in days.
        times : array
            An array of times in days.
        number_of_packets : int
            The number of packets to simulate.
        time_steps : int
            The number of time steps to take between t_start and t_end.
        time_space : string
            The time spacing between time steps.
        seed : int
            The random seed for the simulation.
        fp : str
            Positronium fraction
        spectrum_bins : int
            The number of bins in the spectrum.
        grey_opacity : float
            The grey opacity of the simulation.
        """
        times, effective_times = self.get_times(time_start, time_end, time_space, time_steps)
        total_decays = self.calculate_total_decays(time_start, time_end)
        decay_isotopes = self.decay_isotopes_expanded(total_decays)
        decay_over_time = self.time_evolve_cumulative_decay_expanded(times)


        (self.escape_energy, self.escape_energy_cosi, 
         self.packets_escaped, self.gamma_ray_deposited_energy, 
         self.total_deposited_energy, self.positron_energy,
         self.escape_luminosity
        ) = run_gamma_ray_loop(self.simulation_state, decay_isotopes, decay_over_time,
                                  number_of_packets, times, effective_times, seed, fp,
                                  spectrum_bins, grey_opacity)
        


def main():
    # Load the atom data
    atom_data = AtomData.from_hdf("kurucz_cd23_chianti_H_He.h5")

    # Load the configuration file
    configuration = config_reader.Configuration.from_yaml("6.4_stripped_1foe_0.03Ni/tardis_configs/tardis_config_mesa_stella_1.yml")

    # Create an instance of the TARDISHEWorkflow class
    tardis_he_workflow = TARDISHE(atom_data, configuration)

    # Run the simulation
    tardis_he_workflow.run(
        time_start=0.1,
        time_end=100,
        number_of_packets=1000000,
        time_steps=10,
        time_space="log",
        seed=42,
        fp=0.5,
        spectrum_bins=100,
        grey_opacity=-1
    )

if __name__ == "__main__":
    main()