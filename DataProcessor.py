import numpy as np
from DataLoader import Data
from typing import Dict, List


class Fit:
    """
    class that stored information about fits together
    stored the fit, the covariance matrix,
    and how the original data was cropped
    """
    def __init__(self, fit, cov, crop):
        self.fit = fit
        self.cov = cov
        self.crop = crop



class DataProcessor:
    def __init__(self):
        pass


    # def calculate_stopping_powers(self, atoms_dict, crop=[None, None]):
    def calculate_stopping_powers(self,
                                  all_data: Dict[str, Data],
                                  crop: List[int | None] = [None, None]):
        """
        performs linear fits to kinetic energy data to extract averaged stopping powers
        Parameters:
            all_data (Dict[str, Data])
            crop (List[int | None])

        Returns:
            fits_information (Dict[str, Fit])
        """
        items = all_data.items()
        projectile_positions = {key: value.projectile_positions for key,value in all_data.items()}
        projectile_kinetic_energies = {key: value.projectile_kinetic_energies for key,value in all_data.items()}

        fits_information: Dict[str, Fit] = {}

        for i in range(len(projectile_kinetic_energies.keys())):
            # get raw data
            energy, kinetic_energies = list(projectile_kinetic_energies.items())[i]
            kinetic_energies = np.array(kinetic_energies) * 1e-3  # convert from eV to keV
            _, positions = list(projectile_positions.items())[i]

            initial_position = positions[0]
            distance_travelled = np.array([np.linalg.norm(position - initial_position) for position in positions])

            distance_per_timestep = np.diff(distance_travelled)[0]
            n_timesteps = len(distance_travelled)

            distance_travelled = distance_travelled[crop[0]: n_timesteps - (crop[1] or 0)]
            kinetic_energies = kinetic_energies[crop[0]: n_timesteps - (crop[1] or 0)]
            # (None or 0) is 0, (int or 0) is int
            print(f"""
            For energy {energy}:
            cropping {(crop[0] or 0)} timesteps off the start = {(crop[0] or 0) * distance_per_timestep} Angstroms
            cropping {(crop[1] or 0)} timesteps off the end = {(crop[1] or 0) * distance_per_timestep} Angstroms
            {len(distance_travelled)} timesteps remaining = {len(distance_travelled) * distance_per_timestep:.1f} Angstroms""")

            fit, cov = np.polyfit(distance_travelled, kinetic_energies, 1, cov=True)
            fit_info = Fit(fit, cov, crop)
            fits_information[energy] = fit_info

        return fits_information
