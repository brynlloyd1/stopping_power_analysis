import numpy as np

from dataclasses import field
# imports for typing
from typing import Dict, List, Optional
from DataLoader import Data
from numpy.typing import NDArray
from dataclasses import dataclass

@dataclass
class Fit:
    """stores fit, covariance matrix, and how original data was cropped"""
    fit: NDArray[np.float64]
    cov: NDArray[np.float64]
    crop: List[int | None]

@dataclass
class ChargeStateData:
    density_around_projectile: NDArray[np.float64]
    crop: List[int | None]
    radius: int
    offset: int


class DataProcessor:
    """
    is really just a collection of methods for processing data

    Methods:
        calculate_stopping_powers(all_data: Dict[str, Data], crop: List[int | None]) -> Dict[str, Fit]
            calculates stopping powers by performing a linear fit to (cropped) projectile kinetic energies
        get_electrons_around_proton(all_data: Dict[str, Data]) -> Dict[str, ChargeStateData]
            sums electron density in the region of the projectile for analysis of the projectile charge state
    """

    @staticmethod
    def calculate_stopping_powers(all_data: Dict[str, Data],
                                  crop: List[int | None] = [None, None]) -> Dict[str, Fit]:
        """
        performs linear fits to kinetic energy data to extract averaged stopping powers
        Parameters:
            all_data (Dict[str, Data])
            crop (List[int | None])

        Returns:
            fits_information (Dict[str, Fit])
        """

        # replaced crop: List[int | None] = [None, None] (in case it breaks)

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

    @staticmethod
    def get_electrons_around_proton(all_data: Dict[str, Data],
                                    parameters: Dict[str, int]) -> Dict[str, ChargeStateData]:
        """
        for each energy simulated for a given trajectory, sums electron density around the projectile position
        to get the charge state of the projectile over time

        Parameters:
            all_data (Dict[str, Data])
            parameters (Dict[str, int])

        Returns:
            density_dict (Dict[str, ChargeStateData])
        """

        size = parameters["size"]
        offset = parameters["offset"]
        crop_low = parameters["crop_low"]
        crop_high = parameters["crop_high"]


        charge_state_data_dict = {}
        for energy, data in all_data.items():


            # check that the required data has been loaded in first
            if data.projectile_positions.size == 0:
                raise AttributeError(f"projectile positions not loaded for energy {energy}")
            if data.electron_densities.size == 0:
                raise AttributeError(f"electron densities not loaded for energy {energy}")


            projectile_positions = data.projectile_positions
            electron_densities = data.electron_densities
            electron_density_shape = np.shape(electron_densities[0])
            cell = data.cell

            projectile_position_indices = DataProcessor.find_projectile(projectile_positions, electron_densities[0], cell)

            density_around_projectile = []
            for t, electron_density in enumerate(electron_densities):
                density_sum = 0
                projectile_position_index = projectile_position_indices[t]
                start = projectile_position_index - size
                stop = projectile_position_index + size

                for i in range(start[0], stop[0]):
                    for j in range(start[1], stop[1]):
                        for k in range(start[2], stop[2]):

                            pos = np.array([i, j, k])
                            dist_to_projectile = np.linalg.norm(pos - projectile_position_indices[t])
                            if dist_to_projectile > size:
                                continue

                            pos_inside_cell = pos % electron_density_shape
                            density_sum += electron_density[pos_inside_cell[0], pos_inside_cell[1], pos_inside_cell[2]]

                density_around_projectile.append(density_sum)

            # TODO: remove crop and offset, but for now just set both to zero
            charge_state_data = ChargeStateData(np.array(density_around_projectile), [0, 0], size, offset)
            charge_state_data_dict[energy] = charge_state_data

        return charge_state_data_dict


    @staticmethod
    def find_projectile(projectile_positions: NDArray[np.float64],
                        electron_density: NDArray[np.float64],
                        cell: NDArray[np.float64]) -> NDArray[np.float64]:
        projectile_position_indices = np.array([(projectile_position % cell) / cell * np.shape(electron_density) for projectile_position in projectile_positions], dtype="int64")
        return projectile_position_indices
