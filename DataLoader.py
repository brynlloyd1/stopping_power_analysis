from gpaw import restart, setup_paths

import numpy as np
import pandas as pd

import logging
import os
import re

# imports for typing
from typing import Dict, List, Any
from numpy.typing import NDArray
from ase import Atoms
from gpaw import GPAW
from dataclasses import dataclass, field


logger = logging.getLogger(__name__)


@dataclass
class Data:
    """
    holds anything that can be loaded from data files,
    eg atoms/calc from loading gpw files,
    projectile positions and energies from loading csv's
    or electrond densities from loading npy files
    """

    # needs default values because some attributes will be left empty, depending on which data is being loaded
    atoms_list: List[Atoms] = field(default_factory=list)
    calc_list: List[GPAW] = field(default_factory=list)

    projectile_positions: NDArray[np.float64] = field(default_factory = lambda: np.empty((0,3)))
    projectile_kinetic_energies: NDArray[np.float64] = field(default_factory = lambda: np.array([]))

    electron_densities: NDArray[np.float64] = field(default_factory = lambda: np.array([]))
    cell: NDArray[np.float64] = field(default_factory = lambda: np.array([]))
    supercell_size: NDArray[np.float64] = field(default_factory = lambda: np.array([]))
    starting_position: NDArray[np.float64] = field(default_factory = lambda: np.array([]))
    direction: NDArray[np.float64] = field(default_factory = lambda: np.array([]))

class DataLoader:
    def __init__(self, directory: str):
        """
        (data loader doesn't store any data - all data is passed back to the data handler class)

        Parameters:
            directory (str)
        """
        self.directory: str = directory
        self.which_energies: List[str] = ["all"]
        self.which_timesteps: str = "all"   # not actually supported anymore to set timesteps, but its also not really needed

        self._energy_regex_pattern = re.compile(r"(\d+)k")
        self._timestep_regex_pattern = re.compile(r"k_step(\d+)")
        self._energy_timestep_regex_pattern = re.compile(r"(\d+)k_step(\d+)")

        # tells gpaw where to look for the custom setups that are sometimes used
        setup_paths.insert(0, "./custom_setups")

        self.stopping_data_loaded: bool = False
        self.density_data_loaded: bool = False

    def get_files(self, kind: str) -> Dict[str, List[str]]:
        """
        returns a dictionary of all files of a given type, keys are energies, values are lists of filenames

        Parameters:
            kind (str)

        Returns:
            Dict[str, List[str]]
        """
        all_files = os.listdir(self.directory)
        all_files_kind: List[str] = [f for f in all_files if f.endswith("." + kind)]

        filenames_temp = {}
        if kind == "csv" or kind == "npy":
            for filename in all_files_kind:
                match = self._energy_regex_pattern.search(filename)
                if not match:
                    continue
                energy = match.group(1)
                filenames_temp = self.append_to_dict(filenames_temp, energy, filename)


        elif kind == "gpw":
            timesteps_temp = {}
            for filename in all_files_kind:
                match = self._energy_timestep_regex_pattern.search(filename)
                if not match:
                    continue
                energy, timestep = match.group(1), match.group(2)
                filenames_temp = self.append_to_dict(filenames_temp, energy, filename)
                timesteps_temp = self.append_to_dict(timesteps_temp, energy, timestep)

            # sort files in ascending order according to the timestep
            for energy in filenames_temp.keys():
                paired = sorted(zip(map(int, timesteps_temp[energy]), filenames_temp[energy]))
                _, filenames_sorted = zip(*paired)
                filenames_temp[energy] = list(filenames_sorted)
        else:
            raise ValueError("kind must be 'csv', 'npy' or 'gpw'")

        # rename keys
        filenames_temp = {f"{key} keV": value for key, value in filenames_temp.items()}
        # sort the energies so that they are also in ascending order
        sorted_items = sorted(filenames_temp.items(), key=lambda item: int(item[0].split()[0]))
        filename_dict = dict(sorted_items)

        return filename_dict

    def load_data(self, force_load_gpw: bool = False, force_write_csv: bool = False) -> Dict[str, Data]:
        """
        If .csv files exist, loads position and kinetic energy data from them
        if not, uses gpaw.restart to load in atoms, calc
        Data that is loaded is stored in an instance of Data

        Parameters:
            force_load_gpw (bool): would set to True if you need to access atoms and calc from gpaw.restart
            force_write_csv (bool): would set to True if the data has changed, but the csv files still exist

        Returns:
            all_data (Dict[str, Data]): key is energy, value is Data instance
        """

        all_gpw_files = self.get_files("gpw")
        all_csv_files = self.get_files("csv")

        if self.which_energies == ["all"]:
            self.which_energies = all_gpw_files.keys()

        all_data: Dict[str, Data] = {}
        for energy in self.which_energies:
            write_flag = False

            data = Data()

            if energy not in all_csv_files.keys() or force_load_gpw:
                if not force_load_gpw:
                    logger.info(f"csv not found for {energy} in {os.path.basename(self.directory.rstrip("/"))}, loading gpw files...")
                write_flag = not (energy in all_csv_files.keys()) or force_write_csv  # only want to write if the csv files don't already exist

                for filename in all_gpw_files[energy]:
                    try:
                        atoms, calc = restart(self.directory + filename)
                    except Exception as e:
                        logging.warning(f"failed to read file {filename}, skipping: {e}")
                        continue
                    data.atoms_list.append(atoms)
                    data.calc_list.append(calc)

                    data.projectile_positions = np.vstack([data.projectile_positions, atoms.get_positions()[-1]])


                    # this is the kinetic energy of the entire system
                    # data.projectile_kinetic_energies = np.append(data.projectile_kinetic_energies, atoms.get_kinetic_energy())

                    # this is how ase calculates kinetic energy
                    projectile_momentum = atoms.get_momenta()[-1]
                    projectile_velocity = atoms.get_velocities()[-1]
                    projectile_kinetic_energy = 0.5 * np.vdot(projectile_momentum, projectile_velocity)
                    data.projectile_kinetic_energies = np.append(data.projectile_kinetic_energies, projectile_kinetic_energy)

                    data.cell = np.diag(atoms.get_cell())

            else:
                write_flag = False
                for filename in all_csv_files[energy]:
                    df = pd.read_csv(self.directory + filename, comment="#")
                    data.projectile_positions = df[["projectile x [A]", "projectile y [A]", "projectile z [A]"]].to_numpy()
                    data.projectile_kinetic_energies = df["projectile KE [eV]"].to_numpy()

                    metadata = self.read_metadata(filename)
                    # from DataLoader.write_csv, the second line will contain information about Atoms.cell()
                    data.cell = np.fromstring(metadata[1].strip("# cell: ")[1: -2], sep=" ")

                    data.supercell_size = np.fromstring(metadata[2].strip().replace("# supercell size: [", "").replace("]", ""), sep=" ")
                    data.supercell_size = tuple(data.supercell_size.astype(int))
                    data.starting_position = np.fromstring(metadata[3].strip().replace("# starting_position: [", "").replace("]", ""), sep=" ")
                    data.direction = np.fromstring(metadata[4].strip().replace("# direction: [", "").replace("]", ""), sep=" ")

            if write_flag:
                self.write_csv(energy, data, overwrite_flag = force_write_csv)

            # all_data = self.append_to_dict(all_data, energy, data)
            all_data[energy] = data

        self.stopping_data_loaded = True
        return all_data

    def write_csv(self, energy: str, data: Data, overwrite_flag: bool = False):
        """
        writes projectile positions and kinetic energies to a csv file

        Parameters:
            energy (str): e.g. "40 keV"
            data (Data)
            overwrite_flag (bool): prevents accidental overwriting of existing csv
        """

        filename = f"Al_stopping_{energy.rstrip(' keV')}k.csv"
        if not overwrite_flag:
            if filename in os.listdir(self.directory):
                logger.info(f"File {filename} already exists, must specify overwrite_flag = True to overwrite.")
                return

        logger.info(f"Writing {filename}")

        timesteps = list(range(1, len(data.projectile_positions) + 1))
        df = pd.DataFrame(data = {
            "timestep": timesteps,
            "projectile x [A]": data.projectile_positions[:, 0],
            "projectile y [A]": data.projectile_positions[:, 1],
            "projectile z [A]": data.projectile_positions[:, 2],
            "projectile KE [eV]" : data.projectile_kinetic_energies
        }).to_csv(self.directory + filename, index=False)

        # WRITE ADDITIONAL INFORMATION TO THE BOTTOM OF THE CSV FILE
        # "#" to start the line makes it easy to load with pandas
        metadata = [
            "# --- additional info ---",
            f"# cell: {data.cell}",
            f"# supercell size: {int(data.cell / 4.05)}",   # for trajectory presampling purposes
            f"# starting_position: {data.projectile_positions[0]}",
            f"# direction: {data.atoms_list[0].get_velocities()[-1] / np.linalg.norm(data.atoms_list[0].get_velocities()[-1])}",
        ]

        with open(os.path.join(self.directory, filename), "a") as f:
            for line in metadata:
                f.write(line + "\n")

    def read_metadata(self, filename: str) -> List[str]:
        """
        Parameters:
            path_to_csv_file (str)
        Returns:
            metadata (List[str])
        """
        metadata = []
        with open(os.path.join(self.directory, filename), "r") as f:
            for line in f:
                if line.startswith("#"):
                    metadata.append(line)

        return metadata


    @staticmethod
    def append_to_dict(dictionary: Dict[str, List], key: str, value: List[Any]) -> Dict[str, List]:
        """
        appends to a dictionary, where the value is a list of elements

        Parameters:
            dictionary (dict[str, __])
            key (str)
            value (List[Any])

        Returns:
            dictionary (dict[str, List[Any]])
        """

        if key not in dictionary:
            dictionary[key] = [value]
        else:
            dictionary[key].append(value)

        return dictionary


    def load_densities(self) -> Dict[str, Data]:
        all_gpw_files = self.get_files("gpw")
        all_npy_files = self.get_files("npy")

        if self.which_energies == ["all"]:
            self.which_energies = all_gpw_files.keys()

        all_data = {}
        for energy in self.which_energies:
            write_flag = False

            data = Data()
            electron_densities_temp = []
            if not energy in all_npy_files.keys():
                logger.info(f"npy not found for {energy} in {os.path.basename(self.directory.rstrip("/"))}, loading gpw files...")
                write_flag = True
                for filename in all_gpw_files[energy]:
                    atoms, calc = restart(self.directory + filename)
                    data.atoms_list.append(atoms)
                    data.calc_list.append(calc)
                    electron_densities_temp.append(calc.get_all_electron_density())

                data.electron_densities = np.stack(electron_densities_temp)

            else:
                write_flag = False
                for filename in all_npy_files[energy]:
                    data.electron_densities = np.load(self.directory + filename)

            if write_flag:
                self.write_npy(energy, data)

            all_data[energy] = data

        self.density_data_loaded = True
        return all_data

    def write_npy(self, energy: str, data: Data, overwrite_flag: bool = False):
        """
        writes electron densities to a npy file

        Parameters:
            energy (str): e.g. "40 keV"
            data (Data)
            overwrite_flag (bool): prevents accidental overwriting of existing npy file
        """
        filename = f"Al_stopping_{energy.rstrip(' keV')}k.npy"
        if not overwrite_flag:
            if filename in os.listdir(self.directory):
                logger.info(f"File {filename} already exists, must specify overwrite_flag = True to overwrite.")
                return

        logger.info(f"Writing {filename}")
        np.save(self.directory + filename, data.electron_densities)
