from gpaw import restart

import numpy as np
import pandas as pd

import os
import re
from typing import Dict

class Data:
    def __init__(self):
        self.atoms_list = []
        self.calc_list = []
        self.projectile_positions = np.empty((0,3))
        self.projectile_kinetic_energies = np.array([])
        self.electron_densities = np.array([])

class DataLoader:
    def __init__(self, directory):
        """
        (data loader doesn't store any data - all data is passed back to the data handler class)
        """
        self.directory = directory
        self.which_energies = ["all"]
        self.which_timesteps = "all"

        self._energy_regex_pattern = re.compile(r"(\d+)k")
        self._timestep_regex_pattern = re.compile(r"k_step(\d+)")
        self._energy_timestep_regex_pattern = re.compile(r"(\d+)k_step(\d+)")


    def get_files(self, kind):
        """
        returns a dictionary of all files of a given type, separated by their energies

        Parameters:
            kind (str)
        """
        all_files = os.listdir(self.directory)
        all_files_kind = [f for f in all_files if f.endswith("." + kind)]

        filenames_temp = {}
        if kind == "csv":
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
            raise ValueError("kind must be 'csv' or 'gpw'")

        # rename keys
        filenames_temp = {f"{key} keV": value for key, value in filenames_temp.items()}
        # sort the energies so that they are also in ascending order
        sorted_items = sorted(filenames_temp.items(), key=lambda item: int(item[0].split()[0]))
        filename_dict = dict(sorted_items)

        return filename_dict

    def load_data(self, force_load_gpw=False):
        """
        If .csv files exist, loads position and kinetic energy data from them
        if not, uses gpaw.restart to load in atoms, calc
        Data that is loaded is stored in an instance of Data

        Parameters:
            force_load_gpw (bool)

        Returns:
            data_dict (Dict[str, Data]): key is energy, value is Data instance
        """

        # check if existence of csv files
        all_gpw_files = self.get_files("gpw")
        all_csv_files = self.get_files("csv")

        if self.which_energies == ["all"]:
            self.which_energies = all_gpw_files.keys()

        all_data = {}
        for energy in self.which_energies:
            write_flag = False

            data = Data()

            if not energy in all_csv_files.keys() or force_load_gpw:
                write_flag = not (energy in all_csv_files.keys())  # only want to write if the csv files don't already exist
                for filename in all_gpw_files[energy]:
                    atoms, calc = restart(self.directory + filename)
                    data.atoms_list.append(atoms)
                    data.calc_list.append(calc)
                    data.projectile_positions = np.vstack([data.projectile_positions, atoms.get_positions()[-1]])
                    data.projectile_kinetic_energies = np.append(data.projectile_kinetic_energies, atoms.get_kinetic_energy())

            else:
                write_flag = False
                for filename in all_csv_files[energy]:
                    df = pd.read_csv(self.directory + filename)
                    data.projectile_positions = df[["projectile x [A]", "projectile y [A]", "projectile z [A]"]].to_numpy()
                    data.projectile_kinetic_energies = df["projectile KE [eV]"].to_numpy()

            if write_flag:
                self.write_csv(energy, data)

            # all_data = self.append_to_dict(all_data, energy, data)
            all_data[energy] = data

        return all_data

    def write_csv(self, energy, data):
        """
        writes projectile positions and kinetic energies to a csv file

        Parameters:
            energy (str): e.g. "40 keV"
            data (Data)
        """

        filename = f"Al_stopping_{energy.rstrip(' keV')}k"

        timesteps = list(range(1, len(data.projectile_positions) + 1))
        df = pd.DataFrame(data = {
            "timestep": timesteps,
            "projectile x [A]": data.projectile_positions[:, 0],
            "projectile y [A]": data.projectile_positions[:, 1],
            "projectile z [A]": data.projectile_positions[:, 2],
            "projectile KE [eV]" : data.projectile_kinetic_energies
        }).to_csv(self.directory + filename + ".csv", index=False)


    def append_to_dict(self, dictionary: Dict[str, [any]], key, value):
        """
        appends to a dictionary, where the value is a list of elements

        Parameters:
            dictionary (dict[str, __])

        Returns:
            dictionary (dict[str, __])
        """

        if key not in dictionary:
            dictionary[key] = [value]
        else:
            dictionary[key].append(value)

        return dictionary

    def check_for_npy(self, directory, energy):
        """checks if a .npy containing electron density data exists"""
        trajectory_name = os.path.basename(directory.rstrip("/"))
        npy_filename = f"{trajectory_name}_{energy}"
        return npy_filename in os.listdir(directory)

    def load_from_npy(self, directory, energy):
        """loads .npy containing electron density data exists"""
        trajectory_name = os.path.basename(directory.rstrip("/"))
        npy_filename = f"{trajectory_name}_{energy}"
        return np.load(directory+npy_filename)

    def save_to_npy(self, directory, energy, electron_density_list):
        """writes .npy containing electron density data exists"""
        trajectory_name = os.path.basename(directory.rstrip("/"))
        npy_filename = f"{trajectory_name}_{energy}"
        np.save(f"{trajectory_name}_{energy}", np.array(electron_density_list))

