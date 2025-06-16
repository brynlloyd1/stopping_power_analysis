from gpaw import restart

import numpy as np

import os
import re
import json

class DataLoader:
    def __init__(self, directory):
        """
        (data loader doesnt store any data - all data is passed back to the data handler class)
        """
        self.directory = directory
        self.which_energies = ["all"]
        self.which_timesteps = "all"

        self._energy_regex_pattern = re.compile(r"(\d+)k_step")
        self._timestep_regex_pattern = re.compile(r"k_step(\d+)")
        self._energy_timestep_regex_pattern = re.compile(r"(\d+)k_step(\d+)")

        self.all_gpw_files = self.get_all_gpaw_files()


    def get_all_gpaw_files(self):
        """
        gets .gpw files from a given directory

        Returns:
        dict: key is energy (in keV), value is a list of filenames, sorted by ascending timestep
        """

        # get all files, and throw any that arent .gpw
        all_files = os.listdir(self.directory)
        all_gpw_files = [f for f in all_files if f.endswith(".gpw")]

        filenames_temp = {}
        timesteps_temp = {}

        # regex to extract energy and timestep
        for filename in all_gpw_files:
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

        # rename keys
        filename_dict = {f"{key} keV": value for key, value in filenames_temp.items()}

        # sort the energies so that they are also in ascending order
        sorted_items = sorted(filename_dict.items(), key=lambda item: int(item[0].split()[0]))
        filename_dict = dict(sorted_items)

        return filename_dict

    def set_which_energies(self, which_energies):
        self.which_energies = which_energies

    def set_which_timesteps(self, which_timesteps):
        self.which_timesteps = which_timesteps

    def load(self, print_files=False):
        """
        loads data from .gpw files

        Parameters:
        directory (str): path to file containing .gpw files
        which_energies (list[str]): list of energies to load in. ["all"] will load in all energies
        which_timesteps (str): which timesteps to load in. Options are "all", or "::10" (like the slice)

        Returns:
        atoms_dict (Dict[str, list[Atoms]]): Key is energy (eg. "40 keV"), value is a list of Atoms objects
        calc_dict (Dict[str, list[Calc]]): Key is energy, value is a list of GPAW objects
        """

        if self.which_energies == ["all"]:
            self.which_energies = self.all_gpw_files.keys()

        # check that all specified energies are in all_files.keys()
        for energy in self.which_energies:
            if energy in self.all_gpw_files.keys():
                continue
            else:
                raise KeyError(f"{energy} not found in this directory")

        atoms_dict = {}
        calc_dict = {}
        for energy in self.which_energies:
            for filename in self.all_gpw_files[energy]:

                if self.which_timesteps == "::10":
                    match = self._timestep_regex_pattern.search(filename)
                    if not match:
                        continue
                    timestep = int(match.group(1))
                    if timestep % 10 != 0:
                        continue
                    else:
                        self.read_gpw_to_dict(atoms_dict, calc_dict, energy, filename)

                else:
                    self.read_gpw_to_dict(atoms_dict, calc_dict, energy, filename)

        if print_files:
            print(self.directory)
            print(json.dumps(self.all_gpw_files, indent=4))
        return atoms_dict, calc_dict

    def append_to_dict(self, dictionary, key, value):
        """
        appends to a dictionary, where the value is a list of elements

        Paramters:
        dictionary (dict[str, __])

        Returns:
        dictionary (dict[str, __])
        """

        if key not in dictionary:
            dictionary[key] = [value]
        else:
            dictionary[key].append(value)

        return dictionary

    # can probably move the append_to_dict funtion to be a method of this class
    def read_gpw_to_dict(self, atoms_dict, calc_dict, energy, filename):
        """uses gpaw.restart to load .gpw files. Writes data to dictionaries"""
        atoms, calc = restart(self.directory + filename)
        self.atoms_dict = self.append_to_dict(atoms_dict, energy, atoms)
        self.calc_dict = self.append_to_dict(calc_dict, energy, calc)


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
