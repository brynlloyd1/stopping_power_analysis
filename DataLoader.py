from gpaw import restart
from ase.units import eV, _amu

import numpy as np
import pandas as pd

import os
import re
from typing import Dict, List

class Data:
    def __init__(self):
        self.atoms_list = []
        self.calc_list = []
        self.projectile_positions = np.empty((0,3))
        self.projectile_kinetic_energies = np.array([])
        self.electron_densities = None
        self.cell = None

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
            all_data (Dict[str, Data]): key is energy, value is Data instance
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
                if not force_load_gpw:
                    print(f"csv not found for {energy} in {os.path.basename(self.directory.rstrip("/"))}, loading gpw files...")
                write_flag = not (energy in all_csv_files.keys())  # only want to write if the csv files don't already exist
                for filename in all_gpw_files[energy]:
                    atoms, calc = restart(self.directory + filename)
                    data.atoms_list.append(atoms)
                    data.calc_list.append(calc)

                    data.projectile_positions = np.vstack([data.projectile_positions, atoms.get_positions()[-1]])

                    # TODO: THIS GIVES KINETIC ENERGY OF ENTIRE SYSTEM
                    #       CODE BELOW CALCULATES FROM VELOCITIES, BUT THERE IS A MISTAKE IN THE UNITS SOMEWHERE
                    data.projectile_kinetic_energies = np.append(data.projectile_kinetic_energies, atoms.get_kinetic_energy())
                    # projectile_velocities = atoms.get_velocities()[-1] * 100   # A/ps -> m/s
                    # projectile_mass = atoms.get_masses()[-1] * _amu   # atomic mass units -> kg
                    # kinetic_energy = 0.5 * projectile_mass * np.sum(projectile_velocities**2) / eV   # Joules -> eV
                    # data.projectile_kinetic_energies = np.append(data.projectile_kinetic_energies, kinetic_energy)

                    data.cell = np.diag(atoms.get_cell())

            else:
                write_flag = False
                for filename in all_csv_files[energy]:
                    df = pd.read_csv(self.directory + filename, comment="#")
                    data.projectile_positions = df[["projectile x [A]", "projectile y [A]", "projectile z [A]"]].to_numpy()
                    data.projectile_kinetic_energies = df["projectile KE [eV]"].to_numpy()

                    metadata = []
                    with open(os.path.join(self.directory, filename), "r") as f:
                        for line in f:
                            if line.startswith("#"):
                                metadata.append(line)

                    data.cell = np.fromstring(metadata[1].strip("# cell: ")[1: -2], sep=" ")

            if write_flag:
                self.write_csv(energy, data)

            # all_data = self.append_to_dict(all_data, energy, data)
            all_data[energy] = data

        return all_data

    def write_csv(self, energy, data, overwrite_flag = False):
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
                print(f"File {filename} already exists, must specify overwrite_flag = True to overwrite.")
                return

        print(f"Writing {filename}")

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
            f"# cell: {data.cell}"
        ]

        with open(os.path.join(self.directory, filename), "a") as f:
            for line in metadata:
                f.write(line + "\n")



    def append_to_dict(self, dictionary: Dict[str, List], key, value):
        """
        appends to a dictionary, where the value is a list of elements

        Parameters:
            dictionary (dict[str, __])
            key (str)
            value (list)

        Returns:
            dictionary (dict[str, __])
        """

        if key not in dictionary:
            dictionary[key] = [value]
        else:
            dictionary[key].append(value)

        return dictionary


    def load_densities(self):
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
                print(f"npy not found for {energy} in {os.path.basename(self.directory.rstrip("/"))}, loading gpw files...")
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

        return all_data

    def write_npy(self, energy, data, overwrite_flag = False):
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
                print(f"File {filename} already exists, must specify overwrite_flag = True to overwrite.")
                return

        print(f"Writing {filename}")
        np.save(self.directory + filename, data.electron_densities)
