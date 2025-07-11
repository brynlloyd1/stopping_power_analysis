from unittest import loader

from DataHandler import DataHandler

from dataclasses import field
# imports for typing
from typing import List, Dict, Optional
from DataLoader import Data
from DataProcessor import DataProcessor, Fit, ChargeStateData
from DataVisualiser import DataVisualiser
from DataLoader import DataLoader

import logging

logging.getLogger(__name__)

class EhrenfestAnalysis:
    def __init__(self):
        self.data_handlers = {}

    def initialise_analysis(self, directory):
        """
        creates an instance of DataHandler, which is stored in self.dat_handlers

        Parameters:
        directory (str): path to directory containing data

        Returns:
        DataHandler.trajectory_name (str): name given to this trajectory
        """
        data_handler = DataHandler(directory)
        self.data_handlers[data_handler.trajectory_name] = data_handler
        return data_handler.trajectory_name

    def set_name(self, old_name: str, new_name: str):
        """changes the name of a trajectory, which is stored by the DataHandler instance"""
        if new_name in self.data_handlers.keys():
            logging.warning(f"The name {new_name} already exists, returning without overwriting")
            # raise KeyError(f"{new_name} already exists")
            return

        self.data_handlers[old_name].trajectory_name = new_name
        self.data_handlers[new_name] = self.data_handlers[old_name].pop()

    ################
    # LOADING DATA #
    ################


    def load_data(self,
                  trajectory_name: str,
                  force_load_gpw: bool = False,
                  force_write_csv: bool = False):
        """
        instructs the instance of DataLoader stored in the DataHandler instance to load in the data

        Parameters:
            trajectory_name (str)
            force_load_gpw (bool, optional): Can force loading of gpw files (vs loading csv files), e.g. if you want to use electron_densities
            force_write_csv (bool, optional): Can force writing csv files, eg if stopping data has changed
        """

        data_handler: DataHandler = self.data_handlers[trajectory_name]
        # data_handler.atoms_dict, data_handler.calc_dict = data_handler.data_loader.load_data()
        data_handler.all_data = data_handler.data_loader.load_data(force_load_gpw, force_write_csv)

    def load_densities(self,
                       trajectory_name: str):
        """
        loads electron density data either by using gpaw.restart or by loading in .npy files
        """

        data_handler: DataHandler = self.data_handlers[trajectory_name]
        all_data_temp = data_handler.data_loader.load_densities()

        # have to add this to all_data differently, depending on whether it already contains
        # Data instances and so has attributes to modify, or need to create the instances
        if not data_handler.all_data:
            data_handler.all_data = all_data_temp
        else:
            for energy in data_handler.all_data.keys():
                data_handler.all_data[energy].electron_densities = all_data_temp[energy].electron_densities



    def set_which_energies(self, trajectory_name: str, which_energies: List[str]):
        data_handler: DataHandler = self.data_handlers[trajectory_name]
        data_loader = data_handler.data_loader
        data_loader.which_energies = which_energies



    #########################
    # STOPPING POWER THINGS #
    #########################


    def calculate_stopping_curve(self,
                                 trajectory_name: str,
                                 crop: List[int | None] = [11, None]):
        """tells the DataProcessor instance to perform fits to kinetic energy data"""

        handler: DataHandler = self.data_handlers[trajectory_name]

        if not handler.data_loader.stopping_data_loaded:
            self.load_data(trajectory_name)

        processor: DataProcessor = handler.data_processor
        # fits_information = processor.calculate_stopping_powers(handler.atoms_dict, crop=crop)
        handler.fits = processor.calculate_stopping_powers(handler.all_data, crop=crop)


    def view_fits(self, trajectory_name: str):
        """plots data and fits which comparison to Monte Carlo stopping powers"""
        handler = self.data_handlers[trajectory_name]
        loader = handler.data_loader
        if not loader.stopping_data_loaded:
            loader.load_data()



        visualiser = handler.data_visualiser
        visualiser.plot_all_fits(trajectory_name, handler.all_data, handler.fits)

    def compare_fits(self, trajectory_names: List[str], energy: str):
        comparison_data: Dict[str, Data] = {
            trajectory_name: self.data_handlers[trajectory_name].all_data[energy] for trajectory_name in trajectory_names
        }
        self.data_handlers[trajectory_names[0]].data_visualiser.compare_fits(trajectory_names, energy, comparison_data)



    def compare_to_montecarlo(self, trajectory_names: List[str]):
        """creates a dictionary of trajectory_name : fits for that trajectory, and passes to visualiser"""

        all_fit_info = {trajectory_name: self.data_handlers[trajectory_name].fits for trajectory_name in trajectory_names}
        self.data_handlers[trajectory_names[0]].data_visualiser.montecarlo_comparison(all_fit_info)



    ###########################
    # ELECTRON DENSITY THINGS #
    ###########################


    def visualise_electron_density(self, trajectory_name: str, energy: str):
        """DataVisualiser.visualise_electron_density called"""

        handler = self.data_handlers[trajectory_name]
        loader: DataLoader = handler.data_loader
        visualiser = handler.data_visualiser

        if not loader.density_data_loaded:
            self.load_densities(trajectory_name)

        electron_density_list = handler.all_data[energy].electron_densities
        visualiser.visualise_electron_density(electron_density_list)


    def analyse_proton_charge_state(self, trajectory_name: str, parameters: Dict[str, int]):
        """
        data processor calculates the charge state of the proton by summing the electron density around the proton up to a cutoff radius
        data visualiser plots the proton charge state over time
        and plots time-averaged charge state for each energy in a given trajectory
        Parameters:
            trajectory_name (str)
            parameters (Dict[str, int])
        """
        data_handler: DataHandler = self.data_handlers[trajectory_name]
        processor: DataProcessor = data_handler.data_processor
        visualiser: DataVisualiser = data_handler.data_visualiser

        if not data_handler.data_loader.stopping_data_loaded:
            self.load_data(trajectory_name)
        if not data_handler.data_loader.density_data_loaded:
            self.load_densities(trajectory_name)

        charge_state_data_dict: Dict[str, ChargeStateData] = processor.get_electrons_around_proton(data_handler.all_data, parameters)
        visualiser.plot_charge_state(trajectory_name, charge_state_data_dict, parameters, data_handler.all_data)
