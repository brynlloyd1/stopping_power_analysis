from DataHandler import DataHandler
from typing import List

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

    def set_name(self, old_name, new_name):
        """changes the name of a trajectory, which is stored by the DataHandler instance"""
        if new_name in self.data_handlers.keys():
            raise KeyError(f"{new_name} already exists")

        self.data_handlers[old_name].trajectory_name = new_name
        self.data_handlers[new_name] = self.data_handlers[old_name].pop()

    def load_data(self, trajectory_name):
        """instructs the instance of DataLoaded stored in the DataHandler instance to load in the data"""
        data_handler = self.data_handlers[trajectory_name]
        # data_handler.atoms_dict, data_handler.calc_dict = data_handler.data_loader.load_data()
        data_handler.all_data = data_handler.data_loader.load_data()

    def calculate_stopping_curve(self,
                                 trajectory_name: str,
                                 crop: List[int | None] = [11, None]):
        """tells the DataProcessor instance to perform fits to kinetic energy data"""
        handler = self.data_handlers[trajectory_name]
        processor = handler.data_processor
        # fits_information = processor.calculate_stopping_powers(handler.atoms_dict, crop=crop)
        fits_information = processor.calculate_stopping_powers(handler.all_data, crop=crop)
        handler.fits = fits_information


    def view_fits(self, trajectory_name: str):
        """plots data and fits which comparison to Monte Carlo stopping powers"""
        handler = self.data_handlers[trajectory_name]
        visualiser = handler.data_visualiser
        visualiser.plot_all_fits(handler.all_data, handler.fits)
        # visualiser.plot_all_fits(handler.atoms_dict, handler.fits)


    def compare_to_geant4(self, trajectory_names):
        """comparison of a stopping curve to Monte Carlo stopping power curve"""
        stopping_power_data = {}
        for trajectory_name in trajectory_names:
            handler = self.data_handlers[trajectory_name]
            energies = []
            stopping_powers = []

            for energy, fit_info in handler.fits.items():
                energies.append(int(energy.rstrip(" keV")))
                stopping_powers.append(-fit_info.fit[0]*1e3)   # fit[0] is in keV/Angstrom. convert to eV/Ang

            stopping_power_data[trajectory_name] = {"energies": energies,
                                                    "stopping powers": stopping_powers}
        self.data_handlers[trajectory_names[0]].data_visualiser.geant4_comparison(stopping_power_data)


    def visualise_electron_density(self, trajectory_name, energy):
        """DataVisualiser.visualise_electron_density called"""
        handler = self.data_handlers[trajectory_name]
        loader = handler.data_loader
        visualiser = handler.data_visualiser

        if loader.check_for_npy(handler.directory_name, energy):
            electron_density_list = loader.load_from_npy(handler.directory_name, energy)
        else:
            electron_density_list = [calc.get_all_electron_density() for calc in handler.calc_dict[energy]]
            loader.save_to_npy(handler.directory_name, energy, electron_density_list)

        visualiser.visualise_electron_density(electron_density_list)
