from DataLoader import DataLoader
from DataProcessor import DataProcessor
from DataVisualiser import DataVisualiser

import os

class DataHandler:
    def __init__(self, directory):
        self.directory_name = directory
        self.trajectory_name = os.path.basename(directory.rstrip("/"))
        self.data_loader = DataLoader(directory)
        self.data_processor = DataProcessor()
        self.data_visualiser = DataVisualiser()

        self.atoms_dict = {}
        self.calc_dict = {}
        self.fits = {}


