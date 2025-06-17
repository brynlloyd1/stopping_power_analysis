from DataLoader import DataLoader, Data
from DataProcessor import DataProcessor, Fit
from DataVisualiser import DataVisualiser

import os
from typing import Dict

class DataHandler:
    def __init__(self, directory):
        self.directory_name = directory
        self.trajectory_name = os.path.basename(directory.rstrip("/"))
        self.data_loader = DataLoader(directory)
        self.data_processor = DataProcessor()
        self.data_visualiser = DataVisualiser()

        self.all_data: Dict[str, Data] = {}
        self.fits: Dict[str, Fit] = {}


