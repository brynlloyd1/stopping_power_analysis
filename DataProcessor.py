import numpy as np

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


    def calculate_stopping_powers(self, atoms_dict, crop=[None, None]):
        """
        Parameters:
        atoms_dict (Dict[str, List[Atoms]])
        calc_dict (Dict[str, List[GPAW]])
        """

        projectile_positions = {energy: [atoms.get_positions()
                                            for atoms in atoms_list]
                                for energy, atoms_list in atoms_dict.items()}

        projectile_kinetic_energies = {energy: [atoms.get_kinetic_energy()
                                                for atoms in atoms_list]
                                        for energy, atoms_list in atoms_dict.items()}

        fits_information = {}   # dictionary of energy : Fit

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
