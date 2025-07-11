import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

import utils

# imports for typing
from typing import Dict, List
from DataLoader import Data
from DataProcessor import Fit, ChargeStateData
from numpy.typing import NDArray

import logging

logger = logging.getLogger(__name__)


class DataVisualiser:
    def __init__(self):
        self.geant4_stopping_data = {
            "energies": np.arange(5, 251, 5),
            "stopping powers": np.array([8.41642534,  8.39403791,  9.61722966, 10.47531744, 11.12587244, 11.66446394, 12.3387949 , 12.60404201, 12.70727751, 12.74141819, 12.84432534, 12.96025452, 12.82696325, 12.74379889, 12.69518111, 12.56137771, 12.41201793, 12.33257865, 12.21160185, 12.09496251, 12.01462424, 11.898194, 11.77969121, 11.66621298, 11.50000087, 11.37880456, 11.31399968, 11.22186803, 11.06282411, 10.93625158, 10.81506187, 10.75385659, 10.58539531, 10.49800502, 10.48398759, 10.4024316, 10.27643047, 10.19540826, 10.01933624, 9.98213156, 9.99531111, 9.7817623 , 9.66781196,  9.78197921, 9.71504747, 9.51892708,  9.52356048,  9.36827356,  9.36659151,  9.22948771])
        }

        path_to_srim_data = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/Hydrogen_in_Aluminium_SRIM.txt"
        self.srim_stopping_data = self.load_srim(path_to_srim_data)
        path_to_srim_experimental_data = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/H_in_Aluminium_experimental.csv"
        self.srim_experimental_data = self.load_srim_experimental(path_to_srim_experimental_data)

    ########################################
    # FOR PLOTTING STOPPING POWER ANALYSIS #
    ########################################

    def load_srim(self, path: str) -> pd.DataFrame:
        raw_data = pd.read_csv(path,
                               sep=r"\s+",
                               skiprows = 23,
                               header = None,
                               skipfooter=14,
                               engine="python")

        for i, unit in enumerate(raw_data[1].to_numpy()):
            multiply_by = 1.0  # for converting to keV
            if unit == "eV":
                multiply_by = 1.0e-3
            elif unit == "keV":
                pass
            elif unit == "MeV":
                multiply_by = 1.0e3
            else:
                raise ValueError("Error reading units in SRIM file")

            raw_data.loc[i, 0] *= multiply_by

        data_df = pd.DataFrame(data = {
            "E_k [keV]" : raw_data[0],
            "S_e [eV/A]" : raw_data[2],
            "S_n [eV/A]" : raw_data[3],
            "S [eV/A]" : raw_data[2] + raw_data[3]
        })

        # rounds energy column to 3sf so that it is much easier to search through
        data_df["E_k [keV]"] = np.array([utils.round_sf(x) for x in data_df["E_k [keV]"].to_numpy()])

        return data_df

    def load_srim_experimental(self, path: str) -> pd.DataFrame:
        """ https://www-nds.iaea.org/stopping/tuples/H-Al """

        E15eVcm2atom_to_evA_conversion = 1.6582 # this number is from the bottom of the SRIM data file (not experimental data file)

        raw_data = pd.read_csv(path)
        # some of them are deuterium and tritium
        data = pd.DataFrame(data = {
            "E_k [keV]" : raw_data[raw_data["ion_isotope"] == 1.0]["energy"],
            "S [eV/A]" : raw_data[raw_data["ion_isotope"] == 1.0]["stopping_power"] / E15eVcm2atom_to_evA_conversion,
        })

        return data

    def plot_all_fits(self, trajectory_name: str, all_data: Dict[str, Data], fits):

        """
        for every projectile energy, creates 2 subplots
        1. projectile kinetic energy against distance travelled
        2. gradient of kinetic energy against distance travelled
        """

        projectile_positions = {key: value.projectile_positions for key,value in all_data.items()}
        projectile_kinetic_energies = {key: value.projectile_kinetic_energies for key,value in all_data.items()}

        #################
        # CREATE FIGURE #
        #################

        n_subplots = len(projectile_kinetic_energies.keys())
        fig,axs = plt.subplots(n_subplots, 2, figsize=(15, 5*n_subplots), sharex="col")
        if n_subplots == 1:
            axs = np.array([axs])
        fig.suptitle(f"Stopping Power Fits for {trajectory_name}")
        _ = [ax.set_ylabel("Kinetic Energy [keV]") for ax in axs[:,0]]
        _ = [ax.set_ylabel(r"$\frac{dKE}{dx}$ [keV/$\AA$]") for ax in axs[:,1]]
        _ = [ax.set_xlabel(r"Projectile Position [$\AA$]") for ax in axs[n_subplots-1]]

        #######################################
        # ITERATE OVER EACH PROJECTILE ENERGY #
        #######################################

        for i in range(n_subplots):

            # get raw data
            energy, kinetic_energies = list(projectile_kinetic_energies.items())[i]
            kinetic_energies = np.array(kinetic_energies) * 1e-3  # convert from eV to keV
            _, positions = list(projectile_positions.items())[i]
            initial_position = positions[0]
            distance_travelled = np.array([np.linalg.norm(position - initial_position) for position in positions])

            # plot kinetic energy on left plot
            axs[i, 0].plot(distance_travelled, kinetic_energies, "x")

            # plot gradient of kinetic energy on right plot
            axs[i, 1].plot(distance_travelled[:-1], -np.diff(kinetic_energies)/np.diff(distance_travelled))

            plot_geant4_stopping = False
            if plot_geant4_stopping:
                try:
                    index = np.where(self.geant4_stopping_data["energies"] == int(energy.rstrip(" keV")))[0]
                    pfit = np.poly1d([-1e-3*self.geant4_stopping_data["stopping powers"][index][0], self.geant4_stopping_data["energies"][index][0]])
                    axs[i, 0].plot(distance_travelled, pfit(distance_travelled), label=rf"Geant4: $S_e$ = {self.geant4_stopping_data['stopping powers'][index][0]:.1f} [eV/$\AA$]")
                except Exception as e:
                    logger.info(f"no geant4 fit for {energy}, error: {e}")

            plot_srim_stopping = True
            if plot_srim_stopping:
                try:
                    energy_val = utils.round_sf(float(energy.rstrip(" keV")))
                    row = self.srim_stopping_data.loc[self.srim_stopping_data["E_k [keV]"] == energy_val]
                    stopping_power_SRIM = row["S [eV/A]"].values[0]
                    energy_SRIM = row["E_k [keV]"].values[0]

                    p = np.poly1d([-1e-3 * stopping_power_SRIM, energy_SRIM])
                    axs[i, 0].plot(distance_travelled, p(distance_travelled), label=rf"SRIM: S = {stopping_power_SRIM:.1f} [eV/$\AA$]")
                    axs[i, 1].axhline(1e-3 * stopping_power_SRIM, color="C1", label=rf"SRIM: S = {stopping_power_SRIM:.1f} [eV/$\AA$]")

                except Exception as e:
                    logging.warning(f"no srim fit for {energy}, error: {e}")

            ########################
            # PLOT STOPPING FITS?? #
            ########################
            if fits:
                # extract fit information from dictionary of Fit objects
                fit = fits[energy].fit
                pfit = np.poly1d(fit)
                cov = fits[energy].cov
                crop = fits[energy].crop

                # calculate stopping powers from fits to display on the plots
                stopping_power = -fit[0]
                stopping_power_uncertainty = np.sqrt(cov[0][0])
                label = rf"$S_e$ = {(1e3*stopping_power):.1f} $\pm$ {(1e3*stopping_power_uncertainty):.1f} [eV/$\AA$]"


                # plot fits and uncertainties
                axs[i, 0].plot(distance_travelled, pfit(distance_travelled), color="red", label=label)

                axs[i, 1].plot([distance_travelled[0], distance_travelled[-1]], [-fit[0], -fit[0]], color="red", label=label)
                axs[i, 1].fill_between(distance_travelled, np.ones(len(distance_travelled))*-fit[0] - stopping_power_uncertainty,
                                                        np.ones(len(distance_travelled))*-fit[0] + stopping_power_uncertainty,
                                                        color="red", alpha=0.25)

                # plot points used in the fitting
                n_timesteps = len(distance_travelled)
                axs[i, 0].plot(distance_travelled[crop[0]:n_timesteps - (crop[1] or 0)], kinetic_energies[crop[0]:n_timesteps - (crop[1] or 0)], "x", color="red")
                axs[i, 1].plot(distance_travelled[crop[0] : n_timesteps - (crop[1] or 0) - 1], -np.diff(kinetic_energies)[crop[0] : n_timesteps - (crop[1] or 0) - 1]/np.diff(distance_travelled)[crop[0] : n_timesteps - (crop[1] or 0) - 1], "x", color="red")
                axs[i, 0].legend()
                axs[i, 1].legend()

        plt.show()

    def compare_fits(self, trajectory_names: List[str], energy: str, comparison_data: Dict[str, Data]):
        fig,axs = plt.subplots(1, 2, figsize=(12, 7))
        axs[0].set_title("projectile KE")
        axs[1].set_title("instantaneous stopping power")

        _ = [ax.set_xlabel("distance travelled [Angstrom]") for ax in axs]
        axs[0].set_ylabel("KE [eV]")
        axs[1].set_ylabel("S [eV/$\AA$]")

        for trajectory_name in trajectory_names:
            projectile_positions = comparison_data[trajectory_name].projectile_positions
            projectile_distances = np.array([np.linalg.norm(position - projectile_positions[0]) for position in projectile_positions])
            projectile_kinetic_energies = comparison_data[trajectory_name].projectile_kinetic_energies

            axs[0].plot(projectile_distances, projectile_kinetic_energies, label=trajectory_name)
            axs[1].plot(projectile_distances[:-1], -np.diff(projectile_kinetic_energies)/np.diff(projectile_distances), label=trajectory_name)


        plot_srim_stopping = True
        if plot_srim_stopping:
            try:
                energy_val = utils.round_sf(float(energy.rstrip(" keV")))
                row = self.srim_stopping_data.loc[self.srim_stopping_data["E_k [keV]"] == energy_val]
                stopping_power_SRIM = row["S [eV/A]"].values[0]
                energy_SRIM = row["E_k [keV]"].values[0]

                p = np.poly1d([-1 * stopping_power_SRIM, 1e3*energy_SRIM])
                axs[0].plot(projectile_distances, p(projectile_distances), "--", color="red", alpha=0.75, label=rf"SRIM: S = {stopping_power_SRIM:.1f} [eV/$\AA$]")
                axs[1].axhline(stopping_power_SRIM, ls="--", color="red", alpha=0.75, label=rf"SRIM: S = {stopping_power_SRIM:.1f} [eV/$\AA$]")

            except Exception as e:
                logging.warning(f"no srim fit for {energy}, error: {e}")

        _ = [ax.legend() for ax in axs]
        plt.show()

    def montecarlo_comparison(self, all_fit_info: Dict[str, Dict[str, Fit]]):

        """
        plots comparison to Monte Carlo stopping curve

        Parameters:
        stopping_power_data (Dict[str, Dict[str, Fit]]) where keys are trajectory names, values are a dictionary with energies as keys, Fit instances as values
        """
        fig,ax = plt.subplots(figsize=(15,5))
        ax.set_xlabel("projectile initial kinetic energy [keV]")
        ax.set_ylabel(r"stopping power [eV/$\AA$]")

        # ax.plot(self.geant4_stopping_data["energies"], self.geant4_stopping_data["stopping powers"], "-", label="GEANT4")
        ax.plot(self.srim_stopping_data["E_k [keV]"], self.srim_stopping_data["S [eV/A]"], label="SRIM")
        ax.plot(self.srim_experimental_data["E_k [keV]"], self.srim_experimental_data["S [eV/A]"], ".", color="lightgray", label="experimental data")

        for trajectory_name, fit_info in all_fit_info.items():
            energies = [int(energy_string.rstrip(" keV")) for energy_string in list(fit_info.keys())]
            # stopping_powers = -1e3 * fit_info.fit[0]   # want to plot in eV/Ang
            stopping_powers = [-1e3 * fit.fit[0] for fit in fit_info.values()]
            ax.plot(energies, stopping_powers, "-x", label=trajectory_name)

        ax.set_xlim((-500, 2000))

        # TODO: add option to log scale plot because a lot of papers seem to plot things this way
        # ax.set_xscale("log")
        # ax.set_xlim((0.5, 500))
        ax.legend()
        plt.show()


    ######################################
    # FOR VISUALISING SIMULATION RESULTS #
    ######################################


    def visualise_electron_density(self, electron_densities: NDArray):
        """uses matplotlib widget sliders to create interactive plot to visualise electron density"""
        max_t = np.shape(electron_densities)[0] - 1
        max_slice = np.shape(electron_densities)[2] - 1

        # Initial values
        init_t = 6
        init_slice = 48

        # Create two vertically stacked subplots
        fig, (ax_density, ax_change) = plt.subplots(2, 1, figsize=(12, 10))
        plt.subplots_adjust(left=0.1, right=0.85, bottom=0.25, hspace=0.3)

        # Initial image data
        electron_image = np.rot90(electron_densities[init_t][:, init_slice, :])
        log_electron_image = np.log10(np.clip(electron_image, 1e-10, None))
        electron_image_initial = np.rot90(electron_densities[0][:, init_slice, :])
        electron_density_change_image = electron_image - electron_image_initial

        # Color scale limits
        vlim_change = np.max(np.abs(electron_density_change_image))
        log_vmin_density = np.min(log_electron_image)
        # log_vmin_density = -2
        log_vmax_density = np.max(log_electron_image)

        # Electron density plot
        im_density = ax_density.imshow(log_electron_image, cmap="viridis", vmin=log_vmin_density, vmax=log_vmax_density)
        cb_density = plt.colorbar(im_density, ax=ax_density)
        ax_density.set_title("(log-scale) Electron Density")

        # Change in electron density plot
        im_change = ax_change.imshow(electron_density_change_image, cmap="bwr", vmin=-vlim_change, vmax=vlim_change)
        cb_change = plt.colorbar(im_change, ax=ax_change)
        ax_change.set_title("Change in Electron Density")

        # Slider axes (horizontal, bottom of figure)
        ax_t = plt.axes((0.25, 0.15, 0.5, 0.03))
        ax_slice = plt.axes((0.25, 0.1, 0.5, 0.03))

        # Create sliders
        slider_t = Slider(ax_t, 't (time)', 0, max_t, valinit=init_t, valstep=1)
        slider_slice = Slider(ax_slice, 'slice', 0, max_slice, valinit=init_slice, valstep=1)

        # Colorbar slider axes (vertical, right of top colorbar)
        ax_vmin = plt.axes((0.88, 0.6, 0.015, 0.25))
        ax_vmax = plt.axes((0.91, 0.6, 0.015, 0.25))

        # Create color limit sliders
        slider_vmin = Slider(ax_vmin, 'vmin', -2, log_vmax_density, valinit=log_vmin_density, orientation='vertical')
        slider_vmax = Slider(ax_vmax, 'vmax', -2, log_vmax_density, valinit=log_vmax_density, orientation='vertical')

        # Update function for t/slice sliders
        def update(val):
            t_val = int(slider_t.val)
            slice_val = int(slider_slice.val)

            electron_image = np.rot90(electron_densities[t_val][:, slice_val, :])
            log_electron_image = np.log10(np.clip(electron_image, 1e-10, None)) 
            electron_image_initial = np.rot90(electron_densities[0][:, slice_val, :])
            electron_density_change_image = electron_image - electron_image_initial

            im_density.set_data(log_electron_image)
            im_density.set_clim(slider_vmin.val, slider_vmax.val)

            vlim_change = np.max(np.abs(electron_density_change_image))
            im_change.set_data(electron_density_change_image)
            im_change.set_clim(-vlim_change, vlim_change)

            cb_density.update_normal(im_density)
            cb_change.update_normal(im_change)
            fig.canvas.draw_idle()

        # Update function for vmin/vmax sliders
        def update_vlim(val):
            vmin = slider_vmin.val
            vmax = slider_vmax.val
            if vmin >= vmax:
                return
            im_density.set_clim(vmin, vmax)
            cb_density.update_normal(im_density)
            fig.canvas.draw_idle()

        slider_t.on_changed(update)
        slider_slice.on_changed(update)
        slider_vmin.on_changed(update_vlim)
        slider_vmax.on_changed(update_vlim)

        plt.show()





    def plot_charge_state(self,
                          trajectory_name: str,
                          charge_state_data_dict: Dict[str, List[ChargeStateData]],
                          parameters: Dict[str, int],
                          all_data: Dict[str, Data]):
        """
        plots charge state around the proton for each of the energies simulated for a given trajectory
        Parameters:
            trajectory_name (str)
            charge_state_data_dict (Dict[str, List[ChargeStateData]])
            parameters (Dict[str, int])
            all_data (Dict[str, Data])
        """

        first_data = list(all_data.values())[0]
        pixel_to_angstrom = first_data.cell[0] / np.shape(first_data.electron_densities)[1]

        fig,axs = plt.subplots(len(charge_state_data_dict.keys()), figsize=(16, 8), sharex="col")
        if len(charge_state_data_dict.keys()) == 1:
            axs = np.array([axs])
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle(f"""Charge state for {trajectory_name} trajectory
                         (integrated over {parameters["size"] * pixel_to_angstrom:.1f} Angstrom)""")
        for i, key in enumerate(charge_state_data_dict.keys()):


            crop_low, crop_high = charge_state_data_dict[key].crop

            positions = all_data[key].projectile_positions
            initial_position = positions[0]
            distance_travelled = np.array([np.linalg.norm(position - initial_position) for position in positions])

            # GOTTEN FROM A DFT SIMULATION OF A SINGLE HYDROGEN ATOM
            # IDK WHERE TO PUT THEM...
            DFT_HYDROGEN_CHARGE_STATE = 875   # 99% of electron density
            DFT_HYDROGEN_CHARGE_STATE_RADIUS = 2.2   # radius [Angstroms] which contains 99% of hydrogen's electron density

            charge_state = charge_state_data_dict[key].density_around_projectile / DFT_HYDROGEN_CHARGE_STATE

            axs[i].plot(distance_travelled, charge_state, label=key)
            # axs[i].axhline(np.mean(charge_state), linestyle = "--", color="red", label=f"average charge state = {np.mean(charge_state):.2f}")
            axs[i].set_title(key)
            axs[i].legend()
        plt.show()





