import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from DataLoader import Data
from typing import Dict


class DataVisualiser:
    def __init__(self):
        self.geant4_stopping_data = {
            "energies": np.arange(5, 251, 5),
            "stopping powers": np.array([8.41642534,  8.39403791,  9.61722966, 10.47531744, 11.12587244, 11.66446394, 12.3387949 , 12.60404201, 12.70727751, 12.74141819, 12.84432534, 12.96025452, 12.82696325, 12.74379889, 12.69518111, 12.56137771, 12.41201793, 12.33257865, 12.21160185, 12.09496251, 12.01462424, 11.898194, 11.77969121, 11.66621298, 11.50000087, 11.37880456, 11.31399968, 11.22186803, 11.06282411, 10.93625158, 10.81506187, 10.75385659, 10.58539531, 10.49800502, 10.48398759, 10.4024316, 10.27643047, 10.19540826, 10.01933624, 9.98213156, 9.99531111, 9.7817623 , 9.66781196,  9.78197921, 9.71504747, 9.51892708,  9.52356048,  9.36827356,  9.36659151,  9.22948771])
        }

        path_to_srim_data = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/Hydrogen_in_Aluminium_SRIM.txt"
        self.srim_stopping_data = self.load_srim(path_to_srim_data)

    ########################################
    # FOR PLOTTING STOPPING POWER ANALYSIS #
    ########################################

    def load_srim(self, path):
        raw_data = pd.read_csv(path,
                               sep=r"\s+",
                               skiprows = 23,
                               header = None,
                               skipfooter=14,
                               engine="python")

        data = pd.DataFrame(data = {
            "E_k [keV]" : raw_data[0],
            "S_e [eV/A]" : raw_data[2],
            "S_n [eV/A]" : raw_data[3],
            "S [eV/A]" : raw_data[2] + raw_data[3]
        })

        data.loc[0, "E_k [keV]"] /= 1e3   # just this first one is in eV??

        return data

    def plot_all_fits(self, all_data: Dict[str, Data], fits):

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
        fig.suptitle("Stopping Power Things")
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


            ######################
            # PLOT GEANT4 DATA?? #
            ######################
            plot_geant4_stopping = True
            if plot_geant4_stopping:
                try:
                    index = np.where(self.geant4_stopping_data["energies"] == int(energy.rstrip(" keV")))[0]
                    pfit = np.poly1d([-1e-3*self.geant4_stopping_data["stopping powers"][index][0], self.geant4_stopping_data["energies"][index][0]])
                    axs[i, 0].plot(distance_travelled, pfit(distance_travelled), label=rf"Geant4: $S_e$ = {self.geant4_stopping_data['stopping powers'][index][0]:.1f} [eV/$\AA$]")
                except:
                    print(f"no geant4 fit for {energy}")


            # TODO: WONT FIND ANYTHING AT ALL
            plot_srim_stopping = False
            if plot_srim_stopping:
                try:
                    index = self.srim_stopping_data[np.isclose(self.srim_stopping_data["E_k [keV]"], energy.rstrip(" keV"))]
                    p = np.poly1d(-1e-3 * self.srim_stopping_data.loc[index, "S [ev/A]"], self.srim_stopping_data.loc[index, "E_k [keV]"])
                    axs[i, 0].plot(distance_travelled, p(distance_travelled)) #, label=rf"SRIM: S = {self.}"
                except:
                    print(f"no srim fit for {energy}")

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

                axs[i, 1].plot([distance_travelled[0], distance_travelled[-1]], [-fit[0], -fit[0]], color="red")
                axs[i, 1].fill_between(distance_travelled, np.ones(len(distance_travelled))*-fit[0] - stopping_power_uncertainty,
                                                        np.ones(len(distance_travelled))*-fit[0] + stopping_power_uncertainty,
                                                        color="red", alpha=0.25)

                # plot points used in the fitting
                n_timesteps = len(distance_travelled)
                axs[i, 0].plot(distance_travelled[crop[0]:n_timesteps - (crop[1] or 0)], kinetic_energies[crop[0]:n_timesteps - (crop[1] or 0)], "x", color="red")
                axs[i, 1].plot(distance_travelled[crop[0] : n_timesteps - (crop[1] or 0) - 1], -np.diff(kinetic_energies)[crop[0] : n_timesteps - (crop[1] or 0) - 1]/np.diff(distance_travelled)[crop[0] : n_timesteps - (crop[1] or 0) - 1], "x", color="red")
                axs[i,0].legend()

        plt.show()

    def geant4_comparison(self, stopping_power_data):
        """
        plots comparison to Monte Carlo stopping curve

        Parameters:
        stopping_power_data (Dict[str, np.ndarray()])
        """
        fig,ax = plt.subplots(figsize=(15,5))
        ax.set_xlabel("projectile initial kinetic energy [keV]")
        ax.set_ylabel(r"stopping power [eV/$\AA$]")

        # ax.plot(self.geant4_stopping_data["energies"], self.geant4_stopping_data["stopping powers"], "-", label="GEANT4")
        ax.plot(self.srim_stopping_data["E_k [keV]"], self.srim_stopping_data["S [eV/A]"], label="SRIM")
        for trajectory_name, trajectory in stopping_power_data.items():
            ax.plot(trajectory["energies"], trajectory["stopping powers"], "-x", label=trajectory_name)
        ax.legend()
        plt.show()


    ######################################
    # FOR VISUALISING SIMULATION RESULTS #
    ######################################


    def visualise_electron_density(self, electron_densities):
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
                return  # Don't allow inverted color limits
            im_density.set_clim(vmin, vmax)
            cb_density.update_normal(im_density)
            fig.canvas.draw_idle()

        # Connect sliders
        slider_t.on_changed(update)
        slider_slice.on_changed(update)
        slider_vmin.on_changed(update_vlim)
        slider_vmax.on_changed(update_vlim)

        plt.show()




    # def visualise_electron_density(self, electron_density_list):
    #     """
    #     creates animation using matplotlib sliders to display 2d slice of electron density
    #     has sliders to control slice position, timestep, vmin&vmax
    #     can also toggle logarithmic plotting
    #     """
    #
    #     shape = np.shape(electron_density_list)
    #     timesteps = shape[0]
    #
    #     # Initial values
    #     init_timestep = 0
    #     init_slice = shape[1] // 2
    #     use_log = False
    #
    #     #########################
    #     # CALCULATE VMIN/VMAX'S #
    #     #########################
    #
    #     all_data = np.stack(electron_density_list)
    #     vmin_raw = all_data.min()
    #     vmax_raw = all_data.max()
    #     log_data = np.log10(np.clip(all_data, 1e-10, None))
    #     vmin_log = log_data.min()
    #     vmax_log = log_data.max()
    #
    #     ##############
    #     # SETUP PLOT #
    #     ##############
    #
    #     fig, ax = plt.subplots(figsize=(8, 6))
    #     plt.subplots_adjust(left=0.25, bottom=0.45)  # Extra room for sliders
    #
    #     def get_plot_data(t, slice, log=False):
    #         data = electron_density_list[t][:, slice, :]
    #         data = np.rot90(data)
    #         if log:
    #             data = np.log10(np.clip(data, 1e-10, None))
    #         return data
    #
    #     # Initial image
    #     plot_data = get_plot_data(init_timestep, init_slice, use_log)
    #     im = ax.imshow(plot_data, aspect='auto', vmin=vmin_raw, vmax=vmax_raw)
    #     cb = fig.colorbar(im, ax=ax)
    #     title = ax.set_title(f"Timestep: {init_timestep}, Slice: {init_slice}, Log: {use_log}")
    #
    #     # Axes for sliders and checkbox
    #     ax_timestep = plt.axes((0.25, 0.35, 0.65, 0.03))
    #     ax_slice = plt.axes((0.25, 0.3, 0.65, 0.03))
    #     ax_vmin = plt.axes((0.25, 0.2, 0.65, 0.03))
    #     ax_vmax = plt.axes((0.25, 0.15, 0.65, 0.03))
    #     ax_check = plt.axes((0.05, 0.5, 0.15, 0.1))
    #
    #     # Create sliders and checkbox
    #     slider_timestep = Slider(ax_timestep, 'Timestep', 0, len(electron_density_list) - 1, valinit=init_timestep, valstep=1)
    #     slider_slice = Slider(ax_slice, 'Slice', 0, shape[1] - 1, valinit=init_slice, valstep=1)
    #     slider_vmin = Slider(ax_vmin, 'vmin', vmin_raw, vmax_raw, valinit=vmin_raw)
    #     slider_vmax = Slider(ax_vmax, 'vmax', vmin_raw, vmax_raw, valinit=vmax_raw)
    #     check = CheckButtons(ax_check, ['Log scale'], [use_log])
    #
    #     # update called when a slider is moved
    #     def update(val):
    #         t = int(slider_timestep.val)
    #         s = int(slider_slice.val)
    #         log_flag = check.get_status()[0]
    #
    #         # Determine colour range
    #         if log_flag:
    #             current_vmin = slider_vmin.val
    #             current_vmax = slider_vmax.val
    #         else:
    #             current_vmin = slider_vmin.val
    #             current_vmax = slider_vmax.val
    #
    #         # Update data
    #         new_data = get_plot_data(t, s, log=log_flag)
    #         im.set_data(new_data)
    #         im.set_clim(vmin=current_vmin, vmax=current_vmax)
    #         title.set_text(f"Timestep: {t}, Slice: {s}, Log: {log_flag}")
    #         fig.canvas.draw_idle()
    #
    #     # Log scale toggle — resets vmin/vmax sliders
    #     def toggle_log(label):
    #         log_flag = check.get_status()[0]
    #         if log_flag:
    #             slider_vmin.valmin = vmin_log
    #             slider_vmin.valmax = vmax_log
    #             slider_vmax.valmin = vmin_log
    #             slider_vmax.valmax = vmax_log
    #             slider_vmin.set_val(vmin_log)
    #             slider_vmax.set_val(vmax_log)
    #         else:
    #             slider_vmin.valmin = vmin_raw
    #             slider_vmin.valmax = vmax_raw
    #             slider_vmax.valmin = vmin_raw
    #             slider_vmax.valmax = vmax_raw
    #             slider_vmin.set_val(vmin_raw)
    #             slider_vmax.set_val(vmax_raw)
    #         slider_vmin.ax.set_xlim(slider_vmin.valmin, slider_vmin.valmax)
    #         slider_vmax.ax.set_xlim(slider_vmax.valmin, slider_vmax.valmax)
    #         update(None)
    #
    #     # runs update functions on slider movement
    #     slider_timestep.on_changed(update)
    #     slider_slice.on_changed(update)
    #     slider_vmin.on_changed(update)
    #     slider_vmax.on_changed(update)
    #     check.on_clicked(toggle_log)
    #
    #     plt.show()

