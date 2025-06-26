from EhrenfestAnalysis import EhrenfestAnalysis

from logger_config import setup_logging
import logging

setup_logging()
logger = logging.getLogger(__name__)
logger.info("Starting Ehrenfest Analysis")

# parameters for the charge state analysis
parameters = {
    "size" : 5,
    "offset" : 0,
    "crop_low" : 4,
    "crop_high" : 64,
}

analysis = EhrenfestAnalysis()


###################   hydrogen    ########################
# directory_322_hyperchannelling_hydrogen = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_hyperchannelling_hydrogen/"
# name_322_hyperchannelling_hydrogen = analysis.initialise_analysis(directory_322_hyperchannelling_hydrogen)


directory_622_hyperchannelling_hydrogen = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_hydrogen/"
name_622_hyperchannelling_hydrogen = analysis.initialise_analysis(directory_622_hyperchannelling_hydrogen)
analysis.calculate_stopping_curve(name_622_hyperchannelling_hydrogen)
analysis.view_fits(name_622_hyperchannelling_hydrogen)
# analysis.set_which_energies(name_622_hyperchannelling_hydrogen, which_energies = ["40 keV", "60 keV"])
# analysis.analyse_proton_charge_state(name_622_hyperchannelling_hydrogen, parameters)

# directory_322_presampled1_hydrogen = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_presampled1_hydrogen/"
# name_322_presampled1_hydrogen = analysis.initialise_analysis(directory_322_presampled1_hydrogen)


###################   proton but its still actually hydrogen    ########################
# directory_622_hyperchannelling_proton = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton/"
# name_622_hyperchannelling_proton = analysis.initialise_analysis(directory_622_hyperchannelling_proton)
# analysis.set_which_energies(name_622_hyperchannelling_proton, which_energies=["40 keV"])



# analysis.analyse_proton_charge_state(name_622_hyperchannelling_proton, parameters)
# analysis.calculate_stopping_curve(name_622_hyperchannelling_proton)
# analysis.view_fits(name_622_hyperchannelling_proton)

# directory_333_presampled1_proton = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton/"
# name_333_presampled1_proton = analysis.initialise_analysis(directory_333_presampled1_proton)
# analysis.analyse_proton_charge_state(name_333_presampled1_proton)

###################   proton that might actually be hydrogen    ########################
# directory_622_hyperchannelling_proton2 = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton2/"
# name_622_hyperchannelling_proton2 = analysis.initialise_analysis(directory_622_hyperchannelling_proton2)
# analysis.analyse_proton_charge_state(name_622_hyperchannelling_proton2, parameters)
# analysis.calculate_stopping_curve(name_622_hyperchannelling_proton2)
# # analysis.view_fits(name_622_hyperchannelling_proton2)
# analysis.visualise_electron_density(name_622_hyperchannelling_proton2, energy = "40 keV")


###################   proton that might actually be hydrogen    ########################
directory_422_hyperchannelling_proton2_surface = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton2_surface/"
name_422_hyperchannelling_proton2_surface = analysis.initialise_analysis(directory_422_hyperchannelling_proton2_surface)
analysis.calculate_stopping_curve(name_422_hyperchannelling_proton2_surface, crop = [15, 45])
analysis.view_fits(name_422_hyperchannelling_proton2_surface)
# analysis.analyse_proton_charge_state(name_422_hyperchannelling_proton2_surface, parameters)
# analysis.visualise_electron_density(name_422_hyperchannelling_proton2_surface, energy = "40 keV")
# analysis.visualise_electron_density(name_422_hyperchannelling_proton2_surface, energy = "60 keV")



directory_422_hyperchannelling_proton3_surface = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton3_surface/"
# name_422_hyperchannelling_proton3_surface = analysis.initialise_analysis(directory_422_hyperchannelling_proton3_surface)
# analysis.calculate_stopping_curve(name_422_hyperchannelling_proton3_surface, crop = [15, 49])
# analysis.view_fits(name_422_hyperchannelling_proton3_surface)
# analysis.analyse_proton_charge_state(name_422_hyperchannelling_proton3_surface, parameters)
# analysis.visualise_electron_density(name_422_hyperchannelling_proton3_surface, energy = "1000 keV")












# =====================================================================================================================================
#
# def find_projectile(projectile_positions, electron_density, cell):
#     projectile_position_indices = np.array([(projectile_position % cell) / cell * np.shape(electron_density) for projectile_position in projectile_positions], dtype="int64")
#
#     return projectile_position_indices
#
#
#
# analysis.load_densities(name_622_hyperchannelling_proton)
# energy = "40 keV"
# t = 4
# slice_pos = 40
#
# energy_data = analysis.data_handlers[name_622_hyperchannelling_proton].all_data[energy]
# projectile_positions = energy_data.projectile_positions
# electron_densities = energy_data.electron_densities
# cell = energy_data.cell
#
# proj_indices = find_projectile(projectile_positions, electron_densities[t], cell)
#
# size = 10
# offset = -3
#
#
# density_around_projectile = []
# for t in range(len(electron_densities)):
#     proj_index = proj_indices[t]
#     if proj_index[0] < size or proj_index[0] + size > np.shape(electron_densities)[1]:
#         density_around_projectile.append(None)
#     electron_density = electron_densities[t]
#     density_around_projectile.append(np.sum(electron_density[proj_index[0] - size + offset : proj_index[0] + size + offset,
#                             proj_index[1] - size : proj_index[1] + size,
#                             proj_index[2] - size : proj_index[2] + size]))
#
# plt.plot(density_around_projectile)
# plt.show()
#
#
#
# n_timesteps = len(electron_densities)
# shape = np.shape(electron_densities[0])
# init_t = 0
# init_slice = 25
#
# fig, ax = plt.subplots(figsize=(10, 5))
# plt.subplots_adjust(left=0.1, bottom=0.25)  # Make space for sliders
#
# # Plot initial image
# data = np.rot90(np.log10(electron_densities[init_t][:, init_slice, :]))
# im = ax.imshow(data, cmap='viridis')
# cb = plt.colorbar(im, ax=ax)
#
# # Draw initial circle
# center = (proj_indices[init_t][0] + offset, proj_indices[init_t][1])
# circle = Circle(center, size, color="red", fill=False, linewidth=2)
# ax.add_patch(circle)
#
# ax.axis('off')
#
# # Create sliders
# ax_t = plt.axes((0.1, 0.15, 0.8, 0.03))
# ax_slice = plt.axes((0.1, 0.1, 0.8, 0.03))
#
# slider_t = Slider(ax_t, 'Timestep', 0, n_timesteps - 1, valinit=init_t, valstep=1)
# slider_slice = Slider(ax_slice, 'Slice Pos', 0, shape[1] - 1, valinit=init_slice, valstep=1)
#
#
# # Update function
# def update(val):
#     t = int(slider_t.val)
#     slice_pos = int(slider_slice.val)
#
#     # Update image
#     data = np.rot90(np.log10(electron_densities[t][:, slice_pos, :]))
#     im.set_data(data)
#     im.set_clim(vmin=data.min(), vmax=data.max())  # adjust color scaling if needed
#
#     # Update circle
#     circle.set_center((proj_indices[t][0] + offset, proj_indices[t][1]))
#     radius = size * np.sin(0.5 * np.arccos(np.abs(shape[1]/2 - slider_slice.val)/size))
#     # sketchy
#     try:
#         circle.set_radius(int(radius))
#     except:
#         circle.set_radius(0)
#
#     fig.canvas.draw_idle()
#
#
# # Connect sliders to update function
# slider_t.on_changed(update)
# slider_slice.on_changed(update)
#
# plt.show()
#
