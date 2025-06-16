import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from gpaw import restart


class GroundStateAnalysis:
    def __init__(self, path_to_file):
        self.atoms, self.calc = restart(path_to_file)
        self.lattice_positions = self.atoms.get_positions()   # also includes the projectile...# also includes the projectile...# also includes the projectile...
        self.electron_density = self.calc.get_all_electron_density()

    def visualise_electron_density(self, plot_nuclei=True):
        nx, ny, nz = np.shape(self.electron_density)
        print(f"electron density shape: ({nx}, {ny}, {nz})")
        cell = np.diag(self.atoms.get_cell())

        x = np.linspace(0, cell[0], nx)
        y = np.linspace(0, cell[1], ny)
        z = np.linspace(0, cell[2], nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

        # initial values for the slider
        init_lower = 7
        init_upper = 15

        mask = (self.electron_density >= init_lower) & (self.electron_density <= init_upper)
        X_plot = X[mask]
        Y_plot = Y[mask]
        Z_plot = Z[mask]
        electron_density_plot = self.electron_density[mask]

        fig = plt.figure(figsize=(50, 10))
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlim((0, cell[0]))
        ax.set_ylim((0, cell[1]))
        ax.set_zlim((0, cell[2]))
        # not convinced this line does anything. Is supposed to make equal aspect ratio
        ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))

        if plot_nuclei:
            fcc_xs = self.lattice_positions[:,0]
            fcc_ys = self.lattice_positions[:,1]
            fcc_zs = self.lattice_positions[:,2]

            ax.scatter(fcc_xs, fcc_ys, fcc_zs, s=100, color="C0")


        scat = ax.scatter(X_plot, Y_plot, Z_plot,
                          c=electron_density_plot, cmap="jet")

        fig.colorbar(scat)

        ax_lower = plt.axes((0.25, 0.1, 0.65, 0.03))
        ax_upper = plt.axes((0.25, 0.15, 0.65, 0.03))
        slider_lower = Slider(ax_lower, "lower threshold", 0, 20, valinit=init_lower, valstep=0.01)
        slider_upper = Slider(ax_upper, "upper threshold", 0, 20, valinit=init_upper, valstep=0.01)

        def update(val):
            lower = slider_lower.val
            upper = slider_upper.val
            if upper < lower:
                upper = lower
                slider_upper.set_val(upper)

            mask = (self.electron_density >= lower) & (self.electron_density <= upper)
            for coll in ax.collections: coll.remove()

            ax.scatter(X[mask], Y[mask], Z[mask], c=self.electron_density[mask], cmap="jet")

            if plot_nuclei:
                ax.scatter(fcc_xs, fcc_ys, fcc_zs, s=100, color="C0")

            fig.canvas.draw_idle()

        slider_lower.on_changed(update)
        slider_upper.on_changed(update)

        plt.show()


    def calculate_projectile_charge_state(self):
        projectile_position = self.atoms.get_positions()[-1]
        cutoff_radius = 0.2   # Angstroms
        charge_state = 0

        cell_size = np.diag(self.atoms.get_cell())
        nx, ny, nz = np.shape(self.electron_density)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x_pos = i/nx * cell_size[0]
                    y_pos = j/ny * cell_size[1]
                    z_pos = k/nz * cell_size[2]
                    position_in_space = np.array([x_pos, y_pos, z_pos])

                    distance_to_projectile = np.sqrt(
                        (position_in_space[0] - projectile_position[0])**2 +
                        (position_in_space[1] - projectile_position[1])**2 +
                        (position_in_space[2] - projectile_position[2])**2)

                    if distance_to_projectile <= cutoff_radius:
                        charge_state += self.electron_density[i,j,k]

        return charge_state





if __name__ == "__main__":
    path_to_file = ("/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/dft_data/Al_ground_state_far.gpw")
    a = GroundStateAnalysis(path_to_file)
    print(a.calculate_projectile_charge_state())
    a.visualise_electron_density()
