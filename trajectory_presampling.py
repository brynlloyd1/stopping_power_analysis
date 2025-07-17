from ase.lattice.cubic import FaceCenteredCubic

import numpy as np


class Trajectory:
    def __init__(self, starting_position, direction, dx, n_samples):
        self.starting_position = starting_position
        self.direction = direction
        self._dx = dx
        self._n_samples = n_samples
        self.trajectory = self.create_trajectory_sample()

    def create_trajectory_sample(self):
        trajectory = np.array([self.starting_position + self.direction * self._dx * i for i in range(self._n_samples)])
        return trajectory


class TrajectorySampler:
    """
    should calculate P_ideal using the mini unit cell
    but should calculate P_traj for the exact smame supercell used in the simulations
    use center() or not to use center()
        need to use center() otherwise you get weird results at the edges of the bounding box
        but might need to get rid of center() in the ehrenfest sims
    """

    def __init__(self):
        self.LATTICE_CONSTANT = 4.05

        # attributes for calculating P_ideal
        self.Pideal_nsamples = 50000  # paper uses 50k

        # attributes for calculating P_traj
        self.dx = 0.02  # Angstrom - used in paper
        self.path_length = 150  # Angstrom - used in paper
        self.Ptraj_nsamples = int(self.path_length // self.dx)

        self.create_unitcell()
        self.P_ideal = self.calculate_Pideal()

        self.supercell_size = (3, 3, 3)
        self.create_supercell()

        # attributes for calculating the trajectory metric
        self.bin_width = 0.05  # Angstrom - paper uses bins of 0.01 Angstrom
        # idk if using the lattice constant is correct here
        self.bins = np.arange(0., 0.5 * self.LATTICE_CONSTANT, self.bin_width)
        self.n_bins = len(self.bins)

    # SETTERS #
    # BC MANY THINGS CHANGE WHEN VARIABLES ARE UPDATED
    # TODO: replace the things that change with @property decorator??

    def set_path_length(self, new_path_length):
        """ number of samples must also be updated, so cannot just change the attribute directly """
        self.path_length = new_path_length
        self.Ptraj_nsamples = int(new_path_length // self.dx)

    def set_bin_width(self, new_bin_width):
        """ the default makes nice plots, but you can change it if you really want to """
        self.bin_width = new_bin_width
        self.bins = np.arange(0., 0.5 * self.LATTICE_CONSTANT, self.bin_width)
        self.n_bins = len(self.bins)

    def set_supercell(self, supercell_size):
        self.supercell_size = supercell_size
        self.create_supercell()

    def create_unitcell(self):
        """
        doesn't actually create a unit cell so should probably rename it
        creates a full cube of fcc bc otherwise distances to the nearest nucleus
        can be incorrect - but then keeps the cell size to 4.05 Angstroms
        """

        # fcc cube is 1.5 unit cells...
        self.unit_cell = FaceCenteredCubic("Al",
                                           size=(2, 2, 2),
                                           latticeconstant=self.LATTICE_CONSTANT)

        cell = self.unit_cell.get_cell()
        limit = 1.5 * cell.diagonal() / 2
        mask = np.all(self.unit_cell.positions < limit, axis=1)
        self.unit_cell = self.unit_cell[mask]
        self.unit_cell.set_cell(cell / 2, scale_atoms=False)

        # self.unit_cell.center()
        self.unit_cell_positions = self.unit_cell.get_positions()

    def create_supercell(self):
        self.supercell = FaceCenteredCubic("Al",
                                           size=self.supercell_size,
                                           latticeconstant=self.LATTICE_CONSTANT)
        self.supercell.center()
        self.supercell_positions = self.supercell.get_positions()

    def find_nearest_neighbour(self, point, lattice_positions):
        """
        takes the lattice positions as an argument bc you might want to pass it either a 'unitcell' or a big supercell depending on the situation
        """
        distances = np.array([np.linalg.norm(lattice_position - point)
                              for lattice_position in lattice_positions])

        return np.min(distances), np.argmin(distances)

    def create_trajectory(self, kind="hyperchannelling", params=None) -> Trajectory:
        """
        creates an instance of Trajectory for a given kind of trajectory

        Parameters:
        kind (str): one of "hyperchannelling" (default), "centroid", "random"
        params (Dict[str, np.ndarray[float64]]): parameters for creating a user defined trajectory. Has keys "starting_position" and "direction"

        Returns:
        trajectory (Trajectory)
        """

        if kind == "hyperchannelling":
            # starting_position = np.array([0.0, 0.75 * self.LATTICE_CONSTANT, 0.75 * self.LATTICE_CONSTANT])
            starting_position = np.array([0.0, self.LATTICE_CONSTANT, self.LATTICE_CONSTANT])
            direction = np.array([1.0, 0.0, 0.0])


        elif kind == "centroid":
            print("centroid trajectory will be broken since Ive changed to exclusively using supercells")
            vertex1 = [0., 0.]
            vertex2 = [0.25 * self.LATTICE_CONSTANT, 0]
            vertex3 = [0.25 * self.LATTICE_CONSTANT,
                       0.25 * self.LATTICE_CONSTANT]
            centroid_pos = [(vertex1[0] + vertex2[0] + vertex3[0]) / 3 + self.LATTICE_CONSTANT,
                            (vertex1[1] + vertex2[1] + vertex3[1]) / 3 + self.LATTICE_CONSTANT]
            starting_position = np.array([0.0, centroid_pos[0], centroid_pos[1]])
            direction = np.array([1.0, 0.0, 0.0])

        elif kind == "random":  # else create random trajectory
            starting_position = np.random.uniform(
                0, self.LATTICE_CONSTANT, size=3)
            direction = np.random.uniform(-1, 1, size=3)

        else:  # kind == "custom"
            if not params:
                raise ValueError("cant have a custom trajectory without specifying it in params")

            starting_position = params["starting_position"]
            direction = params["direction"]

        # direction should be a unit vector
        direction /= np.linalg.norm(direction)

        return Trajectory(starting_position, direction, self.dx, self.Ptraj_nsamples)

    def calculate_Pideal(self):
        """
        calculates the distribution of nearest neighbour distances for an indeal trajectory
        where an ideal trajectory is simulated by sampling self.P_ideal_nsamples points inside a unit cell
        """

        cell = self.unit_cell.get_cell()
        x_min = np.min(cell[:, 0])
        x_max = np.max(cell[:, 0])
        y_min = np.min(cell[:, 1])
        y_max = np.max(cell[:, 1])
        z_min = np.min(cell[:, 2])
        z_max = np.max(cell[:, 2])

        # Generate uniformly distributed values per axis
        x = np.random.uniform(x_min, x_max, self.Pideal_nsamples)
        y = np.random.uniform(y_min, y_max, self.Pideal_nsamples)
        z = np.random.uniform(z_min, z_max, self.Pideal_nsamples)

        # Stack into an (n, 3) array of 3D positions
        positions = np.column_stack((x, y, z))

        P_ideal = np.array([])
        for point in positions:
            nn_dist, _ = self.find_nearest_neighbour(
                point, self.unit_cell.positions)
            P_ideal = np.append(P_ideal, nn_dist)
        return P_ideal

    def calculate_Ptraj_supercell(self, kind="hyperchanneling", params=None):
        """
        This function should be used for calculating P_traj
        :param kind:
        :param params:
        :return:
        """
        trajectory = self.create_trajectory(kind=kind, params=params)
        trajectory_supercell = np.mod(trajectory.trajectory, self.supercell.get_cell().diagonal())

        P_traj = np.array([])
        for sample_point in trajectory_supercell:
            nn_dist, _ = self.find_nearest_neighbour(
                sample_point, self.supercell_positions)
            P_traj = np.append(P_traj, nn_dist)

        return P_traj, trajectory

    def calculate_Ptraj_unitcell(self, kind="hyperchanneling", params=None):
        """
        not to be used... use calculate_Ptraj_supercell instead, and make the supercell the same as the one used in the simulation
        """
        trajectory = self.create_trajectory(kind=kind, params=params)
        trajectory_unitcell = np.mod(trajectory.trajectory, self.unit_cell.get_cell().diagonal())

        P_traj = np.array([])
        for sample_point in trajectory_unitcell:
            nn_dist, _ = self.find_nearest_neighbour(
                sample_point, self.unit_cell_positions)
            P_traj = np.append(P_traj, nn_dist)

        return P_traj, trajectory

    def bin_nn_distribution(self, P):
        return np.histogram(P, bins=self.bins)

    def calculate_trajectory_metric(self, P_traj, cutoff=-1):
        """
        probably something wrong with normalising the binned data
        and then not multiplying by bin width in calculation
        not sure tho but it seems to work
        """
        P_ideal_binned, self.bin_edges = self.bin_nn_distribution(self.P_ideal)
        P_traj_binned, _ = self.bin_nn_distribution(P_traj[:cutoff])

        P_ideal_binned = P_ideal_binned.astype(
            np.float64) / np.sum(P_ideal_binned)
        P_traj_binned = P_traj_binned.astype(
            np.float64) / np.sum(P_traj_binned)

        # ----- numerical integration to calculate trajectory metric -----
        D_H_square = 1
        # iterate over bins
        for i in range(len(self.bin_edges) - 1):
            D_H_square -= np.sqrt(P_traj_binned[i] * P_ideal_binned[i])
        # ----- numerical integration to calculate trajectory metric -----

        return P_ideal_binned, P_traj_binned, np.sqrt(D_H_square)
