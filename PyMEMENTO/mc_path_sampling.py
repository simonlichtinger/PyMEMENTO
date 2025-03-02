""" This contains the MCPathSampler class which handles finding the minimum RMSD path through the models. """

import time
import multiprocessing
import shutil
import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd


def run_one_replicate_of_mc(seed, sampler):
    """Wrapper function in order to call run_monte_carlo from a multiprocessing pool.
    It needs to lie outside the class so that it is picklable and can be handled with the
    standard multiprocessing library.

    :param seed: Seed for the random number generator.
    :type seed: int
    :param sampler: The sampler object to run.
    :type sampler: :class:`MCPathSampler`

    :return: The pathsampling output as from run_monte_carlo.
    :rtype: tuple"""
    start_time = time.time()
    print("Replicate is initialised")
    optimum = sampler.run_monte_carlo(
        f"MC_{multiprocessing.current_process().name}_seed{seed}/",
        sampler.MC_STEPS,
        annealing_progression=[
            sampler.STARTING_TEMPERATURE / (float(n) / 500)
            for n in range(500, sampler.MC_STEPS + 500)
        ],
        seed=seed,
    )
    finished_time = time.time()
    print(f"Time taken: {finished_time-start_time}")
    print(optimum)
    return optimum


class MCPathSampler:
    """The MCPathSampler class handles finding a minumum RMSD path through the models which were generated by modeller."""

    def __init__(
        self,
        folder: str,
        number_of_windows: int,
        models_per_window: int,
        starting_temperature: int,
        mc_steps: int,
        number_of_replicates: int,
    ):
        """Initialise the MCPathSampler with all relevant information.

        :param folder: Folder path in which modeller was active.
        :type folder: str
        :param number_of_windows: How many windows are there.
        :type number_of_windows: int
        :param models_per_window: How many models have been built per window.
        :type models_per_window: int
        :param starting_temperature: At what (fictious) temperate should the annealing start.
        :type starting_temperature: int
        :param mc_steps: How many MC steps to perform (half-life of temperature is 500).
        :type mc_steps: int
        :param number_of_replicates: How many MC replicates should be run.
        :type number_of_replicates: int
        """
        # intialise parameters used below
        self.FOLDER = folder
        self.NUMBER_OF_WINDOWS = number_of_windows
        self.MODELS_PER_WINDOW = models_per_window
        self.STARTING_TEMPERATURE = starting_temperature
        self.MC_STEPS = mc_steps
        self.NUMBER_OF_REPLICATES = number_of_replicates

        # Load in all files as MDAnalysis universes
        self.universes = []
        for i in range(self.NUMBER_OF_WINDOWS):
            self.universes.append([])
            for j in range(1, self.MODELS_PER_WINDOW + 1):
                self.universes[-1].append(
                    mda.Universe(
                        self.FOLDER + f"morph{i}/protein.B9999{str(j).zfill(4)}.pdb", in_memory=True
                    )
                )

    def path_energy_function(self, path):
        """Computes the sum of squares of the RMSD between neighbours in a given path. Serves
        as the energy function in the MC sampling.

        :param path: Path of :class:`MDAnalysis.Universe` objects.
        :type path: list

        :return: Sum of squares of RMSD between neighbours, acting as energy function.
        :rtype: float"""
        energy = 0
        for i in range(len(path) - 1):
            energy += np.square(
                rmsd(
                    self.universes[i][path[i]].atoms.positions,
                    self.universes[i + 1][path[i + 1]].atoms.positions,
                )
            )
        return energy

    def __metropolis(self, delta_energy, temperature):
        """Simple metropolis criterion."""
        if delta_energy <= 0:
            return True
        elif np.random.random_sample() <= np.exp(-delta_energy / temperature):
            return True
        else:
            return False

    def __initialize_path(self):
        """Proposes an initial random path."""
        return [
            np.random.randint(self.MODELS_PER_WINDOW)
            for i in range(self.NUMBER_OF_WINDOWS)
        ]

    def __propose_step(self, path):
        """This proposes a step in the path with a uniform transition attempt probability matrix.
        The step is change of a single value in the path list. Returns tuple: (new_path, changed_index)."""
        new_path = path.copy()
        while new_path == path:
            changed_index = np.random.randint(self.NUMBER_OF_WINDOWS)
            new_path[changed_index] = np.random.randint(self.MODELS_PER_WINDOW)
        return (new_path, changed_index)

    def run_monte_carlo(
        self,
        name,
        number_of_steps,
        annealing_progression=None,
        path=None,
        verbose=False,
        seed=None,
    ):
        """This is a simple MC run procedure for our system. Returns the minimum energy configuration as
        a tuple (min_configuration,  min_energy).

        :param name: Name of the folder in which to store the results.
        :type name: str
        :param number_of_steps: How many MC steps to perform.
        :type number_of_steps: int
        :param annealing_progression: Stepwise list of temperatures to use for annealing. If None, no annealing is performed. Defaults to None.
        :type annealing_progression: list, optional
        :param path: Initial path to start from. If None, a random path is generated. Defaults to None.
        :type path: list, optional
        :param verbose: Whether to print out the progress of the MC run. Defaults to False.
        :type verbose: bool, optional
        :param seed: Seed for the random number generator. If None, gets randomly assigned from python random generator. Defaults to None.
        :type seed: int, optional

        :return: Minimum energy configuration (path as list of Universes) and its energy.
        :rtype: tuple"""

        # prepare folder
        local_path = self.FOLDER + name
        os.makedirs(local_path, exist_ok=True)

        # Initialise the system
        if path is None:
            path = self.__initialize_path()
        if not seed is None:
            np.random.seed(seed)
        energy = self.path_energy_function(path)
        temperature = self.STARTING_TEMPERATURE

        trajectory = []
        energy_progression = []
        accepted_moves = []

        # Run loop
        print(f"Starting MC run with {number_of_steps} steps")
        for n in range(number_of_steps):
            # update temperature, if annealing is on
            if not annealing_progression is None:
                temperature = annealing_progression[n]

            # Monte Carlo logic
            (new_path, changed_index) = self.__propose_step(path)
            delta_energy = self.path_energy_function(new_path) - energy
            if self.__metropolis(delta_energy, temperature):
                path = new_path
                energy += delta_energy
                accepted_moves.append(True)
            else:
                accepted_moves.append(False)
            if verbose:
                print(f"Step {n}: {energy}")
            # Handling trajectory output
            trajectory.append(path)
            energy_progression.append(energy)

        print("Done MC run. Dumping trajectory and energies.")
        plt.plot(energy_progression)
        plt.ylabel("Energy proxy for path RMSD of RMSDs")
        plt.xlabel("MC step")
        plt.savefig(local_path + "trajectory.pdf")
        plt.close()

        # Write full trajectory and energies to file
        np.savetxt(local_path + "mc_energies.dat", energy_progression)
        with open(local_path + "mc_trajectory.dat", "w") as f:
            for configuration in trajectory:
                f.write(str(configuration) + "\n")

        # Return minimum configuration and energy, as those are relevant for our optmisation problem
        min_energy = np.min(energy_progression)
        min_index = energy_progression.index(min_energy)
        min_configuration = trajectory[min_index]

        return (min_configuration, min_energy)

    def sample_path(self, poolsize=1):
        """Perform the MC sampling using a multiprocessing pool of the desired size.

        :param poolsize: How many multiprocessing instances to use. Setting 1 bypasses multiprocessing, defaults to 1
        :type poolsize: int, optional
        """
        # Make an run a multiprocessing pool with random seeds.
        # Using starmap because we need to pass on the current instance of the sampler.
        # I'm using a super-ugly workaround here where if poolsize==1 no pool at all
        # is used. This is to avoid a strange memory bug in MDAnalysis for large systems.
        if poolsize > 1:
            pool = multiprocessing.Pool(poolsize)
            solutions = pool.starmap(
                run_one_replicate_of_mc,
                [
                    (np.random.randint(100000), self)
                    for n in range(self.NUMBER_OF_REPLICATES)
                ],
            )
            pool.close()
            pool.join()
        else:
            solutions = [
                run_one_replicate_of_mc(np.random.randint(100000), self)
                for n in range(self.NUMBER_OF_REPLICATES)
            ]

        # Extract overall best solution
        best_solution_energy = min([sol[1] for sol in solutions])
        best_solution_index = [sol[1] for sol in solutions].index(best_solution_energy)
        best_solution_configuration = solutions[best_solution_index][0]

        # Write output to terminal
        print("----------------------------------")
        print("Found best solution:")
        print(best_solution_energy, best_solution_configuration)

        # Make 'best.pdb' files corresponding to the optimal path within the morph directories.
        for i, val in enumerate(best_solution_configuration):
            shutil.copy(
                self.FOLDER + f"morph{i}/protein.B9999{str(val+1).zfill(4)}.pdb",
                self.FOLDER + f"morph{i}/best.pdb",
            )
