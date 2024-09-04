"""
My new proposed discrete model, with a flexible phenotype space.
"""

from dataclasses import dataclass
import numpy as np
from copy import deepcopy
import sim.config.conf as conf
import importlib
import pickle
import time
from sim.phenotype import (
    Phenotype,
    PhenotypeInteractions,
    PhenotypeStructure,
    SequencePhenotypeStructure,
    LatticePhenotypeStructure,
)
import warnings


sequence_chars = tuple("ACDEFGHIKLMNPQRSTVWY")


class UniversalCellParams:
    def __init__(
        self,
        natural_prolif_rate,
        natural_death_rate,
        interaction_induced_rate,
    ):
        # Death rate will be negative for consistency
        self.natural_prolif_base_rate = natural_prolif_rate
        self.natural_death_base_rate = natural_death_rate
        self.interaction_induced_base_rate = interaction_induced_rate


@dataclass
class Selectivities:
    affinity_range: float
    CTL_selectivity: float
    tumour_selectivity: float


class CellBundle:
    def __init__(
        self,
        universal_params: UniversalCellParams,
        phen_struct: PhenotypeStructure,
        cells_at_phenotype: dict,
    ):
        self.cells_at_phenotype = cells_at_phenotype
        self.phen_struct = phen_struct
        self.universal_params = universal_params

    def __len__(self):
        return sum(self.cells_at_phenotype.values())

    def create_cells(self, phenotype: Phenotype, number):

        if phenotype not in self.cells_at_phenotype:
            self.cells_at_phenotype[phenotype] = number
        else:
            self.cells_at_phenotype[phenotype] += number

    def kill_cells(self, phenotype: Phenotype, number):

        if phenotype not in self.cells_at_phenotype:
            raise ValueError(
                "No cells of this phenotype exist. Cannot kill cells. "
                + str(phenotype.id)
            )
        else:
            if self.cells_at_phenotype[phenotype] < number:
                raise ValueError(
                    "Not enough cells of this phenotype exist. Cannot kill cells."
                )
                # self.cells_at_phenotype[phenotype] = 0
            else:
                self.cells_at_phenotype[phenotype] -= number

    def mutate(self, phenotype: Phenotype, number):
        for i in range(number):
            new_phenotype = phenotype.get_random_mutation()
            self.kill_cells(phenotype, 1)
            self.create_cells(new_phenotype, 1)

    @classmethod
    def random(
        self,
        number,
        universal_params: UniversalCellParams,
        phen_struct: PhenotypeStructure,
    ):
        """
        Create a cell bundle with a certain phenotype structure and set of universal cell parameters, randomising the choice of phenotype inside the structure.
        """
        cell_bundle = CellBundle(universal_params, phen_struct, {})
        for i in range(number):
            cell_bundle.create_cells(phen_struct.get_random_phenotype(), 1)
        return cell_bundle

    def verify_phenotype_probabilities(weights):
        """
        TODO: This should no longer be needed, so should be depreciated and removed.
        """
        birth, death, qui = weights
        if qui < 0:
            print("Fixing invalid weights")
            birth += qui / 2
            death += qui / 2
            qui = 0

            if birth < 0:
                death += birth
                birth = 0

            if death < 0:
                birth += death
                death = 0

        return birth, death, qui

    @classmethod
    def evolve_population_loop(
        self,
        cells,
        get_phenotype_probabilities,
    ):
        """
        The default evolve population method. Uses simple logic, but is inefficient.
        """
        new_cells = deepcopy(cells)
        for phenotype, number in cells.cells_at_phenotype.items():

            weights = get_phenotype_probabilities(phenotype)

            rng = np.random.default_rng()
            births, deaths, quiescences = rng.multinomial(number, weights)

            # self.print(births, "|", deaths, "|", quiescences)
            new_cells.create_cells(phenotype, births)
            new_cells.kill_cells(phenotype, deaths)

            # Could just subtract and do this in one step
            new_cells.birth_prob = weights[0]
            new_cells.death_prob = weights[1]
        return new_cells

    @classmethod
    def evolve_population_vec(
        self,
        cells,
        get_phenotype_probabilities,
    ):
        """
        A vectorised version of the evolve population method. Was written for the SequencePhenotypeStructure, but this may no longer work.

        Unused
        TODO: Implement for any phenotype space using more generalised logic.
        """
        new_cells = deepcopy(cells)
        cells_at_phenotype_iterable = list(cells.cells_at_phenotype.items())
        cells_at_phenotype_array = np.empty(
            len(cells_at_phenotype_iterable), dtype=object
        )
        cells_at_phenotype_array[:] = cells_at_phenotype_iterable

        def evolve(tup):
            phenotype, number = tup

            if number != 0:
                weights = get_phenotype_probabilities(phenotype)
                weights = self.verify_phenotype_probabilities(weights)
                rng = np.random.default_rng()
                births, deaths, quiescences = rng.multinomial(number, weights)
                new_cells.create_cells(phenotype, births)
                new_cells.kill_cells(phenotype, deaths)

        vec_evolve = np.vectorize(evolve)
        vec_evolve(cells_at_phenotype_array)
        return new_cells

    @classmethod
    def evolve_population(
        self,
        cells,
        get_phenotype_probabilities,
    ):
        """
        Evolve the population by a single time-step.
        """
        return CellBundle.evolve_population_loop(cells, get_phenotype_probabilities)


class SimulationStateTypes:
    """
    A type specifying what information should be stored in the simulation's history
    """

    @classmethod
    def populations_only(self, state, CTL_cells: CellBundle, tumour_cells: CellBundle):
        """
        Get the lightest SimulationState. Use this if you do not need phenotype proportions, or other details at every time step.
        """
        state.CTL_cells_pop = len(CTL_cells)
        state.tumour_cells_pop = len(tumour_cells)
        return state

    @classmethod
    def whole_cell_bundles(
        self, state, CTL_cells: CellBundle, tumour_cells: CellBundle
    ):
        """
        Get the most detailed SimulationState. Does not throw away any data, but performance may suffer as a result.
        """
        state.CTL_cells_pop = len(CTL_cells)
        state.tumour_cells_pop = len(tumour_cells)
        state.CTL_cells = CTL_cells
        state.tumour_cells = tumour_cells
        return state


class SimulationState:
    """
    A single piece of simulation history.
    """

    type_to_init_dict = {
        "default": SimulationStateTypes.populations_only,
        "detailed": SimulationStateTypes.whole_cell_bundles,
    }

    def __init__(self, CTL_cells: CellBundle, tumour_cells: CellBundle):
        from sim.config.conf import sim_state_init_type

        initialiser = SimulationState.type_to_init_dict[sim_state_init_type]
        self = initialiser(self, CTL_cells, tumour_cells)
        self.init_type = sim_state_init_type


class SimulationHistory:
    """For recording history."""

    def __init__(self, history: list[SimulationState] = []):
        # Do I need to copy the simulation?
        from sim.config.conf import sim_state_init_type

        self.history = history
        self.state_init_type = sim_state_init_type  # We hope this is the same as each sim state, but technically, it could not be

    def update(self, sim_state: SimulationState):
        self.history.append(sim_state)

        # We pickle and unpickle this object to do stuff.

    def __iter__(self):
        return self.history.__iter__()


class Simulation:
    """
    A discrete simulation, functional for any chosen PhenotypeStructures
    """

    def __init__(
        self,
        time_step_size,
        final_time,
        no_init_tumour_cells,
        no_init_CTL_cells,
        tumour_universal_params,
        CTL_universal_params,
        tumour_phenotypic_variation_probability,
        subtype,
        config_name="Unspecified",
        **kwargs,
    ):
        self.config_name = config_name
        self.final_time = final_time
        self.time_step_size = time_step_size
        self.time_step = 0  # An integer describing which time step we're on
        self.final_time_step = int(final_time / time_step_size)

        if subtype == "lattice":
            self.TCR_binding_affinity = kwargs["TCR_binding_affinity"]
            selectivities = kwargs["selectivities"]

            print("Lattice")
            """
            raise NotImplementedError(
                "The vectorisation modification has not yet been made for lattices. Please revert the code to use non `_vec` functions before removing this error statement."
            )
            """
            self.setup_default_lattice_phenotypes(
                kwargs["absolute_max_phenotype"],
                kwargs["no_possible_phenotypes"],
                selectivities,
            )
        elif subtype == "sequence":
            self.sequence_chars = sequence_chars
            self.CTL_sequences = kwargs["CTL_sequences"]
            self.tumour_sequences = kwargs["tumour_sequences"]
            self.sequence_matrix = kwargs["get_sequence_matrix"](self)
            self.affinity_matrix = kwargs["get_affinity_matrix"](self)
            self.binding_scaling = kwargs["binding_scaling"]

            print("Sequence")
            warnings.warn(
                "The sequence matrix functionality has not been implemented, so this will be ignored in the population dynamics."
            )
            self.setup_default_sequence_phenotypes(
                self.tumour_sequences,
                self.CTL_sequences,
                self.sequence_matrix,
                self.affinity_matrix,
                self.binding_scaling,
            )
        else:
            raise NotImplementedError(
                f"The specified subtype {subtype} is not one recognised by the simulation."
            )

        self.tumour_cells = CellBundle.random(
            no_init_tumour_cells, tumour_universal_params, self.tumour_struct
        )
        self.CTL_cells = CellBundle.random(
            no_init_CTL_cells, CTL_universal_params, self.CTL_struct
        )

        self.tumour_phenotypic_variation_probability = (
            tumour_phenotypic_variation_probability
        )

        self.phenotype_tumour_probabilities = {}
        self.phenotype_CTL_probabilities = {}

        self.phenotype_separation_scaling = {}

        self.history = SimulationHistory()

    def setup_default_lattice_phenotypes(
        self,
        absolute_max_phenotype,
        no_possible_phenotypes,
        selectivities: Selectivities,
    ):
        """
        If you are using a lattice phenotype, set up and get the PhenotypeStructure and PhenotypeInteractions
        """
        self.tumour_struct = LatticePhenotypeStructure(
            absolute_max_phenotype, no_possible_phenotypes
        )
        self.CTL_struct = deepcopy(
            self.tumour_struct
        )  # The LatticePhenotypeStructure is designed such that both structures have the same underlying definition, but we have to create two separate instances so that we can specify different parameters for each

        self.phen_int = PhenotypeInteractions.get_default_lattice_interactions(
            self.tumour_struct,
            self.CTL_struct,
            self.TCR_binding_affinity,
            selectivities,
        )

    def setup_default_sequence_phenotypes(
        self,
        possible_tumours: list[str],
        possible_CTLs: list[str],
        sequence_matrix,
        affinity_matrix,
        binding_scaling,
    ):
        """

        Parameters
        ----
        sequence matrix :
            A matrix indicating the compatability/ interactability of two protein sequence letters

        affinity_matrix:
            A matrix indicating the affinities between all the CTLs and tumours present

        binding_scaling:
            A cross-interaction-only additional scaling


        """
        self.tumour_struct = SequencePhenotypeStructure(possible_tumours)
        self.CTL_struct = SequencePhenotypeStructure(possible_CTLs)

        self.phen_int = PhenotypeInteractions.get_default_sequence_interactions(
            self.tumour_struct,
            self.CTL_struct,
            sequence_matrix,
            affinity_matrix,
            binding_scaling,
        )

    def get_immune_score(self):
        """
        Get the immune score. Unused.
        """

        return len(self.CTL_cells.cells) / len(self.tumour_cells.cells)

    def get_average_immune_score(self):
        """
        Get the average immune scoare. Unimplemented
        """
        pass

    def get_phenotype_natural_death_rate_list(
        self, cells: CellBundle, phenotype: Phenotype
    ):
        """
        The default death rate method. Uses simple logic, but is inefficient.
        """

        # Based on death base rate, and a weighted sum of the competition from "close species"
        return cells.universal_params.natural_death_base_rate * sum(
            [
                self.phen_int.get_interaction_scaling(
                    phenotype,
                    other_phenotype,
                )
                * cells_at_phenotype
                for other_phenotype, cells_at_phenotype in cells.cells_at_phenotype.items()
            ]
        )

    def get_phenotype_natural_death_rate_vec(
        self, cells: CellBundle, phenotype: Phenotype
    ):
        """
        A vectorised version of the death rate method. Was written for the SequencePhenotypeStructure, but this may no longer work.

        Unused
        TODO: Implement for any phenotype space using more generalised logic.
        """
        phen_1 = phenotype
        phen_2 = next(iter(cells.cells_at_phenotype.items()))[0]  # assuming non-zero
        struct_tuple = (phen_1.struct, phen_2.struct)
        data = self.phen_int.interaction_data[
            struct_tuple
        ]  # This compares different characters, but isn't directly a matrix of the phenotypes against each other. These will need to be computed at the start, I guess? We assume we'll probably need to use all those datapoints anyway.

        distance_matrix = phen_1.struct.get_distance_matrix(data)
        # print(data)
        vec1 = distance_matrix[phenotype.id, :]
        vec2 = np.zeros(len(phen_2.struct.ids))

        for phen, number in cells.cells_at_phenotype.items():
            vec2[phen.id] = number

        return cells.universal_params.interaction_induced_base_rate * np.dot(vec1, vec2)

    def get_phenotype_natural_death_rate(self, cells: CellBundle, phenotype: Phenotype):
        """
        Compute the current natural death rate of the chosen phenotype due to competition.
        """
        return self.get_phenotype_natural_death_rate_list(cells, phenotype)

    def get_phenotype_interaction_induced_rate_list(
        self,
        cells: CellBundle,
        other_cells: CellBundle,
        phenotype: Phenotype,
    ):
        """
        The default interaction rate method. Uses simple logic, but is inefficient.
        """
        # The rate of growth/ death resulting from the interaction of two sets of cells (tumour and CTL)
        return (
            cells.universal_params.interaction_induced_base_rate
            * self.phen_int.get_cross_term_multipler(
                cells.phen_struct, other_cells.phen_struct
            )
            * sum(
                [
                    self.phen_int.compute_interaction_scaling(
                        phenotype,
                        other_phenotype,
                    )
                    * other_cells_at_phenotype
                    for other_phenotype, other_cells_at_phenotype in other_cells.cells_at_phenotype.items()
                ]
            )
        )

    def get_phenotype_interaction_induced_rate_vec(
        self,
        cells: CellBundle,
        other_cells: CellBundle,
        phenotype: Phenotype,
    ):
        """
        A vectorised version of the interaction rate method. Was written for the SequencePhenotypeStructure, but this may no longer work.

        Unused
        TODO: Implement for any phenotype space using more generalised logic.
        """
        # The rate of growth/ death resulting from the interaction of two sets of cells (tumour and CTL)

        phen_1 = phenotype
        phen_2 = next(iter(other_cells.cells_at_phenotype.items()))[
            0
        ]  # assuming non-zero
        struct_tuple = (phen_1.struct, phen_2.struct)
        data = self.phen_int.interaction_data[struct_tuple]
        # print(data)
        vec1 = data.interaction_matrix[phenotype.id, :]
        vec2 = np.zeros(len(phen_2.struct.ids))

        for phen, number in other_cells.cells_at_phenotype.items():
            vec2[phen.id] = number

        return cells.universal_params.interaction_induced_base_rate * np.dot(vec1, vec2)

    def get_phenotype_interaction_induced_rate(
        self,
        cells: CellBundle,
        other_cells: CellBundle,
        phenotype: Phenotype,
    ):
        """
        Compute the rate of growth/ death of the chosen phenotype induced by interpopulation interaction.
        """

        return self.get_phenotype_interaction_induced_rate_list(
            cells, other_cells, phenotype
        )

    def mutate(self, cells: CellBundle):
        """
        Randomly mutate the cell bundle, with the specified mutation probability.
        """
        new_cells = deepcopy(cells)
        for phenotype, number in cells.cells_at_phenotype.items():
            rng = np.random.default_rng()
            no_mutations = rng.binomial(
                number, self.tumour_phenotypic_variation_probability
            )
            new_cells.mutate(phenotype, no_mutations)
        return new_cells

    def run(self):
        """
        Run the simulation.
        """
        self.print("The simulation is starting.")

        start_time = time.time()
        while self.time_step < self.final_time_step:
            importlib.reload(conf)

            if conf.interrupt:
                print("The simulation has been interrupted and will now safely save.")
                from sim.inputs import reset_interrupt

                reset_interrupt()
                return

            self.time_step += 1

            # Pre-calculation

            # Simulating effects
            self.tumour_cells = self.mutate(self.tumour_cells)

            self.phenotype_tumour_probabilities = {}
            self.phenotype_CTL_probabilities = {}

            self.tumour_cells = CellBundle.evolve_population(
                self.tumour_cells, self.get_phenotype_tumour_probabilities
            )
            self.CTL_cells = CellBundle.evolve_population(
                self.CTL_cells, self.get_phenotype_CTL_probabilities
            )

            self.print("C: ", len(self.tumour_cells), " | T:", len(self.CTL_cells))
            self.print(f"Iteration done after {time.time() - start_time}.")
            self.print("Time step: ", self.time_step, "/", self.final_time_step)

            # Post-calculation
            self.history.update(SimulationState(self.CTL_cells, self.tumour_cells))

            # End it
        self.print("The final time has been reached, so the simulation is over.")

    def print(self, *string):
        """
        Print the string, but in debug mode only.
        """
        if conf.debug:
            print(*string)

    def get_phenotype_tumour_probabilities(self, phenotype):
        if phenotype not in self.phenotype_tumour_probabilities:
            self.phenotype_tumour_probabilities[phenotype] = (
                self.compute_phenotype_tumour_probabilities(phenotype)
            )

        return self.phenotype_tumour_probabilities[phenotype]

    def get_phenotype_CTL_probabilities(self, phenotype):
        if phenotype not in self.phenotype_CTL_probabilities:
            self.phenotype_CTL_probabilities[phenotype] = (
                self.compute_phenotype_CTL_probabilities(phenotype)
            )

        return self.phenotype_CTL_probabilities[phenotype]

    def compute_phenotype_tumour_probabilities(self, phenotype):
        birth = (
            self.time_step_size
            * self.tumour_cells.universal_params.natural_prolif_base_rate
        )
        death = self.time_step_size * (
            self.get_phenotype_natural_death_rate(self.tumour_cells, phenotype)
            + self.get_phenotype_interaction_induced_rate(
                self.tumour_cells, self.CTL_cells, phenotype
            )
        )

        return birth, death, 1 - (birth + death)

    def compute_phenotype_CTL_probabilities(self, phenotype):
        birth = self.time_step_size * (
            self.CTL_cells.universal_params.natural_prolif_base_rate
            + self.get_phenotype_interaction_induced_rate(
                self.CTL_cells, self.tumour_cells, phenotype
            )
        )
        death = self.time_step_size * self.get_phenotype_natural_death_rate(
            self.CTL_cells, phenotype
        )
        return birth, death, 1 - (birth + death)

    def extend(self, additional_time):
        """
        Extend the final time of the simulation.
        """
        self.final_time += additional_time
        self.final_time_step = int(self.final_time / self.time_step_size)

    @classmethod
    def load_simulation(self, path_to_data):
        with open(path_to_data, "rb") as f:
            sim = pickle.load(f)
            print("Successfully opened the previous simulation.")

        return sim

    @classmethod
    def save_simulation(self, path_to_data, sim):
        if input("Append config name?") == "y":
            path_to_data = path_to_data + "_" + sim.config_name

        with open(path_to_data, "wb") as f:
            print("Pickling....")
            pickle.dump(sim, f, pickle.HIGHEST_PROTOCOL)
            print("Pickling done.")


def get_ones_matrix(h: int, w: int):
    """
    Get a generic uniform matrix.
    """
    return np.ones((h, w))


def get_ones_matrix_sequences(sim: Simulation):
    """
    Get a square uniform matrix from a list of sequences.
    """
    l = len(sim.sequence_chars)
    return get_ones_matrix(l, l)


def get_ones_matrix_affinity(sim: Simulation):
    """
    Get a uniform matrix with a height equal to the number of CTL sequences, and a width equal to the number of tumour sequences.
    """
    h = len(sim.CTL_sequences)
    w = len(sim.tumour_sequences)
    return get_ones_matrix(h, w)


def get_random_matrix(h: int, w: int):
    """
    Generate a random matrix with entries in [0,1)
    """
    print(w)
    print(h)
    matrix = np.random.rand(h, w)
    print(matrix)
    return matrix


def get_random_exp_matrix(h: int, w: int, parameter: float = 1):
    """
    Generate an exponential matrix with entries in [0,\infty)
    """
    print(w)
    print(h)
    matrix = np.random.exponential(1 / parameter, (h, w))
    print(matrix)
    return matrix


def get_random_exp_matrix_by_input(h, w):
    """
    Generate a random matrix with entries in the exponential distribution.
    """
    parameter = float(input("What parameter for the distribution would you like?"))
    return get_random_exp_matrix(h, w, parameter)


def get_random_matrix_sequences(sim: Simulation):
    """
    Generate entries from a uniform distribution (this would be a sensible decision).
    """
    l = len(sim.sequence_chars)
    return get_random_matrix(l, l)


def get_random_matrix_affinity(sim: Simulation):
    """
    Generate entries from an exponential distribution (this mimics the idea of not having most phenotypes bind). But we have renormalised and modified to keep within [0,1]
    """
    h = len(sim.CTL_sequences)
    w = len(sim.tumour_sequences)
    exp_matrix = get_random_exp_matrix_by_input(h, w)
    norm_matrix = 1 - (1 / (1 + exp_matrix))
    print(norm_matrix)
    return norm_matrix


def get_identity_matrix_affinity(sim: Simulation):
    """
    Generate a binding probability matrix that is an identity matrix.
    """
    h = len(sim.CTL_sequences)
    w = len(sim.tumour_sequences)
    mat = np.eye(h, w)
    return mat


def get_split_matrix_affinity(sim: Simulation):
    """
    Generate a split matrix as a binding probability matrix.
    """
    k = int(input("How many non-interacting phenotypes do you want?"))
    h = len(sim.CTL_sequences)
    w = len(sim.tumour_sequences)

    row = np.ones(h)
    for i in range(k):
        row[i] = 0
    # First k are 0. Last w-k are 1.

    matrix = np.full((h, w), row)
    print(matrix)
    return matrix
