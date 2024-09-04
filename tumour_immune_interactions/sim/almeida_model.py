"""
An implementation of the model Almeida described, complete with its parameters and interaction logic. This was the model that I based my modifications of off.
"""

import random
import numpy as np
from copy import deepcopy
import sim.config.conf as conf
import importlib
import pickle
from typing import TypeVar


# A set of workarounds for not having Self, from Python 3.11 (change back to self if you upgrade)
TPhenotypeStructure = TypeVar("TPhenotypeStructure", bound="PhenotypeStructure")
TCellBundle = TypeVar("TCellBundle", bound="CellBundle")
TSimulation = TypeVar("TSimulation", bound="Simulation")


class PhenotypeStructure:
    """
    Phenotypes are in the range [0,1], and are floats
    Phenotype IDs are in the range [0, no_possible_values - 1], and are integers
    """

    def __init__(self, abs_max_value, no_possible_values):
        self.abs_max_value = abs_max_value
        self.range = np.linspace(
            -self.abs_max_value, self.abs_max_value, no_possible_values, True
        )
        self.id_range = range(no_possible_values)
        self.step_size = 2 * abs_max_value / no_possible_values  # check this
        self.no_possible_values = no_possible_values

    """
    def shift(self, phen, no_steps, direction):
        print(phen)
        phen += no_steps * direction * self.step_size
        if phen > self.abs_max_value:
            phen = self.abs_max_value
        if phen < -self.abs_max_value:
            phen = -self.abs_max_value
        return phen
    """

    def shift(self, phen_id, no_steps, direction):
        phen_id += no_steps * direction
        if phen_id > self.no_possible_values - 1:
            phen_id = self.no_possible_values - 1
        elif phen_id < 0:
            phen_id = 0
        return phen_id

    def get_phenotype_by_id(self, id):
        return (id * self.step_size) - self.abs_max_value

    def get_random_phenotype_id(self):
        return random.randint(0, self.no_possible_values - 1)

    @classmethod
    def is_excluded_phenotype(
        self, phen_struct: TPhenotypeStructure, phenotype_id: int, exclude_percent=0.1
    ):
        no_values = phen_struct.no_possible_values
        return (
            phenotype_id < no_values * exclude_percent
            or phenotype_id > no_values * (1 - exclude_percent)
        )

    """
    def get_random_phenotype(self):
        return (
            random.randint(0, self.no_possible_values - 1) * self.step_size
        ) - self.abs_max_value
    """


class UniversalCellParams:
    def __init__(
        self,
        natural_prolif_rate,
        natural_death_rate,
        interaction_induced_rate,
        selectivity,
    ):
        # Death rate will be negative for consistency
        self.natural_prolif_base_rate = natural_prolif_rate
        self.natural_death_base_rate = natural_death_rate
        self.interaction_induced_base_rate = interaction_induced_rate
        self.selectivity = selectivity


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

    def create_cells(self, phenotype_id, number):
        """
        if PhenotypeStructure.is_excluded_phenotype(self.phen_struct, phenotype_id):
            return
        """

        if phenotype_id not in self.cells_at_phenotype:
            self.cells_at_phenotype[phenotype_id] = number
        else:
            self.cells_at_phenotype[phenotype_id] += number

    def kill_cells(self, phenotype_id, number):
        if phenotype_id not in self.cells_at_phenotype:
            raise ValueError("No cells of this phenotype exist. Cannot kill cells.")
        else:
            if self.cells_at_phenotype[phenotype_id] < number:
                raise ValueError(
                    "Not enough cells of this phenotype exist. Cannot kill cells."
                )
            else:
                self.cells_at_phenotype[phenotype_id] -= number

    def mutate_int(self, phenotype_id, number, no_steps, direction):
        new_phenotype_id = self.phen_struct.shift(
            phenotype_id, no_steps, direction
        )  # If on the boundary, the new phenotype is the _same_ as the old phenotype. Indeed, this aligns with Almeida requiring "abortion" of the mutation if it moves out of the region.

        self.kill_cells(phenotype_id, number)
        """
        if PhenotypeStructure.is_excluded_phenotype(self.phen_struct, new_phenotype_id):
            return
        """
        self.create_cells(new_phenotype_id, number)

    def mutate_left(self, phenotype_id, number):
        self.mutate_int(phenotype_id, number, 1, -1)

    def mutate_right(self, phenotype_id, number):
        self.mutate_int(phenotype_id, number, 1, 1)

    @classmethod
    def random(
        self,
        number,
        universal_params: UniversalCellParams,
        phen_struct: PhenotypeStructure,
    ):
        cell_bundle = CellBundle(universal_params, phen_struct, {})
        for i in range(number):
            cell_bundle.create_cells(phen_struct.get_random_phenotype_id(), 1)
        return cell_bundle

    @classmethod
    def evolve_population(
        self,
        cells: TCellBundle,
        get_phenotype_probabilities,
    ):
        new_cells = deepcopy(cells)
        for phenotype_id, number in cells.cells_at_phenotype.items():
            # print(number)
            # number = 100 * number
            """
            if PhenotypeStructure.is_excluded_phenotype(
                cells.phen_struct, phenotype_id
            ):
                continue
            """
            weights = get_phenotype_probabilities(phenotype_id)
            # print(weights)
            rng = np.random.default_rng()
            births, deaths, quiescences = rng.multinomial(number, weights)
            # print(births, "|", deaths, "|", quiescences)
            new_cells.create_cells(phenotype_id, births)
            new_cells.kill_cells(phenotype_id, deaths)
            # print(len(new_cells))
            # Could just subtract and do this in one step
        return new_cells


class SimulationStateTypes:
    @classmethod
    def populations_only(
        self, state: TCellBundle, CTL_cells: CellBundle, tumour_cells: CellBundle
    ):
        state.CTL_cells_pop = len(CTL_cells)
        state.tumour_cells_pop = len(tumour_cells)
        return state

    @classmethod
    def whole_cell_bundles(
        self, state, CTL_cells: CellBundle, tumour_cells: CellBundle
    ):
        state.CTL_cells_pop = len(CTL_cells)
        state.tumour_cells_pop = len(tumour_cells)
        state.CTL_cells = CTL_cells
        state.tumour_cells = tumour_cells
        return state


class SimulationState:
    type_to_init_dict = {
        "default": SimulationStateTypes.populations_only,
        "detailed": SimulationStateTypes.whole_cell_bundles,
    }

    def __init__(
        self,
        CTL_cells: CellBundle,
        tumour_cells: CellBundle,
    ):
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
    def __init__(
        self,
        time_step_size,
        final_time,
        no_possible_phenotypes,
        absolute_max_phenotype,
        no_init_tumour_cells,
        no_init_CTL_cells,
        tumour_universal_params,
        CTL_universal_params,
        TCR_affinity_range,
        TCR_binding_affinity,
        tumour_phenotypic_variation_probability,
        config_name="Unspecified",
    ):
        self.config_name = config_name
        self.final_time = final_time
        if conf.m_adjustment:
            scalar = 1  # 1.15, 2
            time_step_size *= scalar
            final_time *= scalar
        self.time_step_size = time_step_size
        self.time_step = 0  # An integer describing which time step we're on
        self.final_time_step = int(final_time / time_step_size)

        self.phen_struct = PhenotypeStructure(
            absolute_max_phenotype, no_possible_phenotypes
        )
        self.tumour_cells = CellBundle.random(
            no_init_tumour_cells, tumour_universal_params, self.phen_struct
        )
        self.CTL_cells = CellBundle.random(
            no_init_CTL_cells, CTL_universal_params, self.phen_struct
        )

        self.TCR_affinity_range = TCR_affinity_range
        self.TCR_binding_affinity = TCR_binding_affinity
        self.tumour_phenotypic_variation_probability = (
            tumour_phenotypic_variation_probability
        )

        self.phenotype_tumour_probabilities = {}
        self.phenotype_CTL_probabilities = {}
        self.phenotype_separation_scaling = {}

        self.history = SimulationHistory()
        self.temp_scalar = 1

    def get_immune_score(self):
        return len(self.CTL_cells.cells) / len(self.tumour_cells.cells)

    def get_average_immune_score(self):
        pass

    def compute_interaction_scaling(
        self, phenotype_1_id, phenotype_2_id, range, distance_type="line"
    ):
        """
        Return 0 if phenotypes are out of the selectivity range; otherwise return a constant, modified to account for the boundary.
        """
        phenotype_1 = self.phen_struct.get_phenotype_by_id(phenotype_1_id)
        phenotype_2 = self.phen_struct.get_phenotype_by_id(phenotype_2_id)
        if distance_type == "line":
            distance = abs(phenotype_1 - phenotype_2)
        if distance_type == "circular":
            distance = Simulation.get_circular_distance(phenotype_1, phenotype_2)
        if distance <= range:
            return 1 / (
                min(phenotype_1 + range, self.phen_struct.abs_max_value)
                - max(phenotype_1 - range, -self.phen_struct.abs_max_value)
            )
        else:
            return 0

    @classmethod
    def get_circular_distance(self, a, b):
        """
        Get the circular distance on [0,1] between a and b.
        """
        if b < a:
            # Switch so a < b
            temp = a
            a = b
            b = temp
        return min(abs(b - a), abs(a + 1 - b))

    def get_interaction_scaling(
        self,
        phenotype_1_id,
        phenotype_2_id,
        range,
    ):
        if (
            phenotype_1_id,
            phenotype_2_id,
        ) not in self.phenotype_separation_scaling:
            self.phenotype_separation_scaling[(phenotype_1_id, phenotype_2_id)] = (
                self.compute_interaction_scaling(phenotype_1_id, phenotype_2_id, range)
            )
        return self.phenotype_separation_scaling[(phenotype_1_id, phenotype_2_id)]

    def get_phenotype_natural_death_rate(self, cells: CellBundle, phenotype_id):
        # Based on death base rate, and a weighted sum of the competition from "close species"
        return cells.universal_params.natural_death_base_rate * sum(
            [
                self.get_interaction_scaling(
                    phenotype_id,
                    other_phenotype_id,
                    cells.universal_params.selectivity,
                )
                * cells_at_phenotype
                for other_phenotype_id, cells_at_phenotype in cells.cells_at_phenotype.items()
            ]
        )

    def get_phenotype_interaction_induced_rate(
        self,
        cells: CellBundle,
        other_cells: CellBundle,
        phenotype_id,
    ):
        # The rate of growth/ death resulting from the interaction of two sets of cells (tumour and CTL)
        return (
            cells.universal_params.interaction_induced_base_rate
            * self.TCR_binding_affinity
            * sum(
                [
                    self.get_interaction_scaling(
                        phenotype_id,
                        other_phenotype_id,
                        self.TCR_affinity_range,
                    )
                    * other_cells_at_phenotype
                    for other_phenotype_id, other_cells_at_phenotype in other_cells.cells_at_phenotype.items()
                ]
            )
        )

    def mutate(self, cells: CellBundle):
        new_cells = deepcopy(cells)
        for phenotype_id, number in cells.cells_at_phenotype.items():
            rng = np.random.default_rng()
            mutations = rng.binomial(
                number, self.tumour_phenotypic_variation_probability
            )
            mutate_lefts, mutate_rights = rng.multinomial(mutations, [1 / 2.0] * 2)
            new_cells.mutate_left(phenotype_id, mutate_lefts)
            new_cells.mutate_right(phenotype_id, mutate_rights)
        return new_cells

    def run(self):
        self.print("The simulation is starting.")
        while self.time_step < self.final_time_step:
            importlib.reload(conf)

            if conf.interrupt:
                print("The simulation has been interrupted and will now safely save.")
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
            self.print("Iteration done.")
            self.print("Time step: ", self.time_step, "/", self.final_time_step)
            # Post-calculation

            self.history.update(SimulationState(self.CTL_cells, self.tumour_cells))

            # End it
        self.print("The final time has been reached, so the simulation is over.")

    def print(self, *string):
        if conf.debug:
            print(*string)

    def get_phenotype_tumour_probabilities(self, phenotype_id):
        if phenotype_id not in self.phenotype_tumour_probabilities:
            self.phenotype_tumour_probabilities[phenotype_id] = (
                self.compute_phenotype_tumour_probabilities(phenotype_id)
            )

        return self.phenotype_tumour_probabilities[phenotype_id]

    def get_phenotype_CTL_probabilities(self, phenotype_id):
        if phenotype_id not in self.phenotype_CTL_probabilities:
            self.phenotype_CTL_probabilities[phenotype_id] = (
                self.compute_phenotype_CTL_probabilities(phenotype_id)
            )

        return self.phenotype_CTL_probabilities[phenotype_id]

    def compute_phenotype_tumour_probabilities(self, phenotype_id):
        birth = (
            self.time_step_size
            * self.tumour_cells.universal_params.natural_prolif_base_rate
        )
        death = self.time_step_size * (
            self.get_phenotype_natural_death_rate(self.tumour_cells, phenotype_id)
            + self.get_phenotype_interaction_induced_rate(
                self.tumour_cells, self.CTL_cells, phenotype_id
            )
        )
        birth = self.temp_scalar * birth
        death = self.temp_scalar * death
        return birth, death, 1 - (birth + death)

    def compute_phenotype_CTL_probabilities(self, phenotype_id):
        birth = self.time_step_size * (
            self.CTL_cells.universal_params.natural_prolif_base_rate
            + self.get_phenotype_interaction_induced_rate(
                self.CTL_cells, self.tumour_cells, phenotype_id
            )
        )
        death = self.time_step_size * self.get_phenotype_natural_death_rate(
            self.CTL_cells, phenotype_id
        )
        birth = self.temp_scalar * birth
        death = self.temp_scalar * death
        return birth, death, 1 - (birth + death)

    def extend(self, additional_time):
        self.final_time += additional_time
        self.final_time_step = int(self.final_time / self.time_step_size)

    @classmethod
    def load_simulation(self, path_to_data) -> TSimulation:
        with open(path_to_data, "rb") as f:
            sim = pickle.load(f)
            print("Successfully opened the previous simulation.")

        return sim

    @classmethod
    def save_simulation(self, path_to_data, sim):
        with open(path_to_data, "wb") as f:
            print("Pickling....")
            pickle.dump(sim, f, pickle.HIGHEST_PROTOCOL)
            print("Pickling done.")
