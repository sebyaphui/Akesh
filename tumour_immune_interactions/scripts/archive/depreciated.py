"""
A depreciated implementation of cells that involved storing lists of many cell objects, rather than simply storing a count. 
"""

from ..discrete_model import PhenotypeStructure, UniversalCellParams
import random


class Cell:
    """
    Depreciated.
    """

    def __init__(self, phenotype_id, phen_struct: PhenotypeStructure):
        self.phenotype_id = phenotype_id
        self.phen_struct = phen_struct

    def mutate_int(self, no_steps, direction):
        self.phenotype_id = self.phen_struct.shift(
            self.phenotype_id, no_steps, direction
        )

    def mutate_left(self):
        self.mutate_int(1, -1)

    def mutate_right(self):
        self.mutate_int(1, 1)

    @classmethod
    def random(self, phenotype_structure: PhenotypeStructure):
        return Cell(phenotype_structure.get_random_phenotype_id(), phenotype_structure)


class Cells:  # Has a 1000x scaling
    """
    Depreciated and inefficient cell storing object
    """

    def __init__(
        self,
        cells: set[Cell],
        universal_params: UniversalCellParams,
    ):
        self.cells = cells
        self.no_cells_at_phenotype = {}
        self.universal_params = universal_params

    def compute_cells_at_each_phenotype(self):
        self.no_cells_at_phenotype = {}  # Reset dictionary (inefficient, but clear)

        for cell in self.cells:
            if cell.phenotype_id in self.no_cells_at_phenotype:
                self.no_cells_at_phenotype[cell.phenotype_id] += 1
            else:
                self.no_cells_at_phenotype[cell.phenotype_id] = 1

    def get_no_cells_at_phenotype(self, phenotype_id):
        if phenotype_id in self.no_cells_at_phenotype:
            return (
                1000 * self.no_cells_at_phenotype[phenotype_id]
            )  # This should be removed
        else:
            return 0

    def __len__(self):
        return len(self.cells)

    @classmethod
    def random(
        self,
        number,
        universal_params: UniversalCellParams,
        phen_struct: PhenotypeStructure,
    ):
        # Create some number of random cells
        cells = set()
        for i in range(number):
            cell = Cell.random(phen_struct)
            cells.add(cell)
        return Cells(cells, universal_params)

    @classmethod
    def evolve_population(
        self,
        cells,
        get_phenotype_probabilities,
    ):
        """
        phenotype_probabilities has the birth, death and quiescence probabilities of the population
        """
        new_cells = set()
        dead_cells = set()
        for cell in cells.cells:
            weights = get_phenotype_probabilities(cell.phenotype_id)

            action_name = random.choices(
                population=["birth", "death", "quiescence"],
                weights=weights,
            )[0]
            """
            if action_name != "quiescence":
                print(action_name)
            """
            if action_name == "birth":
                new_cell = Cell(cell.phenotype_id, cell.phen_struct)
                new_cells.add(new_cell)
            elif action_name == "death":
                dead_cells.discard(cell)
        for cell in dead_cells:
            cells.cells.discard(cell)
        for cell in new_cells:
            if cell in cells.cells:
                print("Already here")
            cells.cells.add(cell)
