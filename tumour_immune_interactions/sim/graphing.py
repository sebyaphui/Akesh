"""
Graphing functionality, including line graphs, histograms and fish plots. Used inside the simulation, and can also be utilised externally.
"""

from sim.discrete_model import Simulation, SimulationState, CellBundle
from sim.config.conf import path_to_data, path_to_output
import numpy as np
import pandas as pd
from typing import Callable
import os

import warnings


def get_sim(path=path_to_data) -> Simulation:
    try:
        sim = Simulation.load_simulation(path)
        return sim
    except IOError:
        print("Simulation not found.")
        return


def get_pops(sim: Simulation):
    tumour_cells_pop = []
    CTL_cells_pop = []
    for state in sim.history.history:
        tumour_cells_pop.append(state.tumour_cells_pop)
        CTL_cells_pop.append(state.CTL_cells_pop)
    return tumour_cells_pop, CTL_cells_pop


def graph(sim: Simulation, CTL_label="CTL Cells", tumour_label="Tumour Cells"):
    import matplotlib.pyplot as plt

    tumour_cells_pop, CTL_cells_pop = get_pops(sim)
    times = np.linspace(0, sim.time_step * sim.time_step_size, sim.time_step)
    plt.plot(times, CTL_cells_pop, label=CTL_label)
    plt.plot(times, tumour_cells_pop, label=tumour_label)
    plt.legend()
    return plt


def report_graph(sim: Simulation):
    import matplotlib.pyplot as plt

    plt.figure()
    plt.rcParams["figure.figsize"] = (10, 8)
    plt.rcParams["font.size"] = 25
    plt.rcParams["lines.linewidth"] = 4
    plt.rcParams["xtick.labelbottom"] = True
    plt = graph(sim, "T-cells", "Tumour cells")
    plt.xlabel("Time (days)", labelpad=17)
    plt.ylabel("Population Size", labelpad=17)
    return plt


def graph_tumour(sim: Simulation):
    import matplotlib.pyplot as plt

    tumour_cells_pop, CTL_cells_pop = get_pops(sim)
    times = np.linspace(0, sim.time_step * sim.time_step_size, sim.time_step)
    # plt.plot(times, CTL_cells_pop, label="CTL Cells")
    plt.plot(times, tumour_cells_pop, label="Tumour Cells")
    plt.legend()
    return plt


def fish_tumour(sim: Simulation):
    return fish(sim, "tumour")


def fish_CTL(sim: Simulation):
    return fish(sim, "CTL")


def fish(sim: Simulation, bundle_name, absolute=False):
    from pyfish import fish_plot, process_data, setup_figure
    import matplotlib.pyplot as plt

    warnings.simplefilter(action="ignore", category=FutureWarning)

    if sim.history.state_init_type == "detailed":
        bundles = [get_bundle(state, bundle_name) for state in sim.history]
        populations = []
        phenotypes = set()
        for step in range(len(bundles)):
            bundle = bundles[step]
            for phenotype, pop in bundle.cells_at_phenotype.items():
                populations.append([phenotype.id, step, pop])
                phenotypes.add(phenotype.id)
        populations_df = pd.DataFrame(
            np.array(populations), columns=["Id", "Step", "Pop"]
        )
        parent_tree = []
        """
        We need to construct the parent tree. I'm not sure quite how this works with mutations, since the parent could be any from the left/ right.
        This means no phenotypes have no parents.. And hmm, one phenotype can have multiple phenotypes.
        To generate this, I'd want to collect all present phenotypes, then order chronologically, and assign parents in descending order.
        """
        ordered_phenotypes = list(phenotypes)
        ordered_phenotypes.sort()
        parent = None
        child = None
        for i in range(len(ordered_phenotypes)):
            child = ordered_phenotypes[i]
            if parent is not None:
                parent_tree.append([parent, child])
            parent = child
        parent_tree_df = pd.DataFrame(
            np.array(parent_tree), columns=["ParentId", "ChildId"]
        )
        data = process_data(populations_df, parent_tree_df, absolute=absolute)
        sequences = bundles[-1].phen_struct.sequences
        # dist = max([len(seq) for seq in sequences])
        setup_figure(width=1200, height=1000 * len(sequences) / 12)
        times = np.linspace(0, sim.time_step * sim.time_step_size, sim.time_step)
        data = list(data)
        data[1] = times
        data[0].columns = times

        fish_plot(*data)
        return plt


def report_fish(sim: Simulation, bundle_name, absolute=False):
    import matplotlib.pyplot as plt

    plt.rcParams["font.size"] = 25
    plt.rcParams["xtick.labelbottom"] = True

    plot = fish(sim, bundle_name, absolute)
    plot.ylim([0.5, 1])
    plot.yticks([])
    plot.ylabel("")
    plot.xlabel("Time (days)")

    state = next(sim.history.__iter__())
    bundle = get_bundle(state, bundle_name)
    sequences = bundle.phen_struct.sequences

    dist = max([len(seq) for seq in sequences])

    initdtop = (1 - 0.965) * 12 / len(sequences)
    top = 1 - initdtop
    left = -dist * (sim.time_step * sim.time_step_size / 30) * 1.2
    dtop = 0.502 / len(sequences)

    for seq in sequences:
        plot.text(left, top, seq)
        top -= dtop

    return plt


def report_fish_tumour(sim: Simulation):
    return report_fish(sim, "tumour")


def report_fish_tumour_absolute(sim: Simulation):
    return report_fish(sim, "tumour", True)


def report_fish_CTL(sim: Simulation):
    return report_fish(sim, "CTL")


def get_bundle(state: SimulationState, bundle_name) -> CellBundle:
    if bundle_name == "tumour":
        return state.tumour_cells
    elif bundle_name == "CTL":
        return state.CTL_cells
    else:
        return None


def graph_from_path(path=path_to_data):
    sim = get_sim(path)
    graph(sim).show()


def flatten_dict(dict: dict):
    return [key for key, val in dict.items() for _ in range(val)]


def to_ids(phens):
    return [phen.id for phen in phens]


def hist_base(tumour_cells: CellBundle, CTL_cells: CellBundle, phen_struct):
    import matplotlib.pyplot as plt

    tumour_cell_phenotypes = to_ids(flatten_dict(tumour_cells.cells_at_phenotype))
    CTL_cell_phenotypes = to_ids(flatten_dict(CTL_cells.cells_at_phenotype))

    plt.hist(
        CTL_cell_phenotypes,
        label="CTL Cells",
        bins=phen_struct.no_possible_values,
        alpha=0.6,
    )
    plt.hist(
        tumour_cell_phenotypes,
        label="Tumour Cells",
        bins=phen_struct.no_possible_values,
        alpha=0.6,
    )
    plt.legend()
    return plt


def hist_from_state(state: SimulationState):
    return hist_base(
        state.tumour_cells, state.CTL_cells, state.tumour_cells.phen_struct
    )


def hist(sim: Simulation):
    import matplotlib.pyplot as plt

    tumour_cell_phenotypes = flatten_dict(sim.tumour_cells.cells_at_phenotype)
    CTL_cell_phenotypes = flatten_dict(sim.CTL_cells.cells_at_phenotype)

    plt.hist(
        CTL_cell_phenotypes,
        label="CTL Cells",
        bins=sim.phen_struct.no_possible_values,
        alpha=0.6,
    )
    plt.hist(
        tumour_cell_phenotypes,
        label="Tumour Cells",
        bins=sim.phen_struct.no_possible_values,
        alpha=0.6,
    )
    plt.legend()
    return plt


def report_hist_general(data, log=False):
    import matplotlib.pyplot as plt

    plt.rcParams["figure.figsize"] = (10, 8)
    plt.rcParams["font.size"] = 25
    plt.rcParams["lines.linewidth"] = 4
    plt.rcParams["xtick.labelbottom"] = True

    plt.hist(data, log=log)
    plt.xlabel("", labelpad=20)
    plt.ylabel("", labelpad=20)

    return plt


plt_fn_label = {
    graph: "graph",
    hist: "hist",
    fish_tumour: "fish_tumour",
    fish_CTL: "fish_CTL",
    report_graph: "graph_report",
    report_fish_CTL: "fish_CTL_report",
    report_fish_tumour: "fish_tumour_report",
    report_fish_tumour_absolute: "fish_tumour_absolute_report",
}


def savefig(
    sim: Simulation = None,
    plt_fn: Callable = report_graph,
    path_to_output=path_to_output,
):
    if sim is None:
        sim = get_sim()
    plt = plt_fn(sim)
    out_path = f"{path_to_output}output_{sim.config_name}_{plt_fn_label[plt_fn]}.png"
    # I'm not going to worry about removing the inputs here, because this should never have a file in it, because it's a completely empty directory!
    while os.path.exists(out_path):
        overwrite = input("Would you like to overwrite the existing file?")
        if overwrite == "y":
            break
        else:
            still_save = input("Would you still like to save the file?")
            if still_save == "y":
                custom = input("Would you like to enter a fully custom path?")
                if custom == "y":
                    out_path = input("Enter the new path:")
                    continue
                else:
                    version = input("Enter which version this figure is:")
                    out_path = f"{path_to_output}output_{sim.config_name}_{plt_fn_label[plt_fn]}_{version}.png"
                    continue
            else:
                print("The plot has not been saved.")
                return
    plt.savefig(out_path, bbox_inches="tight")


def savefig_unsafe(
    sim: Simulation = None,
    plt_fn: Callable = report_graph,
    path_to_output=path_to_output,
    extra="",
):
    """
    Save (without any safety) the file in the designated location. This will overwrite all previous files.
    """
    if sim is None:
        sim = get_sim()
    plt = plt_fn(sim)
    out_path = f"{path_to_output}output_{sim.config_name}_{plt_fn_label[plt_fn]}.png"
    if extra != "":
        out_path = f"{path_to_output}output_{sim.config_name}_{plt_fn_label[plt_fn]}_{extra}.png"
    plt.savefig(out_path, bbox_inches="tight")


def hist_from_path(path=path_to_data):
    sim = get_sim(path)
    hist(sim).show()
