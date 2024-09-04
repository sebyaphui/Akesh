"""
The entry point and interface for the program.
"""

from sim.run import run, extend
from sim.inputs import set_up_and_get_arguments
from cProfile import Profile
from pstats import SortKey, Stats
from sim.graphing import graph_from_path, savefig


save_fig, overwrite, config_name = set_up_and_get_arguments()

# extend(100)
run(overwrite, config_name)


def profile():
    with Profile() as profile:
        run(overwrite, config_name)
        Stats(profile).strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats()


"""
profile()
"""

if save_fig is None:
    save_fig = input("Would you like to save the graph?")

# This naming _is_ confusing, but alas
if save_fig == "y":
    savefig()
