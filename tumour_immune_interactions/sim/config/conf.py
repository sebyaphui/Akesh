from sim.inputs import get_file_dir


interrupt = False
debug = True
path_to_data = get_file_dir() + "/../sim_data/sim.pickle"
path_to_output = get_file_dir() + "/../outputs/sim"
m_adjustment = True  # Make my simulation seem more like Marta's (But this can likely be removed, since it's just for small temporary changes)
sim_state_init_type = "detailed"


def hello_world(sim):
    print("hello")
