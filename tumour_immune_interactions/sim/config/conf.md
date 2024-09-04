## The most important configuration file


Item | Effect
--- | ---
interrupt | Set to true to safely interrupt at the next time step. This will save the simulation as a `.pickle`, and end the execution. (The file will also be reset so `interrupt = False` once more)
debug | Set to true to see output triggered by `self.print` inside `Simulation`
path_to_data | Where you read from and save to for simulation `sim.pickle` files
path_to_output | Where graph outputs and more go.
sim_state_init_type | `detailed` or ==?==. Specifies how much data is stored inside each `SimulationState`. `detailed` will be more performance heavy, but is certainly necessary for e.g. fish plotting