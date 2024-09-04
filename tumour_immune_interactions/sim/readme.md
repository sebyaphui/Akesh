## All the files used for the simulation


Name | Usage
-- | --
`almeida_model.py` | The original discrete lattice simulation (though it should still be possible to run this simulation using `discrete_model.py`). Create a `Simulation` object `sim`, and then do `sim.run()`. It is not recommended to manually set up a simulation this way however. Refer to `discrete_model.py` instead.
`continuous_model.py` | The continuous limit of Almeida's model. Use `continuous_run.py` to set up and run this simulation.
`discrete_model.py` | The full discrete model, supporting any `PhenotypeSpace` include lattice and sequence. Use `main.py` or `run.py` to set up and run this simulation. Otherwise, initialise a `Simulation` object manually, and do `sim.run()`.
`extend.py` | For extending a simulation that has been saved. Will use _specifically_ the path specified in `conf.py`, so usually `sim_data/sim.pickle`.
`graphing.py` | For creating line graphs, histograms and fish plots of results. Standard functions will use _specifically_ the path specified in `conf.py`, so usually `sim_data/sim.pickle`.
`inputs.py` | For any I/O related to the simulation. Loading simulations and configurations, verifying, and more.
`main.py` | The default entry into the program. Run the simulation, letting you create a new one if there is not one pre-existing. You are then allowed to generate a graph, which will be placed in the default output directory in `conf.py`, so usually `outputs/sim`.
`phenotype.py` | All the phenotype-related classes. Explore this if you wish to make a custom phenotype.
`run.py` | The file that runs the discrete model. Use over `main.py` if you do not wish to create graphs immediately after the simulation is run.

## General input notes

`y` is always a yes, and any other input is regarded as a no.

## Notes on input for `discrete_model.py`

At the end, you are asked about whether you want to "append the config name?". Generally say no to this. This setting is most useful when you are creating lots of simulations and want to save each by a different filename. In particular, graph viewing won't work by default since the file will no longer be `sim.pickle`

If you just want to run a simulation one at a time, and then see the graph, indeed, say no to appending the config name.

## See `config/readme.md` for more important configuration notes