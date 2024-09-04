# A set of simulations modelling the behaviour and interaction of T-cells and tumour cells in the body.


## Setting up

Begin by installing the `tumour_immune_interactions` package. 

```py
python -m pip install -e .
```

You will then for example be able to do:

```py
from sim.discrete_model import Simulation
```

## Configuring

### First create or find the configuration you would like within `configurations.csv`
See `sim/config/readme.md`.

### The primary simulation configuration file `conf.py`
This is used for:
1) Interrupting the simulation
2) Specifying where outputs and simulation files go
3) Changing how much simulation data is stored

See `sim/config/conf.md` for details.


### Caveats

You can specify a `sequence_matrix_config` and `Simulation` will load this in, but the simulation logic does not currently support using it (currently it remains unused and an "identity" matrix is used instead).

The code, and input files refer to an `affinity_matrix`. This is a relic that should be removed, because this matrix is utilised as a `binding_probability_matrix` instead. Renaming the codebase can be done, but would require some time.

## Running a simulation

Running `main.py` will run the discrete simulation. 
Running `continuous_run.py` will run the continuous simulation (not developed by me, but adapted by me)


# What simulations can be run, and what features exist?

Three types are fully implemented:
- Discrete lattice simulations
- Discrete sequence simulations
- Continuous "lattice" simulations

The discrete sequence simulation was the type made for this project.

Functionality
- A configuration manager
- Argument validation
- The ability to interrupt simulations, save and load them
- Flexible storage of simulation history
- The ability to implement new `PhenotypeStructure` types
- A graphing suite which supports line graphs, histograms and fish plots


Note that we refer to T-cells as `CTL`s in the code.