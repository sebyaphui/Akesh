## Configurations
This section describes what parameters need to be specified to run the various simulations.

The configuration loaded is dependent on:
- `simtype` (this will come from the file you run: `continuous_run.py` or `run.py` i.e. whether you run Marta's model, or the discrete model)
`discrete`
`continuous`

defaults to `discrete`

- `subtype` (specified in the config file, and defaulting to `lattice`)
`lattice`
`sequence`

`name` must be unique.


All arguments are either:
- required (0)
- optional (1)
- ignored (2)

In `arguments.md`, we list (excluding `simtype`, `subtype` and `name`) all the arguments in each category for all the currently supported simulations.


## Matrix configuration
This will be done as a json file. 

### Path
- "from" : `"path"`
- "path" : the actual path from the source directory (`sim`)
- "delimiter" : If empty, "" or unspecified, then we assume a space " " delimiting

Internally, this will load the matrix, have it stored as the output of a function, and pass this function into the simulation.

### Function
- "from" : `"function"`
- "where" : if empty "", then we assume globally in the source file (`Simulation`); otherwise we import that relative reference
- "function" : the name of the function (cannot be inside a class). This function takes in the `Simulation` object after all initialisations have been done (so you may need to inspect the init file), and returns a matrix 