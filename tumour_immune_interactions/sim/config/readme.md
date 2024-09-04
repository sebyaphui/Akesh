# All the configuration and data files needed for all the simulations

Name | Use 
-- | --
`configurations.csv` | For tracking every possible set of parameters. The simulation reads a particular row from this file, when given a certain `name` to get all the parameters for the simulation it is running.
`arguments.csv ` | Here we define the `simtype` and `subtype` of a simulation that can run. This file specifies for each argument, whether it is required (`1`), optional (`2`) or ignored (`3`). Used for verification purposes. Aligned with the headers of `configurations.csv` (done manually).
`conf.py` | Small pythonic configuration variables, often used for debugging purposes (or anything that would be more difficult/ unintuitive than justified to store in a text file/ as a `configurations.csv` parameter).
`peptides_from_TULIP.txt` | A set of possible tumour cell sequences, extracted from TULIP.
`TCRs_from_TULIP.txt` | A set of possible CTL cell seqeuences, extracted from TULIP.

# See `configurations.md` for information on how to set up a new simulation i.e. a new row in `configurations.csv`
# See `arguments.md` for which parameters are optional, required or unused for each sequence type
# See `../info/sequence_phenotype_structure.md` for details on what types of matrices are allowed