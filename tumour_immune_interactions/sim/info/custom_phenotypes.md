## 1. How do we go about creating a custom phenotype?

First, we decide what `PhenotypeStructure`s will be used (`ps1`, `ps2`, ...)
These structures can interact with each other, in pairs. For example, we could have an interaction between `ps1` and `ps4`. Hence, we need to define _how_ this interaction occurs, and do so for every possible set of pair.
By design (though possibly inconvenient), you have to define _all_ unordered pairs, since (`ps1`,`ps2`) might require a different computation to (`ps2`,`ps1`). We do this through a `PhenotypeInteractions` object.

Each interaction consists of a function, and a piece of data/ global arguments that are fed into the function.
In particular, a function looks like this:

`f(phen_1 : Phenotype, phen_2 : Phenotype, range : float, *data)`

where `data` can be any information/ object whatsoever.

## 2. How is this done programatically?

If you wish to use pre-existing `PhenotypeStructure` classes, then it only suffices to follow this section (and there are still lots of possibilities for interactions!). If you would like to create a _new_ class, see Section 3 below.

We set up our phenotypes before the `Simulation` starts, in its `__init__` function.

Start by defining all the structures you would like to use:
```py
self.struct_A = SomePhenotypeStructure(*args_A)
self.struct_B = SomeOtherPhenotypeStructure(*args_B)
```
(as pseudocode)

These structures should be assigned to `CellBundle`s later in this initialisation.

Now you can create your phenotype interactions object
```py
self.phen_int = PhenotypeInteractions()
```
For each pair, choose or define your function and data. We set these interactions up as follows:
```py
self.phen_int.set_up_interactions(
    (struct_A, struct_A),
    data_AA,
    f_AA
)
...
```
(Rinse and repeat)

Where should this all be placed in the code? Strictly, you can just replace the if statement, for the default subtype, `lattice`. If you want to do this properly, however, see Section 4.

## 3. Creating a new `PhenotypeStructure` class
You must inherit the `ABC`. 
It may be useful to create class methods detailing "generic data" or an "interaction function" that this structure will use when interacting with itself. 

Must be hashable, and comparable (with `__eq__`). This can simply be done by hashing or comparing the arguments resepectively.

## 4. Creating a new `subtype`
For my simulation, the `simtype` is always `discrete`.

To create a `subtype`, you just have to specify in `arguments.csv` what parameters are required, optional and ignored (`0`, `1`, `2`).

To implement this subtype into the simulation, find the if statement that checks the subtype. Here you want to:
1) Load the relevant keyword arguments (e.g. `kwargs["my_parameter"]`)
2) Run the phenotype setup outlined in Section 2


