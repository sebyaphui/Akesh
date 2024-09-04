### How do you add a new parameter?
- Add it to `configurations.csv`
- Make it a _necessary_ or _optional_ argument for your desired `subtype`/ `simtype` (though I assume the latter is likely `discrete`)
- See `run.py` and add your parameter, so it gets loaded under the relevant `subtype` if statement. e.g. `yourparam = cf.yourparam`
- Also in `run.py`, add it as a keyword argument to the `Simulation` constructor
- Inside the `__init__` method of `Simulation`, read out the keyword argument, under the loading for the relevant `subtype`