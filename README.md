# Pazy Wing Model for SHARPy

This repository contains the [SHARPy](http://github.com/imperialcollegelondon/sharpy) version of the equivalent beam
model of the Pazy wing described in [1]. The model
has been slightly modified to model the mass and inertia in a distributed fashion per element rather than lumped at the 
finite-element nodes.

## Installation

Clone the repository to your local computer. It is intended for use within SHARPy so you must ensure that SHARPy is 
properly installed in your system, including the conda environment provided.

With SHARPy and its environment installed, the only required step to use this model is to run from the terminal where
you are running your scripts

```bash
source <path-to-pazy-model>/bin/pazy_vars.sh
```
 
This will append the location of the `pazy-model` folder to the system's path, such that your code will be able to find
the modules.

## Using the model

In your SHARPy model generation script, import the Pazy wing:

```python
from pazy_wing_model import PazyWing
```

The class `PazyWing` is all you need as interface, and it is documented in its source file.

The initialisation call takes the form of

```python
wing = PazyWing(case_name, case_route, pazy_settings)
```

where the `pazy_settings` refer to whether the skin should be modelled, discretisation refinement etc. A list of these
settings is provided in the source file.

For a purely structural simulation (no aerodynamic grid generated):

```python
wing.generate_structure()
wing.save_files()
```

That will generate the `.fem.h5` file that is used by SHARPy.

The `PazyWing` class has an attribute `config` that contains the ConfigObj with the simulation settings. You can
complete it with the desired simulation parameters. When you are done:

```python
wing.config.write()
```

That will save the `.sharpy` file containing the simulation settings. You can now run SHARPy!

## References

[1] Riso, C., & Cesnik, C. E. S.. Equivalent Beam Distributions of the Pazy Wing. University of Michigan, 2020.

## Contact

If you have any questions and want to get in touch, 
[contact us](https://www.imperial.ac.uk/aeroelastics/people/goizueta/).

If you have any questions on how to use the model or find any bugs please file an issue. 