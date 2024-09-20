# Solar Sim

A solar system orbit simulator based on Kepler and Newton's orbital geometry and mechanics.

![GitHub License](https://img.shields.io/github/license/RandomKiddo/SolarSim)


___

### Command-Line Arguments

During runtime, run the following to see the most up-to-date arguments:
```shell
python3 sim.py -h
```

The following output is a direct pasting of the above shell command:

usage: SolarSim [-h] [--dt DT] [-u UNITS] [-s] [-p PAUSE] [--fp FP] [--sp SP] [-cut CUTOFF] [-g GIF] [--pt PT] [-v]

A python-based celestial body system simulator

optional arguments: <br>
  `-h, --help`           show this help message and exit<br>
  `--dt DT`               differential time step to use in the simulation<br>
  `-u UNITS, --units UNITS`
                        system of units to use for the simulation<br>
  `-s, --solar`           use the solar system for the simulation<br>
  `--fp FP `              the filepath to a json file representing the system<br>
  `--sp SP`             the path to save the simulation paths as a csv file<br>
  `-cut CUTOFF, --cutoff CUTOFF`
                        the zoom factor to zoom-in by<br>
  `-g GIF, --gif GIF`     store the simulation as a gif with the provided filepath<br>
  `--pt PT`               plot pause timing to use in the simulation<br>
  `-v, --verbose `        log more information to the console<br>

Functionality warnings (This does not include simulation warnings like "large dt" or "fp not found":
* Gif functionality not currently integrated-- will be added in next update.
* 'Save-path functionality working, but reading from csv does not save all metadata.

___

### Running a Simulation

Running a simulation can be done in two ways: (i) using the command line, or (ii) hard-coding the simulation in python.

I. Using the Command Line

Navigate to the directory the `sim.py` file is located. Ensure that its dependencies are installed or located in the same directory, namely the `cgen.py` file.

```shell
python3 sim.py [args]
# for arguments, see above section
```

If not using the solar system and instead using a user-defined system, make sure to include the `--fp FP` argument to the JSON file with the system information. Use the `solar_system.json` file as reference on how to write such files or if you have the python `List[CelestialBody]` instance, use the `json_creator.py` file.

If using the solar system, be sure to include the `solar_system.json` file in the same directory as the `sim.py` file.

II. Hard-Coding in Python

If hard coding in python, download all `.py` files and store in a findable place in order to import them. If the solar system model is going to be use, download the `solar_system.json` file as well.

```py
"""
Some function and constructor arguments have been left out.
This has been done for ease of reading this section.
To see all function and constructor arguments, inspect the source code files.
"""

from cgen import *
from sim import *

system = 'au' # AU system of units

# Using the solar system
solar_system = System.solar_system()
simulation = Simulation(c_system=solar_system, system=system)
simulation.sim()

# Creating a custom system
bodies = [
  CelestialBody(
    m=Constants.solar_mass(system=system), 
    r=Constants.earth_radius(system=system)
  ),
  CelestialBody(
    m=Constants.earth_mass(system=system),
    r=5e-5,
    a=1.0
  )
]
c_system = System(bodies=bodies)
simulation = Simulation(c_system=c_system, system=system)
simulation.sim()
```

___

### Simulation Warnings

This simulation emits warnings for the following run-time events:

* Differential Time Step Warning: If the `dt` parameter is relatively large (`>0.01`), a warning is raised saying that the simulation could yield lossy results.
* File Path Warning: If in the command line the solar system argument is not provided and a filepath is not provided or the solar system argument is not provided and the given filepath could not be found, then a warning is raised stating that the solar system will be used against the arguments from the command line.
* Functionality Warning: If using the `-g` or `--sp` arguments in the command line, this warning is raised saying that the functionality of these arguments aren't completed or aren't functioning completely as-intended.
* Unit Warning: The provided system of units is not hard-coded and therefore not recognized. Raised when the given system of units is not one of `'mks'`, `'au'` or `'cgs'`.
* No-Body System Warning: When initializing a `System` instance, the `bodies` argument is `None`, raising a warning that the solar system would be used by default.

___

[Back to Top](#solar-sim)

<sub>This page was last edited on 09.19.2024</sub>