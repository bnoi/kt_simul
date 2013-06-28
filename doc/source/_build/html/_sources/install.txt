How to install kt_simul
=======================

Dependencies
------------

Package manager versions of the python libraries should work.
 
* Python >= 2.5
* Numpy >= 1.4 and Scipy >= 0.9
* Cython >= 0.14
* Qt4 and PySide (optional)

Installation
------------

You should first clone the github version of this code, then
use the setup script, whether via:

  # python setup.py build

for a local install, or
  
  # python setup.py install

for a system install.
You will need a C compiler for the cython part.

Usage
-----
 
>>> import kt_simul.simul_spindle as sim
>>> metaph = sim.Metaphase()
>>> metaph.simul()
>>> metaph.show_one()
>>> help(metaph)

Should provide a first view of the simulation.
Using the GUI:

    # python [path_to_package]/kt_smul/gui/kineto_simulation.py

Further details can be found by looking at the code.
