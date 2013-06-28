kt_simul
========

Python model of chromosome mouvements during mitosis in Fission Yeast http://mitotic-machine.org.

This python module provides the simulation of mitotic spindle elements (for now, the
kinetochore and the spindle pole bodies), during cell division in
fission yeast.

If you're interested in this work, please contact :

- Sylvie Tournier (sylvie.tournier-gachet@univ-tlse3.fr) for
  all aspects related to the biology and for academical collaboration purposes.
- Guillaume Gay (gllm.gay@gmail.com) for all aspects related to the code

The underlying model is fully described in:

G. Gay, T.Courthéoux, C. Reyes, S. Tournier, Y. Gachet. *A stochastic model of kinetochore–microtubule attachment
accurately describes fission yeast chromosome segregation* J. Cell Biol 2012 ``doi: 10.1083/jcb.201107124``

This article should be used for any citation of this work.

FUNDING
-------

This project is funded by the French National Research Agency as:
   *ANR- BLAN 1206 01 Chromocatch, programme blanc 2010*

LICENCE
-------

This code is provided under the GPL compatible CeCILL licence (see
LICENCE for full details).

DEPENDENCIES
------------

Package manager versions of the python libraries should work.

- Python >= 2.5
- Numpy >= 1.4 and Scipy >= 0.9
- Cython >= 0.14
- Qt4 and PySide (optional)

INSTALLATION
------------

You should first clone the github version of this code, then
use the setup script, whether via:

```
python setup.py build
```

for a local install, or

```
python setup.py install
```

for a system install.
You will need a C compiler for the cython part.

USAGE:
-----

```python
import kt_simul.simul_spindle as sim
metaph = sim.Metaphase()
metaph.simul()
metaph.show_one()
help(metaph)
```

Should provide a first view of the simulation.
Using the GUI:

```
python [path_to_package]/kt_smul/gui/kineto_simulation.py
```

Further details can be found by looking at the code.

