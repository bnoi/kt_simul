# `kt_simul` - Python model of chromosome mouvements during mitosis
[![Build Status](https://travis-ci.org/bnoi/kt_simul.svg)](https://travis-ci.org/bnoi/kt_simul)

Python model of chromosome mouvements during mitosis in Fission Yeast (http://damcb.org).

This python module provides the simulation of mitotic spindle elements (for now,
the kinetochore and the spindle pole bodies), during cell division in fission
yeast. For more details see [online documentation](http://bnoi.github.io/kt_simul/ "kt_simul documentation").

If you're interested in this work, please contact :

- Sylvie Tournier (sylvie.tournier-gachet@univ-tlse3.fr) for
  all aspects related to the biology and for academical collaboration purposes.
- Guillaume Gay (guillaume@damcb.com) for all aspects related to the code

The underlying model is fully described in:

G. Gay, T.Courthéoux, C. Reyes, S. Tournier, Y. Gachet. *A stochastic model of
kinetochore–microtubule attachment accurately describes fission yeast chromosome
segregation* J. Cell Biol 2012 [(link to article)](http://jcb.rupress.org/content/196/6/757.abstract).

This article should be used for any citation of this work.

Funding
-------

This project is funded by the French National Research Agency as:
   *ANR- BLAN 1206 01 Chromocatch, programme blanc 2010*

Contributors
------------

- Guillaume Gay <gllm.gay@gmail.com> : main author
- Hadrien Mary <hadrien.mary@gmail.com> : contributor

`kt_simul` is part of the BNOI Project <https://github.com/bnoi>.


Licence
-------

This code is provided under the GPL compatible CeCILL licence (see
LICENCE for full details).

Dependencies
------------

- Python 3.4
- numpy
- scipy
- pandas
- matplotlib
- PyQt4 (optional)
- vispy (optional)

Installation
------------

You can install `kt_simul` via pip :

    pip install -e git+https://github.com/bnoi/kt_simul.git#egg=master

Usage
-----

```python
from kt_simul.core.simul_spindle import Metaphase
from kt_simul.io.simuio import SimuIO
from kt_simul.core import parameters

paramtree = parameters.get_default_paramtree()
paramtree['dt'] = 10
paramtree['span'] = 2000
paramtree['t_A'] = 1750

measuretree = parameters.get_default_measuretree()
measuretree['mean_metaph_k_dist'] = 0.3  # 0.3

# Init simu
meta = Metaphase(verbose=True,
                 paramtree=paramtree,
                 measuretree=measuretree,
                 initial_plug='random',
                 keep_same_random_seed=False,
                 force_parameters=[])

# Launch simu
meta.simul()

# Save results
SimuIO(meta).save("simu.h5")

# Show trajectories (matplotlib needed)
fig = meta.show()
```

![Chromosomes trajectories](examples/trajectories.png "Chromosomes trajectories")

Simulation can be "played" with Qt based GUI:

```python
from kt_simul.gui.animation import Animator

anim = Animator(meta)
anim.play()
```

![GUI](examples/gui.gif "GUI")
