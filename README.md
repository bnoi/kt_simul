# `kt_simul` - Python model of chromosome mouvements during mitosis
[![Build Status](https://travis-ci.org/bnoi/kt_simul.svg)](https://travis-ci.org/bnoi/kt_simul)

Python model of chromosome mouvements during mitosis in Fission Yeast.

This python module provides the simulation of mitotic spindle elements, during cell division.

If you're interested in this work, please contact :

- __Sylvie Tournier__ (sylvie.tournier-gachet@univ-tlse3.fr) for all aspects related to the biology and for academical collaboration purposes.
- __Guillaume Gay__ (guillaume@damcb.com and http://damcb.org) for all aspects related to the physics and the code.

The underlying model is fully described in:

[G. Gay, T.Courthéoux, C. Reyes, S. Tournier, Y. Gachet. A stochastic model of kinetochore–microtubule attachment accurately describes fission yeast chromosome segregation J. Cell Biol 2012](http://jcb.rupress.org/content/196/6/757.abstract)

This article should be used for any citation of this work.

The 1D version has been used to describe chromosome congression in the following article : [Mary, H., Fouchard, J., Gay, G., Reyes, C., Gauthier, T., Gruget, C., Pecreaux, J., Tournier, S. and Gachet, Y. (2015). Fission yeast kinesin-8 controls chromosome congression independently of oscillations. J. Cell Sci](http://jcs.biologists.org/content/128/20/3720)

The original version in 1D is available in the git branch `simu_1D`. The main branch `master` uses a new geometry system to simulate chromosomes in 3D.

Some code example and documentation con be found here : http://nbviewer.ipython.org/github/bnoi/kt_simul/tree/master/doc/notebooks/.

Mathematical documentation about the 1D model is available at [doc/simu_spindle_1d.pdf](doc/simu_spindle_1d.pdf).

## Installation

The best way to install the dependencies is to use the [Anaconda Python distribution](https://www.continuum.io/downloads). Then use the following command lines :

```bash
conda config --add channels conda-forge
conda install numpy scipy pandas matplotlib
```

To install the `kt_simul` model, do :

```bash
pip install -e git+https://github.com/bnoi/kt_simul.git#egg=master
```

## Funding

This project is funded by the French National Research Agency as : *ANR- BLAN 1206 01 Chromocatch, programme blanc 2010*

## Contributors

- Guillaume Gay <guillaume@damcb.com> : main author
- Hadrien Mary <hadrien.mary@gmail.com> : contributor

## Licence

This code is provided under the GPL compatible CeCILL licence (see
[LICENCE](LICENSE) for full details).
