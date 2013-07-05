from kt_simul.io.xml_handler import ParamTree
from kt_simul.core.simul_spindle import Metaphase
from kt_simul.io.simuio import SimuIO
from kt_simul.core import parameters

PARAMFILE = parameters.PARAMFILE
MEASUREFILE = parameters.MEASUREFILE

# Change some parameters
paramtree = ParamTree(PARAMFILE)
paramtree.change_dic('dt', 10)
paramtree.change_dic('span', 2000)
paramtree.change_dic('t_A', 1750)

measuretree = ParamTree(MEASUREFILE, adimentionalized=False)

# Init simu
meta = Metaphase(verbose=True, paramtree=paramtree, measuretree=measuretree, initial_plug='random')

# Launch simu
meta.simul()

# Save results
SimuIO(meta).save("simu.h5")

# Show trajectories (matplotlib needed)
meta.show()
