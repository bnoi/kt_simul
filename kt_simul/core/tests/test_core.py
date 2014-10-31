import numpy as np

from kt_simul.io.xml_handler import ParamTree
from kt_simul.core.simul_spindle import Metaphase
from kt_simul.core import parameters


PARAMFILE = parameters.PARAMFILE
MEASUREFILE = parameters.MEASUREFILE

# Change some parameters
measuretree = ParamTree(MEASUREFILE, adimentionalized=False)


from kt_simul.core.components import Spindle, Spb, Chromosome, Organite
from kt_simul.core.components import Centromere

class DumbKD:
    def __init__(self, paramtree):
        self.num_steps = 10
        self.params = paramtree.relative_dic
        self.prng = np.random.RandomState()
        self.initial_plug = 'random'

def test_components():

    paramtree = ParamTree(PARAMFILE)
    dKD = DumbKD(paramtree)

    spindle = Spindle(dKD)
    organite = Organite(spindle, 0, 0)
    assert organite.pos == 0
    assert organite.normal == 0

    spb = Spb(spindle, 1, 2)
    assert spb.side == 1

    chromosome = Chromosome(spindle, 0)
    assert chromosome.id == 0

    cenA = Centromere(chromosome, 'A')
    assert cenA.pos == chromosome.pos