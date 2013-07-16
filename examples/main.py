from __future__ import unicode_literals
from __future__ import division

import argparse
import sys
sys.path.append('..')

from kt_simul.io.xml_handler import ParamTree
from kt_simul.core.simul_spindle import Metaphase
from kt_simul.io.simuio import SimuIO
from kt_simul.core import parameters

if __name__ == '__main__':

    PARAMFILE = parameters.PARAMFILE
    MEASUREFILE = parameters.MEASUREFILE

    # Change some parameters
    paramtree = ParamTree(PARAMFILE)
    paramtree.change_dic('dt', 10)
    paramtree.change_dic('span', 2000)
    paramtree.change_dic('t_A', 1750)

    measuretree = ParamTree(MEASUREFILE, adimentionalized=False)

    # Arguments parser
    parser = argparse.ArgumentParser(description='KtSimu Launcher')
    parser.add_argument('--new', '-n', action='store_true', default=False,
                         help="Run a new simulation")
    parser.add_argument('--results', '-r', type=str, default="simu.h5",
                         help="Specified results hdf5 file for new or already runned simulation")
    args = parser.parse_args()

    results_file = args.results

    if args.new:

        # Init simu
        meta = Metaphase(verbose=True, paramtree=paramtree, measuretree=measuretree, initial_plug='random')

        # Launch simu
        meta.simul()

        # Save results
        SimuIO(meta).save(results_file)

        # Show trajectories (matplotlib needed)
        meta.show()

    else:

        meta = SimuIO().read(results_file)
        meta.show()
