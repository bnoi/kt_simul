from __future__ import unicode_literals
from __future__ import division

import argparse
import sys

sys.path.append('..')

from kt_simul.io.parameters_io import ParamTree
from kt_simul.core.simul_spindle import Metaphase
from kt_simul.io.simuio import SimuIO
from kt_simul.core import parameters

if __name__ == '__main__':

    PARAMFILE = parameters.PARAMFILE
    MEASUREFILE = parameters.MEASUREFILE

    # Change some parameters
    param_tree = ParamTree(PARAMFILE)
    # param_tree["dt"] = 10
    # param_tree["span"] = 2000
    # param_tree["t_A"] = 17500
    # print(param_tree.params)

    measure_tree = ParamTree(MEASUREFILE, adimentionalized=False)
    # print(measure_tree.params)

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
        meta = Metaphase(verbose=True, paramtree=param_tree, measuretree=measure_tree, initial_plug='random')

        # Launch simu
        meta.simul()

        # Save results
        SimuIO(meta).save(results_file)

        # Show trajectories (matplotlib needed)
        meta.show().savefig("trajectories.png")

    else:

        meta = SimuIO().read(results_file)
        meta.show()
