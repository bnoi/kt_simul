import sys
sys.path.append('..')

import argparse

from kt_simul.gui.animation import Animator
from kt_simul.io.simuio import SimuIO
from kt_simul.core.simul_spindle import Metaphase
from kt_simul.core import parameters


if __name__ == '__main__':

    PARAMFILE = parameters.PARAMFILE
    MEASUREFILE = parameters.MEASUREFILE

    # Arguments parser
    parser = argparse.ArgumentParser(description='KtSimu Animator')
    parser.add_argument('--new', '-n', action='store_true', default=False,
                         help="Run a new simulation")
    parser.add_argument('--results', '-r', type=str, default="simu.h5",
                         help="Specified results hdf5 file for new or already runned simulation")
    args = parser.parse_args()

    results_file = args.results

    if args.new:

        meta = Metaphase(paramfile=PARAMFILE,
                         measurefile=MEASUREFILE,
                         verbose=True)
        meta.simul()
        io = SimuIO(meta)
        io.save(results_file)

    else:

        meta = SimuIO().read(results_file)

    # Run animation
    anim = Animator(meta)
    anim.play()
