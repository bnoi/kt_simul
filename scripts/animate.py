import argparse

from kt_simul.gui.animation import Animator
from kt_simul.io.simuio import SimuIO
from kt_simul.core.simul_spindle import Metaphase

import config


if __name__ == '__main__':

    # Try to import local config.py variables
    params_file = config.params_file
    measures_file = config.measures_file
    results_file = config.results_file

    # Arguments parser
    parser = argparse.ArgumentParser(description='KtSimu Animator')
    parser.add_argument('--new', '-n', action='store_true', default=False,
                         help="Run a new simulation")
    parser.add_argument('--params', '-p', type=str, default=params_file,
                         help="For new simulation, specified params.xml file")
    parser.add_argument('--measures', '-m', type=str, default=measures_file,
                         help="For new simulation, specified measures.xml file")
    parser.add_argument('--results', '-r', type=str, default=results_file,
                         help="Specified results Zip file for new or already runned simulation")
    args = parser.parse_args()

    results_file = args.results

    if args.new:

        params_file = args.params
        measures_file = args.measures

        meta = Metaphase(paramfile=params_file,
                         measurefile=measures_file,
                         verbose=True)
        meta.simul()
        io = SimuIO(meta)
        io.save(results_file)

    else:

        meta = SimuIO().read(results_file)

    # Run animation
    anim = Animator(meta)
    anim.play()
