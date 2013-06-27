# -*- coding: utf-8 -*-

import argparse

from kt_simul.cluster import Launcher
from kt_simul.cluster.process import Process

import config

if __name__ == '__main__':

    # Try to import local config.py variables
    params_file = config.params_file
    measures_file = config.measures_file
    results_file = config.results_file
    nsimu = config.nsimu
    results_path = config.results_path
    simu_name = config.simu_name
    evaluations = config.evaluations

    # Arguments parser
    parser = argparse.ArgumentParser(description='KtSimu Launcher')

    parser.add_argument('--nsimu', "-n", type=int, default=nsimu,
                        help='Number of simulations to launch (default = 10)')
    parser.add_argument("--results", "-r", type=str, default=results_path,
                        help='Directory to store results')
    parser.add_argument("--name", "-a", type=str, default="",
                        help='Name of the simulations')

    parser.add_argument("--eval", '-e', default=False,
                        action="store_true", help='Launch pool evaluations after simulations has runned')
    parser.add_argument("--only-eval", '-o', dest="only_eval", type=str, default="",
                        help='Launch only pool evaluations and specified ktsimu results path')

    parser.add_argument('--params', '-p', type=str, default=params_file,
                         help="For new simulation, specified params.xml file")
    parser.add_argument('--measures', '-m', type=str, default=measures_file,
                         help="For new simulation, specified measures.xml file")

    args = parser.parse_args()

    PARAMFILE = params_file
    MEASUREFILE = measures_file

    result_path = args.results
    number_simu = args.nsimu
    name = args.name
    pool_eval = args.eval
    only_eval = args.only_eval
    ncore = 4

    if name == "":
        name = "%s_n%i" % (config.simu_name, number_simu)


    if only_eval:
        p = Process(results_path=only_eval)
        resu = p.evaluate(run_all=True, debug=True)

    else:
        l = Launcher(result_path,
                     number_simu,
                     name=name,
                     ncore=ncore,
                     paramfile=PARAMFILE,
                     measurefile=MEASUREFILE)

        l.run()

        if pool_eval:
            p = Process(results_path=l.results_path)
            resu = p.evaluate(groups=evaluations, debug=True)
