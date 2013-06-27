# -*- coding: utf-8 -*-

import argparse

import numpy as np

from kt_simul.cluster import Explorator
from kt_simul.cluster.process_exploration import ProcessExploration


import config

def main():

    # Try to import local config.py variables
    params_file = config.params_file
    measures_file = config.measures_file
    results_file = config.results_file
    nsimu = config.nsimu
    results_path = config.results_path
    explo_simu_name = config.explo_simu_name
    evaluations = config.evaluations

    # Arguments parser
    parser = argparse.ArgumentParser(description='KtSimu Explorator')

    parser.add_argument('--nsimu', "-n", type=int, default=nsimu,
                        help='Number of simulations to launch (default = 10)')
    parser.add_argument("--results", "-r", type=str, default=results_path,
                        help='Directory to store results')
    parser.add_argument("--name", "-a", type=str, default="",
                        help='Name of the simulations')

    parser.add_argument("--eval", '-e', default=True,
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

    # Here we explore ld_slope parameter
    parameter_to_explore = {'name': 'ldep', 'vector': np.arange(0, 1, 0.1)}

    if name == "":
        name = "%s_n%i_p%i_%s" % (config.explo_simu_name, number_simu,
                                  len(parameter_to_explore['vector']),
                                  parameter_to_explore['name'])

    if only_eval:
        p = ProcessExploration(results_path=only_eval)
        resu = p.evaluate(run_all=True, debug=True)

    else:

        explorator = Explorator(result_path,
                                number_simu,
                                name=name,
                                ncore=ncore,
                                paramfile=PARAMFILE,
                                measurefile=MEASUREFILE,
                                parameter_to_explore=parameter_to_explore,
                                pool_eval=pool_eval,
                                name_without_date=True
                                )

        explorator.run()

        if pool_eval:
            p = ProcessExploration(results_path=explorator.results_path)
            resu = p.evaluate(run_all=True, debug=True)

if __name__ == '__main__':
    main()
