# -*- coding: utf-8 -*-
import gc

import numpy as np

from kt_simul.cluster import Explorator
from kt_simul.cluster.process_exploration import ProcessExploration
from kt_simul.io.xml_handler import ParamTree


def one_explo(param, param_to_explore):

    PARAMFILE = "params.xml"
    MEASUREFILE = "measures.xml"

    result_path = "/media/thor/data/ktsimu/small_explo_ld_slope_ld0"
    result_path = "/home/hadim/local/data/ktsimu/explo_ld_slope_ld0"
    number_simu = 500
    name = "explo_%s_%s" % (param['name'], str(param['value']))
    ncore = 4

    paramtree = ParamTree(PARAMFILE)
    paramtree.change_dic(param['name'], param['value'])

    explorator = Explorator(result_path,
                            number_simu,
                            name=name,
                            ncore=ncore,
                            paramtree=paramtree,
                            measurefile=MEASUREFILE,
                            parameter_to_explore=parameter_to_explore,
                            pool_eval=True,
                            name_without_date=True
                            )
    explorator.run()

    p = ProcessExploration(results_path=explorator.results_path)
    p.evaluate(run_all=True, debug=True)

    del explorator
    del p
    gc.collect()

if __name__ == '__main__':

    pvec = np.arange(0, 2, 0.1)

    parameter_to_explore = {'name': 'ld_slope', 'vector': pvec}
    param = {}
    param['name'] = 'ld0'

    args = []
    for ld0 in pvec:
        param['value'] = ld0
        one_explo(param, parameter_to_explore)
