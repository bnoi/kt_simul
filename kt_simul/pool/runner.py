import logging
import gc
import itertools
import os

import pandas as pd
import joblib as jb

from ..core.simul_spindle import Metaphase
from . import analyzers


class Runner:
    """
    """

    def __init__(self):
        """
        """
        self.paramtree = None
        self.measuretree = None
        self.n_simus = None
        self.n_simus_total = None
        self.params_matrix = None

        self.path = None
        self.results = None

        self.force_parameters = []
        self.initial_plug = 'random'

        self.analyzer = None

    def set_path(self, path):
        """
        """
        self.path = path

        if not os.path.exists(path):
            os.makedirs(path)

        logging.basicConfig(filename=os.path.join(path, "progress.log"),
                            level=logging.INFO,
                            format='%(asctime)s:%(levelname)s:%(message)s')
        logging.getLogger('kt_simul').setLevel(logging.WARNING)

        simu_path = os.path.join(path, 'simus')
        if not os.path.exists(simu_path):
            os.makedirs(simu_path)

    def set_params_from_list(self, params):
        """
        """

        self.params_matrix = pd.DataFrame(params)
        self.params_matrix.to_hdf(os.path.join(self.path, 'reference.h5'), 'params_matrix')

        self.n_simus_total = len(self.params_matrix) * self.n_simus

    def set_params_from_vector(self, params):
        """
        """

        labels = [p['name'] for p in params]
        self.params_matrix = [p['range'] for p in params]
        self.params_matrix = pd.DataFrame(list(itertools.product(*self.params_matrix)),
                                          columns=labels)

        self.n_simus_total = len(self.params_matrix) * self.n_simus

    def run(self, analyzer=analyzers.simu_analyzer, n_jobs=-1, backend="multiprocessing"):
        """
        """

        jobs = []
        i = 1
        for param_id, params in self.params_matrix.iterrows():
            for _ in range(self.n_simus):
                kwargs = {"i": i,
                          "n_total": self.n_simus_total,
                          "params": params,
                          "paramtree": self.paramtree,
                          "measuretree": self.measuretree,
                          "force_parameters": self.force_parameters,
                          "initial_plug": self.initial_plug}

                jobs.append(jb.delayed(run_simu)(analyzer=analyzer, **kwargs))
                i += 1

        p = jb.Parallel(n_jobs=n_jobs, verbose=11, backend=backend)
        results = p(jobs)

        results = list(itertools.chain(*results))
        self.results = pd.DataFrame(results)
        self.results.to_hdf(os.path.join(self.path, "results.h5"), 'df')

        logging.info('Done')


def run_simu(i, n_total,
             params,
             paramtree,
             measuretree,
             force_parameters,
             initial_plug,
             analyzer=analyzers.simu_analyzer):
    """
    """

    logging.info("Processing simulation: {}/{}".format(i, n_total))

    current_paramtree = paramtree.copy()
    current_measuretree = measuretree.copy()

    for k, v in params.iteritems():
        current_paramtree[k] = v

    simu = Metaphase(verbose=False,
                     paramtree=current_paramtree,
                     measuretree=current_measuretree,
                     force_parameters=force_parameters,
                     initial_plug=initial_plug)

    simu.simul()

    results = analyzer(simu, params)

    del simu
    gc.collect()

    return results
