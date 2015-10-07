import logging
import gc
import itertools
import os
import imp

import pandas as pd
import joblib as jb

from ..core.simul_spindle import Metaphase
from ..io.simuio import SimuIO
from ..utils.progress import print_progress
from . import helper
from ..io import ParamTree


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

        self.simu_path = None
        self.ref_path = None

        self.force_parameters = []
        self.initial_plug = 'random'

        self.analyzer = None

    def set_path(self, path):
        """
        """
        self.path = path

        if not os.path.exists(path):
            os.makedirs(path)

        import logging
        imp.reload(logging)
        logging.basicConfig(filename=os.path.join(path, "progress.log"),
                            level=logging.INFO,
                            format='%(asctime)s:%(levelname)s:%(message)s')

        logging.getLogger('kt_simul').setLevel(logging.WARNING)

        self.simu_path = os.path.join(path, 'simus')
        if not os.path.exists(self.simu_path):
            os.makedirs(self.simu_path)

        self.ref_path = os.path.join(self.path, 'reference.h5')

    def load_params(self):
        """
        """

        self.params_matrix = pd.read_hdf(self.ref_path, 'params_matrix')

        self.paramtree = pd.read_hdf(self.ref_path, 'params')
        self.paramtree = ParamTree(df=self.paramtree)

        self.measuretree = pd.read_hdf(self.ref_path, 'measures')
        self.measuretree = ParamTree(df=self.measuretree, adimentionalized=False)

    def load_simus(self):
        """
        """

        self.load_params()

        all_metas = {}

        for label in self.params_matrix['label']:

            simus_path = os.path.join(self.simu_path, '{}.h5'.format(label))
            print(label)

            for i in range(self.n_simus):
                print_progress(int(i * 100 / self.n_simus))

                ioo = SimuIO()
                simu = ioo.read(simus_path, simu_id="simu_{}".format(i + 1),
                                paramtree=self.paramtree,
                                measuretree=self.measuretree)
                del ioo
                yield label, simu

            print_progress(-1)

    def set_params_from_list(self, params):
        """
        """

        if os.path.isfile(self.ref_path):
            os.remove(self.ref_path)

        self.params_matrix = pd.DataFrame(params)
        self.params_matrix.to_hdf(self.ref_path, 'params_matrix')

        self.n_simus_total = len(self.params_matrix) * self.n_simus

    def set_params_from_vector(self, params):
        """
        """

        if os.path.isfile(self.ref_path):
            os.remove(self.ref_path)

        labels = [p['label'] for p in params]
        self.params_matrix = [p['range'] for p in params]
        self.params_matrix = pd.DataFrame(list(itertools.product(*self.params_matrix)),
                                          columns=labels)

        self.n_simus_total = len(self.params_matrix) * self.n_simus

    def run(self, analyzer=helper.simu_analyzer, n_jobs=-1,
            backend="multiprocessing", save=False):
        """
        """

        pd.DataFrame(self.paramtree.params).to_hdf(self.ref_path, 'params')
        pd.DataFrame(self.measuretree.params).to_hdf(self.ref_path, 'measures')

        jobs = []
        i = 1
        for param_id, params in self.params_matrix.iterrows():

            for simu_id in range(self.n_simus):

                path = os.path.join(self.simu_path, "{}_simu_{}.h5".format(params['label'],
                                                                           simu_id + 1))
                if os.path.isfile(path):
                    os.remove(path)

                if save:
                    save = path

                kwargs = {"i": i,
                          "n_total": self.n_simus_total,
                          "params": params.dropna(),
                          "paramtree": self.paramtree.copy(),
                          "measuretree": self.measuretree.copy(),
                          "force_parameters": self.force_parameters,
                          "initial_plug": self.initial_plug,
                          "analyzer": analyzer,
                          "save": save
                          }

                jobs.append(jb.delayed(run_simu)(**kwargs))
                i += 1

        p = jb.Parallel(n_jobs=n_jobs, verbose=11, backend=backend)
        results = p(jobs)

        results = list(itertools.chain(*results))
        self.results = pd.DataFrame(results)
        self.results.set_index('label', inplace=True)
        self.results.to_hdf(os.path.join(self.path, "results.h5"), 'df')

        # Merge simu files
        if save:
            for param_id, params in self.params_matrix.iterrows():

                simu_path = os.path.join(self.simu_path, "{}.h5".format(params['label']))
                store = pd.HDFStore(simu_path)

                for simu_id in range(self.n_simus):

                    path = os.path.join(self.simu_path, "{}_simu_{}.h5".format(params['label'],
                                                                               simu_id + 1))

                    local_store = pd.HDFStore(path)

                    for k in local_store.keys():
                        store['/simu_{}{}'.format(simu_id + 1, k)] = local_store[k]

                    local_store.close()
                    os.remove(path)

                store.close()


        logging.info('Done')


def run_simu(i, n_total,
             params,
             paramtree,
             measuretree,
             force_parameters,
             initial_plug,
             save=False,
             analyzer=analyzers.simu_analyzer):
    """
    """

    logging.info("Processing simulation: {}/{}".format(i, n_total))

    current_paramtree = paramtree.copy()
    current_measuretree = measuretree.copy()

    for k, v in params.iteritems():
        current_paramtree[k] = v
        if k in ['k_d0', 'kappa_c', 'Fmz', 'Vmz']:
            force_parameters.append(k)

    simu = Metaphase(verbose=False,
                     paramtree=current_paramtree,
                     measuretree=current_measuretree,
                     force_parameters=force_parameters,
                     initial_plug=initial_plug)

    simu.simul()

    if save:
        SimuIO(simu).save(save, save_tree=False)

    results = analyzer(simu, params)

    del simu
    gc.collect()

    return results
