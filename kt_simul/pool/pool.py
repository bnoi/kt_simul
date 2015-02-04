"""
Launcher can run several simulation with same parameters.
"""

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import logging
import os
import datetime
import math
import multiprocessing
import itertools
import tempfile
import uuid
import gc

import pandas as pd

from ..core.simul_spindle import Metaphase
from ..io.simuio import SimuIO
from ..io import ParamTree
from ..utils.progress import pprogress
from ..utils.dict import sanitize_dict
from ..utils.dict import guess_number_dict

log = logging.getLogger(__name__)


class Pool:

    """
    Pool launchs a pool of simulation with same parameters.
    """

    def __init__(self, simu_path,
                 paramtree=None,
                 measuretree=None,
                 load=False,
                 n_simu=100,
                 initial_plug='random',
                 force_parameters=[],
                 parallel=True,
                 verbose=True):

        self.verbose = verbose
        if not self.verbose:
            log.disabled = True
        else:
            log.disabled = False

        self.simu_path = simu_path
        self.simu_full_path = self.simu_path

        if not load:
            self.paramtree = paramtree
            self.measuretree = measuretree
            self.initial_plug = initial_plug
            self.force_parameters = force_parameters
            self.parallel = parallel
            self.n_simu = n_simu

            if os.path.isfile(self.simu_path):
                raise Exception("%s exists." % self.simu_path)

            self.metaphases_path = []
            self.simus_run = False
            self.digits = int(math.log10(self.n_simu)) + 1

            # Create metadata.h5 files
            # with paramtree, measuretree, intial_plug,
            # parallel, date/time, n_simu

            store = pd.HDFStore(self.simu_full_path)
            store['params'] = paramtree.params
            store['measures'] = measuretree.params

            metadata = {'n_simu': self.n_simu,
                        'parallel': self.parallel,
                        'initial_plug': initial_plug,
                        'force_parameters': str(force_parameters),
                        'datetime': str(datetime.datetime.now())}
            metadata = sanitize_dict(metadata)
            store['metadata'] = pd.Series(metadata)
            store.close()

        else:
            store = pd.HDFStore(self.simu_full_path)
            self.paramtree = ParamTree(df=store['params'])
            self.measuretree = ParamTree(df=store['measures'], adimentionalized=False)

            metadata = guess_number_dict(store['metadata'].to_dict())
            self.n_simu = metadata['n_simu']
            self.parallel = metadata['parallel']
            self.initial_plug = metadata['initial_plug']
            try:
                self.force_parameters = metadata['force_parameters']
            except KeyError:
                self.force_parameters = []
            store.close()

            self.metaphases_path = []
            self.digits = int(math.log10(self.n_simu)) + 1

            self.simus_run = True



    def run(self):
        """
        Run simulations
        """

        if self.simus_run:
            log.error('Pool has already been simulated.')
            return False

        if self.parallel:
            def init_worker():
                import signal
                signal.signal(signal.SIGINT, signal.SIG_IGN)

            ncore = multiprocessing.cpu_count() + 1
            log.info('Parallel mode enabled: %i cores will be used to run %i simulations' %
                       (ncore, self.n_simu))
            pool = multiprocessing.Pool(
                processes=ncore, initializer=init_worker)

        i = 0

        # Build arguments list
        simu_parameters = {'paramtree': self.paramtree,
                           'measuretree': self.measuretree,
                           'initial_plug': self.initial_plug,
                           'force_parameters': self.force_parameters,
                           'verbose': False}

        if self.parallel:
            simus_path = []
            unique_id = uuid.uuid4()
            for i in range(self.n_simu):
                name = os.path.join(tempfile.gettempdir(), "simu_{}_{}.h5")
                simus_path.append(name.format(unique_id, i))
        else:
            simus_path = itertools.repeat(self.simu_full_path)

        arguments = zip(itertools.repeat(simu_parameters),
                        simus_path,
                        range(self.n_simu),
                        itertools.repeat(self.digits))

        try:
            # Launch simulation
            if self.parallel:
                results = pool.imap_unordered(_run_one_simulation, arguments)
            else:
                results = map(_run_one_simulation, arguments)

            # Get unordered results and log progress
            for i in range(self.n_simu):
                next(results)
                if self.verbose:
                    pprogress((i + 1) / self.n_simu * 100, "(%i / %i)" %
                              (i + 1, self.n_simu))

            if self.verbose:
                pprogress(-1)

            pool.terminate()
            pool.join()

        except KeyboardInterrupt:
            pool.terminate()
            pool.join()
            raise CanceledByUserException(
                'Simulation has been canceled by user')

        if self.parallel:
            store = pd.HDFStore(self.simu_full_path)
            for simu_tmp in simus_path:
                st = pd.HDFStore(simu_tmp)
                for key in st.keys():
                    store[key] = st[key]
                st.close()
                os.remove(simu_tmp)
            store.close()

            del pool

        log.info("Pool simulations are done")
        self.simus_run = True

    def load_metaphases(self):
        """
        """
        all_simu_id = [get_simu_id(i, self.digits) for i in range(self.n_simu)]

        for simu_id in all_simu_id:

            ioo = SimuIO()
            meta = ioo.read(self.simu_full_path,
                            simu_id=simu_id,
                            paramtree=self.paramtree,
                            measuretree=self.measuretree)
            del ioo
            yield meta

    def load_metaphase_parallel(self, pre_processing_func=None , verbose=True):
        """
        """
        all_fpath = [os.path.join(self.simu_path, fname) for i, fname in enumerate(self.metaphases_path)]
        args = zip(all_fpath, itertools.repeat(self.paramtree), itertools.repeat(self.measuretree))

        proc_pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        results = []
        for i, meta in enumerate(proc_pool.imap_unordered(_load_metaphase_single, args)):
            if verbose:
                pprogress(int(i * 100 / self.n_simu))
            print(meta.KD.chromosomes)
            if pre_processing_func:
                results.append(pre_processing_func(meta))
            else:
                results.append(meta)

        if verbose:
            pprogress(-1)

        return results


def _run_one_simulation(args):
    """
    """
    simu_parameters, simu_path, i, digits = args
    simu_id = get_simu_id(i, digits)
    meta = Metaphase(**simu_parameters)
    meta.simul()
    SimuIO(meta).save(simu_path, simu_id, save_tree=True)
    del meta.KD
    del meta
    gc.collect()
    return (i, simu_id)

def _load_metaphase_single(args):
    fpath, paramtree, measuretree = args
    meta = SimuIO().read(fpath,
                         paramtree=paramtree,
                         measuretree=measuretree)
    return meta


def get_simu_id(i, digits):
    return "simu_%s" % (str(i).zfill(digits))
