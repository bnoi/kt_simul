"""
Launcher can run several simulation with same parameters.
"""

from __future__ import division

import logging
import os
import sys
import datetime
import math
import multiprocessing
import itertools

import pandas as pd

from kt_simul.core.simul_spindle import Metaphase
from kt_simul.io.simuio import SimuIO
from kt_simul.io.simuio import build_tree
from kt_simul.io.xml_handler import ParamTree
from kt_simul.utils.progress import pprogress

logger = logging.getLogger(__name__)


class FolderExistException(Exception):
    pass


class FolderNotExistException(Exception):
    pass


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
                 parallel=True,
                 verbose=True):

        self.verbose = verbose
        if not self.verbose:
            logger.disabled = True
        else:
            logger.disabled = False

        self.simu_path = simu_path

        if not load:
            self.paramtree = paramtree
            self.measuretree = measuretree
            self.initial_plug = initial_plug
            self.parallel = parallel
            self.n_simu = n_simu

            # Create a folder. Raise an exeception if it exists.
            if os.path.isdir(self.simu_path):
                raise FolderExistException("%s exists." % self.simu_path)
            else:
                os.makedirs(self.simu_path)

            self.metaphases_path = []
            self.simus_run = False
            self.digits = int(math.log10(self.n_simu)) + 1

            # Create metadata.h5 files
            # with paramtree, measuretree, intial_plug,
            # parallel, date/time, n_simu

            store = pd.HDFStore(os.path.join(self.simu_path, "metadata.h5"))
            store['params'] = paramtree.to_df()
            store['measures'] = measuretree.to_df()

            metadata = pd.Series({'n_simu': self.n_simu,
                                  'parallel': self.parallel,
                                  'initial_plug': initial_plug,
                                  'datetime': str(datetime.datetime.now())})
            store['metadata'] = metadata
            store.close()

        else:
            if not os.path.isdir(self.simu_path):
                raise FolderNotExistException(
                    "%s does not exists." % self.simu_path)

            store = pd.HDFStore(os.path.join(self.simu_path, "metadata.h5"))
            self.paramtree = ParamTree(root=build_tree(store['params']))
            self.measuretree = ParamTree(root=build_tree(store['measures']),
                                         adimentionalized=False)

            self.n_simu = store['metadata']['n_simu']
            self.parallel = store['metadata']['parallel']
            self.initial_plug = store['metadata']['initial_plug']
            store.close()

            self.metaphases_path = []
            self.digits = int(math.log10(self.n_simu)) + 1

            for i in range(self.n_simu):
                fname = "simu_%s.h5" % (str(i).zfill(self.digits))
                self.metaphases_path.append(
                    os.path.join(self.simu_path, fname))

            self.simus_run = True

    def run(self):
        """
        Run simulations
        """

        if self.simus_run:
            logger.error('Pool has already been simulated.')
            return False

        if self.parallel:
            def init_worker():
                import signal
                signal.signal(signal.SIGINT, signal.SIG_IGN)

            ncore = multiprocessing.cpu_count() + 1
            logger.info('Parallel mode enabled: %i cores will be used to run %i simulations' %
                       (ncore, self.n_simu))
            pool = multiprocessing.Pool(
                processes=ncore, initializer=init_worker)

        i = 0

        # Build arguments list
        simu_parameters = {'paramtree': self.paramtree,
                           'measuretree': self.measuretree,
                           'initial_plug': self.initial_plug,
                           'verbose': False,
                           'reduce_p': True}

        arguments = itertools.izip(itertools.repeat(simu_parameters),
                                   itertools.repeat(self.simu_path),
                                   range(self.n_simu),
                                   itertools.repeat(self.digits))

        try:
            # Launch simulation
            if self.parallel:
                results = pool.imap_unordered(_run_one_simulation, arguments)
            else:
                results = itertools.imap(_run_one_simulation, arguments)

            # Get unordered results and log progress
            for i in range(self.n_simu):
                result = results.next()
                if self.verbose:
                    pprogress((i + 1) / self.n_simu * 100, "(%i / %i)" %
                              (i + 1, self.n_simu))

            if self.verbose:
                pprogress(-1)

        except KeyboardInterrupt:
            pool.terminate()
            pool.join()
            raise CanceledByUserException(
                'Simulation has been canceled by user')

        logger.info("Pool simulations are done")
        self.simus_run = True

    def load_metaphases(self):
        """
        """
        logger.info("Load Metaphase in memory")
        metaphases = []
        for i, fname in enumerate(self.metaphases_path):
            if self.verbose:
                    pprogress((i + 1) / self.n_simu * 100)
            fpath = os.path.join(self.simu_path, fname)
            metaphases.append(SimuIO().read(fpath))
        if self.verbose:
                pprogress(-1)

        return metaphases

def _run_one_simulation(args):
    """
    """
    simu_parameters, simu_path, i, digits = args

    fname = "simu_%s.h5" % (str(i).zfill(digits))
    fpath = os.path.join(simu_path, fname)
    meta = Metaphase(**simu_parameters)
    meta.simul()
    SimuIO(meta).save(fpath)
    return (i, fname)
