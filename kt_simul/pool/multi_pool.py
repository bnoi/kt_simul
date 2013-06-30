"""
MultiPool module
"""

from __future__ import division

import os
import sys
import datetime
import logging
import itertools
import gc

import pandas as pd

from kt_simul.io.simuio import build_tree
from kt_simul.utils.progress import pprogress
from kt_simul.io.xml_handler import ParamTree
from kt_simul.core import parameters
from kt_simul.pool import Pool

PARAMFILE = parameters.PARAMFILE
MEASUREFILE = parameters.MEASUREFILE

logger = logging.getLogger(__name__)


class FolderExistException(Exception):
    pass


class FolderNotExistException(Exception):
    pass


class MultiPool:

    def __init__(self, multi_pool_path,
                 parameters,
                 trees,
                 paramtree=None,
                 measuretree=None,
                 load=False,
                 n_simu=10,
                 initial_plug='random',
                 parallel=True,
                 verbose=True):

        self.verbose = verbose
        if not self.verbose:
            logger.disabled = True
        else:
            logger.disabled = False

        self.multi_pool_path = multi_pool_path

        if not load:
            # Create a folder. Raise an exeception if it exists.
            if os.path.isdir(self.multi_pool_path):
                raise FolderExistException("%s exists." % self.multi_pool_path)
            else:
                os.makedirs(self.multi_pool_path)

            self.parameters = parameters
            self.paramtree = paramtree
            self.measuretree = measuretree
            self.n_simu = n_simu
            self.parallel = parallel
            self.initial_plug = initial_plug
            self.trees = trees

            self.simus_run = False

            self._build_path()

            store = pd.HDFStore(
                os.path.join(self.multi_pool_path, "metadata.h5"))
            store['simus_path'] = self.simus_path
            store['params'] = paramtree.to_df()
            store['measures'] = measuretree.to_df()

            metadata = pd.Series({'n_simu': self.n_simu,
                                  'parallel': self.parallel,
                                  'initial_plug': initial_plug,
                                  'datetime': str(datetime.datetime.now())})
            store['metadata'] = metadata
            store.close()
        else:
            if not os.path.isdir(self.multi_pool_path):
                raise FolderNotExistException(
                    "%s does not exists." % self.multi_pool_path)

            self.simus_run = True

            store = pd.HDFStore(
                os.path.join(self.multi_pool_path, "metadata.h5"))
            self.simus_path = store['simus_path']
            self.paramtree = ParamTree(root=build_tree(store['params']))
            self.measuretree = ParamTree(root=build_tree(store['measures']),
                                         adimentionalized=False)

            self.n_simu = store['metadata']['n_simu']
            self.parallel = store['metadata']['parallel']
            self.initial_plug = store['metadata']['initial_plug']
            store.close()

    def run(self,):
        """
        """
        if self.simus_run:
            logger.error('Pool has already been simulated.')
            return False

        n = self.simus_path.shape[0]
        logger.info('Run simulations for %i different set of parameters' % n)
        logger.info('Each set of parameters runs %s simulations' % self.n_simu)
        sys.stdout.flush()
        for i, (parameters, rows) in enumerate(self.simus_path.iterrows()):
            if self.verbose:
                pprogress((i + 1) / n * 100)
                sys.stdout.flush()

            paramtree = self.paramtree
            measuretree = self.measuretree

            names = map(lambda x: x[0], self.parameters)
            for i, parameter in enumerate(parameters):
                if self.trees[i] == 'paramtree':
                    paramtree.change_dic(names[i], parameter)
                elif self.trees[i] == 'measuretree':
                    measuretree.change_dic(names[i], parameter)

            simu_path = os.path.join(self.multi_pool_path, rows['relpath'])

            pool_params = {'load': False,
                           'paramtree': paramtree,
                           'measuretree': measuretree,
                           'n_simu': self.n_simu,
                           'simu_path': simu_path,
                           'parallel': self.parallel,
                           'verbose': False
                           }

            pool = Pool(**pool_params)
            pool.run()

            del pool
            gc.collect()

        if self.verbose:
            pprogress(-1)

        logger.info("Simulations are done")
        self.simus_run = True

    def iter_pools(self, func, verbose=True):
        """
        func will recieves a list of metaphases and list of parameters for this
        pool. This function is usefull because it provides some cleanup memory
        and progressbar.
        """

        results = []
        n = self.simus_path.shape[0]
        for i, (parameters, rows) in enumerate(self.simus_path.iterrows()):
            if verbose:
                pprogress((i + 1) / n * 100)

            relpath = rows['relpath']
            pool_path = os.path.join(self.multi_pool_path, relpath)
            pool_params = {'load': True,
                           'simu_path': pool_path,
                           'verbose': False
                           }
            pool = Pool(**pool_params)
            metas = pool.load_metaphases()

            results.append(func(parameters, metas))

        if verbose:
            pprogress(-1)

        return results

    def _build_path( self):
        names_set = itertools.repeat(map(lambda x: x[0], self.parameters))
        values_set = list(
            itertools.product(*map(lambda x: x[1], self.parameters)))
        n_floats_set = itertools.repeat(map(lambda x: x[2], self.parameters))

        index = []
        pool_path = []

        for names, values, nfloats in zip(names_set, values_set, n_floats_set):
            relpath = ""
            for name, value, nfloat in zip(names, values, nfloats):
                relpath = os.path.join(
                    relpath, "%s_%.*f" % (name, nfloat, value))
            abspath = os.path.join(self.multi_pool_path, relpath)
            index.append(values)
            pool_path.append(relpath)

        self.simus_path = pd.DataFrame(pool_path, columns=['relpath'])
        self.simus_path.index = pd.MultiIndex.from_tuples(
            index, names=map(lambda x: x[0], self.parameters))
