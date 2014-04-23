"""
MultiPool module
"""

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import datetime
import logging
import itertools
import gc

import pandas as pd
import numpy as np

from kt_simul.io.simuio import build_tree
from kt_simul.utils.progress import pprogress
from kt_simul.io.xml_handler import ParamTree
from kt_simul.core import parameters
from kt_simul.pool import Pool

PARAMFILE = parameters.PARAMFILE
MEASUREFILE = parameters.MEASUREFILE

log = logging.getLogger(__name__)


class FolderExistException(Exception):
    pass


class FolderNotExistException(Exception):
    pass


class MultiPool:
    """
    MultiPool

    TODO: Add doc
    """

    def __init__(self, multi_pool_path,
                 parameters=None,
                 trees=None,
                 paramtree=None,
                 measuretree=None,
                 load=False,
                 n_simu=10,
                 initial_plug='random',
                 force_parameters=[],
                 parallel=True,
                 verbose=True):

        self.verbose = verbose
        if not self.verbose:
            log.disabled = True
        else:
            log.disabled = False

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
            self.force_parameters = force_parameters
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
                                  'force_parameters': str(force_parameters),
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
            try:
                self.force_parameters = store['metadata']['force_parameters']
            except KeyError:
                self.force_parameters = []
            store.close()

            self.load_pools()

        self.parameter_labels = self.simus_path.index.values

    def run(self,):
        """
        """
        if self.simus_run:
            log.error('Pool has already been simulated.')
            return False

        n = self.simus_path.shape[0]
        log.info('Run simulations for %i different set of parameters' % n)
        log.info('Each set of parameters runs %s simulations' % self.n_simu)

        sys.stdout.flush()
        for i, (parameters, rows) in enumerate(self.simus_path.iterrows()):
            if self.verbose:
                pprogress((i) / n * 100)
                sys.stdout.flush()

            paramtree = self.paramtree
            measuretree = self.measuretree

            if isinstance(parameters, np.ndarray):
                parameters = parameters.tolist()

            names = list(map(lambda x: x[0], self.parameters))
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
                           'initial_plug': self.initial_plug,
                           'force_parameters': self.force_parameters,
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

        self.load_pools()
        log.info("Simulations are done")
        self.simus_run = True

    def load_pools(self):
        """
        """
        self.pools = []
        self.parameters = []
        for parameters, pool_path in self.simus_path.iterrows():
            simu_path = os.path.join(self.multi_pool_path, pool_path['relpath'])
            pool_params = {'load': True,
                           'simu_path': simu_path,
                           'verbose': True
                           }
            self.pools.append(Pool(**pool_params))
        return self.pools

    def _build_path(self):
        names_set = itertools.repeat(list(map(lambda x: x[0], self.parameters)))
        values_set = list(
            itertools.product(*map(lambda x: x[1], self.parameters)))
        n_floats_set = itertools.repeat(list(map(lambda x: x[2], self.parameters)))

        index = []
        pool_path = []

        for names, values, nfloats in zip(names_set, values_set, n_floats_set):
            relpath = ""
            for name, value, nfloat in zip(names, values, nfloats):
                relpath = os.path.join(
                    relpath, "%s_%.*f" % (name, nfloat, value))
            # abspath = os.path.join(self.multi_pool_path, relpath)
            index.append(values)
            pool_path.append(relpath)

        self.simus_path = pd.DataFrame(pool_path, columns=['relpath'])
        self.simus_path.index = pd.MultiIndex.from_tuples(
            index, names=list(map(lambda x: x[0], self.parameters)))

    def set_labels(self, labels):
        """
        """

        self.parameter_labels = labels

    def reorder(self, new_order):
        """
        """

        self.parameter_labels = [self.parameter_labels[i] for i in new_order]
        self.pools = [self.pools[i] for i in new_order]
