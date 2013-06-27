# -*- coding: utf-8 -*-
"""
Pool evaluations are small plugin that can be automatically
launch after severals simulations.

PoolEvaluation take a path where .zip files are stored. To open simulations
uses kt_simul.simuio.SimuIO().read()
"""

import sys
import os
import logging
import re
import json

from kt_simul.utils.filesystem import tree
from kt_simul.io.simuio import SimuIO
from kt_simul.core.simul_spindle import Metaphase

EVAL_PATH = os.path.abspath(os.path.dirname(__file__))

logger = logging.getLogger(__name__)


def find_pool_evaluations(name = None, groups=[], run_all = False):
    """
    This function return a list of PoolEvaluation classes that are enabled. In the future, it will have the capability to filter evaluations.

    .. note:: Source from http://www.luckydonkey.com/2008/01/02/python-style-plugins-made-easy/

    .. todo:: add selector and filter

    :param groups: Select plugins from these groups
    :type groups: str

    :rtype: list
    :return: a list if classes that are subclasses of cls
    """

    path = EVAL_PATH
    cls = PoolEvaluation

    subclasses = []

    def look_for_subclass(modulename):
        module = __import__(modulename)

        # walk the dictionaries to get to the last one
        d = module.__dict__
        for m in modulename.split('.')[1:]:
            try:
                d = d[m].__dict__
            except:
                pass

        #look through this dictionary for things
        #that are subclass of Job
        #but are not Job itself
        for key, entry in d.items():
            if key == cls.__name__:
                continue
            try:
                if issubclass(entry, cls):
                    subclasses.append(entry)
            except TypeError:
                #this happens when a non-type is passed in to issubclass. We
                #don't care as it can't be a subclass of Job if it isn't a
                #type
                continue

    sys.path.insert(0, path)

    plugin_files = []
    for x in tree(path):
        if x.endswith(".py"):
            mod = x[:-3].replace(os.sep, ".")
            plugin_files.append(mod)

    sys.path.insert(0, path)

    for plugin in plugin_files:
        look_for_subclass("kt_simul.analysis.pool_evaluations." + plugin)

    # Filter plugin here
    evaluations = []
    for cls in subclasses:
        if cls.enable:
            if not run_all:
                if not groups and cls.group == None and not name:
                    evaluations.append(cls)
                elif name and cls.name == name:
                    evaluations.append(cls)
                else:
                    if cls.group in groups:
                        evaluations.append(cls)
            else:
                evaluations.append(cls)

    return evaluations


class RunFunctionNotImplemented(Exception):
    pass


class PoolEvaluation(object):
    """
    Abstract class. All pool evaluation plugins need to write a
    class that inherit from PoolEvaluation.

    .. warning:: run() method need to be implemented.
    """

    def __init__(self,):
        pass

    def run(self, *args):
        raise RunFunctionNotImplemented("run() method need to be implemented \
            in your PoolEvaluation plugin.")

    def get_nsimu(self, path):
        """
        Return the number of simulations inside a path. This method look to
        the file simu.log
        """

        logpath = os.path.join(path, "simu.log")
        log = json.load(open(logpath))

        return log["number_of_simulations"]

    def get_name(self, path):
        """
        Return the name of simulations inside a path. This method look to
        the file simu.log
        """

        logpath = os.path.join(path, "simu.log")
        log = json.load(open(logpath))

        return log["name"]

    def get_simus_params(self, path):
        """
        Return empty Metaphase object to give simulation parameters
        before iterating on all simulations
        """

        params = os.path.join(path, "params.xml")
        measures = os.path.join(path, "measures.xml")
        meta = Metaphase(paramfile=params, measurefile=measures, verbose=False)
        return meta

    def get_params(self, path):
        """
        """

        params = os.path.join(path, "params.xml")
        measures = os.path.join(path, "measures.xml")
        meta = Metaphase(paramfile=params, measurefile=measures, verbose=False)
        return meta.KD.params



    def iter_simulations(self, path, nsimu=-1, print_progress=True):
        """
        Iterator on zip file which contain simulationr results.

        :return simu_id, meta: Simulation id and an instance of
                               Metaphase object
        """

        i = 0
        last = -1
        for simu in os.listdir(path):
           # Find simu id with zip filename
            simu_id = ""
            try:
                simu_id = int(re.search('simu\_(.*)\.zip', simu).group(1))
            except:
                logger.info("Filename %s is not well formated" % simu)
                break

            # Create Metpahase object
            meta = SimuIO().read(os.path.join(path, simu))

            if meta:

                if print_progress:
                    progress = "%0.0f" % ((i * 100.) / nsimu)
                    if last != progress:
                        logger.info("Progression: %s%%" % progress)
                        last = progress
                    i += 1

                yield simu_id, meta

            del meta

__all__ = []
for ev in find_pool_evaluations():
    __all__.append(ev)

def get():
    """
    Return all available pool_evaluations
    """
    evals = [e.name for e in find_pool_evaluations(run_all = True)]
    return evals

__all__ += ['find_pool_evaluations', 'PoolEvaluation', 'RunFunctionNotImplemented', 'get']
