"""
Launcher can run several simulation with same parameters.
It handle data storing.
"""

import os
import time
import datetime
import json
import shutil
import gc
import logging
import multiprocessing
from multiprocessing import Pool, Queue

from kt_simul.io.xml_handler import ParamTree
from kt_simul.core.simul_spindle import Metaphase, PARAMFILE, MEASUREFILE
from kt_simul.core import parameters
from kt_simul.utils.size import get_folder_size
from kt_simul.io.simuio import SimuIO
from ..utils.ipython import in_ipython

logger = logging.getLogger(__name__)


class Launcher:
    """
    """

    def __init__(self, results_path,
                 nsimu,
                 name="",
                 ncore=None,
                 paramtree=None,
                 measuretree=None,
                 paramfile=PARAMFILE,
                 measurefile=MEASUREFILE,
                 verbose=True,
                 name_without_date=False):
        """

        :results_path: The path where simulation results are stored
        :type results_path: string

        :nsimu: Number of simulation to run
        :type nsimu: int

        :ncore: Number of process to launch at the same time
        :type ncore: int

        :param paramtree: The paramtree contains the parameters for the simulation
            if paramtree is None, the parameters are read
            from the file paramfile. Defaults to None.
        :type paramtree: ParamTree instance or None

        :param measuretree: The measuretree contains the observed characteristics
            of the mitosis e.g. metaphase spindle elongation rate, etc.
            if measuretree is None, the measures are read from the file
            indicated by the measurefile argument. Defaults to None.
        :type measuretree: ParamTree instance or None

        :param paramfile: Path to a xml file to read the parameters from. Defaults to the
            file params.xml in the module's default directory. Other parameter
            files can be produced by editing and changing the default one.
            If the paramtree argument is not None,  paramfile is ignored
        :type paramfile: string

        :param measurefile: Path to a xml file to read the measures from. Defaults to the
            file measures.xml in the module's default directory.
            Other measure files can be produced by editing and changing
            the default one. If the measuretree argument is not None, measurefile
            is ignored
        :type measurefile: string


        """

        # Enable or disable log console
        self.verbose = verbose
        logger = logging.getLogger(__name__)
        if not self.verbose:
            logger.disabled = True
        else:
            logger.disabled = False

        if paramtree is None:
            self.paramtree = ParamTree(paramfile)
        else:
            self.paramtree = paramtree
        if measuretree is None:
            self.measuretree = ParamTree(measurefile, adimentionalized=False)
        else:
            self.measuretree = measuretree

        # Reduce parameters
        parameters.reduce_params(self.paramtree, self.measuretree)

        self.results_path = results_path
        self.nsimu = nsimu
        if not ncore:
            self.ncore = multiprocessing.cpu_count()
        else:
            self.ncore = ncore

        self.last_progress = 0

        self.name = name

        # Make result directory
        if not os.path.exists(self.results_path):
            os.makedirs(self.results_path)

        # Results directory according to date and time
        if name_without_date and name:
            dirname = name
        else:
            now = datetime.datetime.now()
            if name:
                dirname = now.strftime("%Y.%m.%d") + "_" + name
            else:
                dirname = now.strftime("%Y.%m.%d")

        self.results_path = os.path.join(self.results_path, dirname)

        # TODO: Ugly fix
        if os.path.isdir(self.results_path):
            self.results_path += "_2"

        # Remove existing directory of it exist
        if os.path.exists(self.results_path):
            shutil.rmtree(self.results_path)
        os.makedirs(self.results_path)

        # Save params.xml and measures.xml in results dir
        self.paramtree.save(os.path.join(self.results_path, "params.xml"))
        self.measuretree.save(os.path.join(self.results_path, "measures.xml"))

        self.raw_results_path = os.path.join(self.results_path, "raw")
        os.makedirs(self.raw_results_path)

        # Redirect log to run.log
        logfile = os.path.join(self.results_path, "run.log")
        handler = logging.FileHandler(logfile)
        logging.getLogger(__name__).addHandler(handler)

    def run(self):
        """
        """

        logger.info("Starting %i simulations on %i core" % \
            (self.nsimu, self.ncore))

        if in_ipython():
            parallel = False
        else:
            parallel = True

        if parallel and os.name == 'posix':
            os.system("taskset -p 0xff %d" % os.getpid())

        params = [self.paramtree, self.measuretree, self.verbose]
        all_params = []
        for i in range(self.nsimu):
            all_params.append([i, self.raw_results_path] + params)

        queue = Queue()
        self.start_time = time.time()
        p = Pool(self.ncore, run_init, [queue], maxtasksperchild=4)
        p.map_async(run_one, all_params, chunksize=1)

        # Monitoring simulations
        done = False
        simu_ids = range(self.nsimu)
        while not done:
            mess = queue.get()
            if mess["state"] == "stop":
                simu_ids.remove(mess["id"])
                self.log_progress(len(simu_ids))

            if not simu_ids:
                done = True
        self.total_time = time.time() - self.start_time

        p.close()

        results_size = get_folder_size(self.results_path)
        logger.info("Simulations are done")
        logger.info("Results are stored in %s (%s MB)" % (self.results_path, results_size))

        self.create_log()

        del p
        del queue
        gc.collect()

    def log_progress(self, simu_left, precision=2):
        """
        Log the progression of the simulations

        :param simu_left: Number of simulation remaining
        :type simu_left: int

        """
        progress = 100.0 - ((simu_left * 100.0) / self.nsimu)
        progress = round(progress, precision)
        # Don't print the same progression !
        if progress > self.last_progress:
            time_remaining, spent_time = self.estimate_remaining_time(simu_left)
            output = "Progression : %0.2f%% | " % progress
            output += "ETA : %s | " % time_remaining
            output += "Spent time : %s" % spent_time
            logger.info(output)
            self.last_progress = progress

    def estimate_remaining_time(self, simu_left):
        """
        Estimate remaining time
        """
        delta_time = time.time() - self.start_time
        delta_simu = self.nsimu - simu_left
        estimate_time = (simu_left * delta_time) / delta_simu
        estimate_time = time.strftime('%H:%M:%S', time.gmtime(estimate_time))
        spent_time = time.strftime('%H:%M:%S', time.gmtime(delta_time))
        return estimate_time, spent_time

    def create_log(self):
        """
        Create logfile in results folder
        """

        log = {}
        log["name"] = self.name
        log["duration"] = self.paramtree.absolute_dic['span']
        log["dt"] = self.paramtree.absolute_dic['dt']
        log["num_steps"] = int(log["duration"] / log["dt"])
        log["number_of_simulations"] = self.nsimu
        log["number_of_simulations"] = self.nsimu
        log["spent_time"] = time.strftime('%H:%M:%S', time.gmtime(self.total_time))
        log["results_folder_size_in_MB"] = get_folder_size(self.results_path)

        f = open(os.path.join(self.results_path, "simu.log"), 'w')
        f.write(json.dumps(log, sort_keys=True, indent=4))
        f.close()


def run_one(args):
    return _run_one(*args)


def _run_one(simu_id, result_path, paramtree, measuretree, verbose):
    """
    This function need to be outisde Launcher class to allow
    multiprocessing module to work
    """
    queue = _run_one.q

    queue.put({"id": simu_id, "state": "start"})
    meta = Metaphase(verbose=False,
                     paramtree=paramtree,
                     measuretree=measuretree,
                     reduce_p=True)
    meta.simul()
    meta.evaluate()

    # Build filename
    simu_path = os.path.join(result_path, "simu_%06d.zip" % simu_id)

    # Write simulation result
    io = SimuIO(meta)
    io.save(simufname=simu_path, verbose=False)

    queue.put({"id" : simu_id, "state" : "stop" })

    del meta
    del io
    gc.collect()
    return True


def run_init(q):
    """
    Trick function to pass Queue in _run_one
    """
    _run_one.q = q
