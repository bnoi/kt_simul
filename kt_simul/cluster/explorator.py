import os
import datetime
import shutil
import logging
import copy
import json
import gc

from kt_simul.io.xml_handler import ParamTree
from kt_simul.core import parameters
from kt_simul.cluster import Launcher
from kt_simul.core.simul_spindle import PARAMFILE, MEASUREFILE
from kt_simul.cluster.process import Process

logger = logging.getLogger(__name__)


class Explorator:
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
                 parameter_to_explore={},
                 pool_eval=False,
                 name_without_date=False):
        """
        """

        if not parameter_to_explore:
            raise Exception("No parameter to explore found !")

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

        self.parameter_to_explore = parameter_to_explore

        self.results_path = results_path
        self.nsimu = nsimu
        self.nparam = len(self.parameter_to_explore['vector'])
        self.total_simu = self.nsimu * self.nparam

        self.last_progress = 0

        self.name = name
        self.ncore = ncore
        self.pool_eval = pool_eval

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

        # Remove existing directory of it exist
        if os.path.exists(self.results_path):
            shutil.rmtree(self.results_path)
        os.makedirs(self.results_path)

        # Redirect log to run.log
        logfile = os.path.join(self.results_path, "run.log")
        handler = logging.FileHandler(logfile)
        logging.getLogger(__name__).addHandler(handler)

    def run(self):
        """
        """

        logger.info("Starting %i simulations for %s parameters (%s total simulations)"\
            % (self.nsimu, self.nparam, self.total_simu))

        param_name = self.parameter_to_explore['name']
        for param_value in self.parameter_to_explore['vector']:

            name = "%s_%.2f" % (param_name, param_value)
            fullpath = os.path.join(self.results_path, name)

            paramtree_tmp = copy.deepcopy(self.paramtree)
            paramtree_tmp.change_dic(param_name, param_value)

            l = Launcher(self.results_path,
                         self.nsimu,
                         name=name,
                         ncore=self.ncore,
                         paramtree=paramtree_tmp,
                         measuretree=self.measuretree,
                         name_without_date=True)

            l.run()

            del paramtree_tmp
            del l
            gc.collect()

            if self.pool_eval:
                p = Process(results_path=fullpath)
                p.evaluate(debug=True, run_all=True, verbose=False)
                del p
                gc.collect()

        self.create_log()

    def create_log(self):
        """
        Create logfile in results folder
        """

        log = {}
        log["name"] = self.name
        log["nsimu"] = self.nsimu
        log["parameter vector length"] = len(self.parameter_to_explore['vector'])
        log["total simu number"] = self.total_simu
        log["parameter to explore name"] = self.parameter_to_explore['name']
        log["parameter vector"] = self.parameter_to_explore['vector'].tolist()

        log["duration"] = self.paramtree.absolute_dic['span']
        log["dt"] = self.paramtree.absolute_dic['dt']
        log["num_steps"] = int(log["duration"] / log["dt"])

        f = open(os.path.join(self.results_path, "simu.log"), 'w')
        f.write(json.dumps(log, sort_keys=True, indent=4))
        f.close()


