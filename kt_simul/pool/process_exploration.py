import logging
import os

import numpy as np

from kt_simul.analysis.explo_pool_evaluations import find_explo_pool_evaluations

logger = logging.getLogger(__name__)


class ProcessExploration:
    """
    """

    def __init__(self, results_path):
        """
        """
        self.results_path = results_path

        # Retrieve pool simulations folder
        self.pool_folder = self.get_pool_folder(self.results_path)

        # Where to store analysis files (plot, etc)
        self.eval_results = os.path.join(self.results_path, 'analysis')

        # Create eval_results_path if it does not exist
        if not os.path.isdir(self.eval_results):
            os.makedirs(self.eval_results)

    def get_pool_folder(self, path):
        """
        """

        pool = []
        for d in os.listdir(path):
            simud = os.path.join(path, d)

            if os.path.isdir(simud):

                files = set(os.listdir(simud))
                to_check = set(['simu.log', 'raw', 'measures.xml', 'params.xml'])
                if to_check.issubset(files):
                    pool.append(simud)
        return pool

    def evaluate(self, groups=[], debug=False, run_all=False):
        """
        """

        logger.info("Starting exploration pool evaluations")
        all_explo_explo_pool_evaluations = find_explo_pool_evaluations(groups=groups, run_all=run_all)

        if not all_explo_explo_pool_evaluations:
            logger.info("No pool evaluations found")
            return False

        for explo_pool_evaluation in all_explo_explo_pool_evaluations:
            logger.info("Running %s" % explo_pool_evaluation.name)
            if debug:
                result = explo_pool_evaluation().run(self.results_path,
                                                self.pool_folder,
                                                self.eval_results)
                logger.info("%s done" % explo_pool_evaluation.name)
            else:
                try:
                    result = explo_pool_evaluation().run(self.results_path,
                                                self.pool_folder,
                                                self.eval_results)
                    logger.info("%s done" % explo_pool_evaluation.name)
                except Exception as e:
                    result = np.nan
                    logger.info("%s returns errors : %s" % (explo_pool_evaluation.name, e))

        logger.info("All exploration pool evaluations processed")

        del result
        return True
