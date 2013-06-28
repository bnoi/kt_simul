import logging
import os

import numpy as np

from kt_simul.analysis.pool_evaluations import find_pool_evaluations

logger = logging.getLogger(__name__)


class Process:
    """
    """

    def __init__(self, results_path):
        """
        """
        self.results_path = results_path

        # Where simulations raw data are stored
        self.raw_path = os.path.join(self.results_path, 'raw')

        # Where to store analysis files (plot, etc)
        self.eval_results = os.path.join(self.results_path, 'analysis')

        self.observations = {}

    def evaluate(self, name=None, groups=[], draw=True, debug=False, run_all=False, verbose=True):
        """
        """

        logger.info("Starting pool evaluations")
        all_pool_evaluations = find_pool_evaluations(name=name, groups=groups, run_all=run_all)

        if not all_pool_evaluations:
            logger.info("No pool evaluations found")
            return False

        for pool_evaluation in all_pool_evaluations:
            logger.info("Running %s" % pool_evaluation.name)
            if debug:
                result = pool_evaluation().run(self.results_path,
                                                self.raw_path,
                                                self.eval_results,
                                                draw=draw,
                                                verbose=verbose)
                logger.info("%s done" % pool_evaluation.name)
            else:
                try:
                    result = pool_evaluation().run(self.results_path,
                                                self.raw_path,
                                                self.eval_results,
                                                draw=draw,
                                                verbose=verbose)
                    logger.info("%s done" % pool_evaluation.name)
                except Exception as e:
                    result = np.nan
                    logger.info("%s returns errors : %s" % (pool_evaluation.name, e))

            if name and not run_all:
                return result
            else:
                current_name = pool_evaluation.name
                self.observations[current_name] = result

        logger.info("All pool evaluations processed")

        del result
        return self.observations
