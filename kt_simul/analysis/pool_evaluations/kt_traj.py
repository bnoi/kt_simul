import numpy as np
import os
import logging

import matplotlib.pyplot as plt

from kt_simul.analysis.pool_evaluations import PoolEvaluation
from kt_simul.draw.plot import dic_plot

logger = logging.getLogger(__name__)


class KtTrajectories(PoolEvaluation):
    """
    """

    name = "Kinetochores Trajectories"
    description = "Make profile of the simulation drawing kinetochores trajectories"
    group = "trajectories"
    enable = True

    def __init__(self,):
        pass

    def run(self, simu_path, raw_path, eval_results_path, draw=True, verbose=True):
        """
        """

        # Create eval_results_path if it does not exist
        if not os.path.isdir(eval_results_path):
            os.makedirs(eval_results_path)

        params = self.get_params(simu_path)
        meta_infos = self.get_simus_params(simu_path)
        num_steps = meta_infos.num_steps
        nsimu = int(self.get_nsimu(simu_path))
        name = self.get_name(simu_path)
        ana_onset = int(params["t_A"])
        N = int(params["N"])

        trajectories = {'spbR': np.zeros((nsimu, num_steps), dtype=float),
                        'spbL': np.zeros((nsimu, num_steps), dtype=float),
                        'chromosomes': [],
                       }

        for i in range(N):
            cen_traj = [np.zeros((nsimu, num_steps), dtype=float),
                        np.zeros((nsimu, num_steps), dtype=float)]
            trajectories['chromosomes'].append(cen_traj)

        if verbose:
            logger.info("Loading data from simulations files")
        for i, (simu_id, meta) in enumerate(self.iter_simulations(raw_path,
                                                        nsimu=nsimu,
                                                        print_progress=verbose)):
            results = meta.evaluate(name="Kinetochores Trajectories", verbose=False)

            trajectories['spbR'][i] = results['spbR']
            trajectories['spbL'][i] = results['spbL']

            for j, ch in enumerate(results['chromosomes']):
                trajectories["chromosomes"][j][0][i] = ch[0]
                trajectories["chromosomes"][j][1][i] = ch[1]

        logger.disabled = False

        if not draw:
            return trajectories

        # Configure and plot the graph
        if verbose:
            logger.info("Plotting results")
        timelapse = meta_infos.timelapse

        plot_data = {}
        plot_data['title'] = "Chromosome kinetochores trajectories"
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time'}
        plot_data['yaxis'] = {'label': 'Relative position to SPB axes center', 'axis': []}
        plot_data['error'] = False
        plot_data['legend'] = False

        # Draw parameters box
        plot_data["params_box"] = [{'name': "Name", 'data': name},
                                   {'name': "Simulations number", 'data': nsimu},
                             ]

        # Add annotation about anaphase onset
        plot_data["annotations"] = []
        plot_data["annotations"].append({'s': 'Anaphase onset: %i' % ana_onset,
                                         'xy': (ana_onset, -2.5),
                                         'xytext': (0, 20),
                                         'textcoords': 'offset points',
                                         'ha': 'center',
                                         'va': 'bottom',
                                         'bbox': dict(boxstyle='round,pad=0.2',
                                                      fc='yellow',
                                                      alpha=0.3),
                                         'arrowprops': dict(arrowstyle='->',
                                                            # connectionstyle='arc3,rad=0.5',
                                                            color='black')
                                         })

        plot_data['yaxis']['axis'].append({'data': trajectories['spbR'].mean(axis=0),
                                            'color': 'black',
                                            'error': trajectories['spbR'].std(axis=0)
                                            })

        plot_data['yaxis']['axis'].append({'data': trajectories['spbL'].mean(axis=0),
                                            'color': 'black',
                                            'error': trajectories['spbL'].std(axis=0)
                                            })

        # User matplotlib to get color gradient
        cm = plt.get_cmap('gist_rainbow')

        for i, ch in enumerate(trajectories['chromosomes']):

            color = cm(1.*i/N)

            plot_data['yaxis']['axis'].append({'data': ch[0].mean(axis=0),
                                   'color': color,
                                   'error': ch[0].std(axis=0)
                                   })
            plot_data['yaxis']['axis'].append({'data': ch[1].mean(axis=0),
                                   'color': color,
                                   'error': ch[1].std(axis=0)
                                   })

        dic_plot(plot_data, os.path.join(eval_results_path, "Kt_Trajectories.svg"))

        return trajectories
