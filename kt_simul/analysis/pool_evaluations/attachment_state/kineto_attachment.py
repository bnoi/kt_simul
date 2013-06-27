import numpy as np
import os
import logging

from kt_simul.analysis.pool_evaluations import PoolEvaluation
from kt_simul.draw.plot import dic_plot

logger = logging.getLogger(__name__)


class KinetoAttachment(PoolEvaluation):
    """
    """

    name = "Kineto Attachment"
    description = """Process and plot mean and std of the returned values of
    Evaluation called 'Kineto Attachment'"""
    group = "attachment_state"
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
        Mk = int(params["Mk"])

        attach_state = {'correct_attached': np.zeros((nsimu, num_steps)),
                        'incorrect_attached': np.zeros((nsimu, num_steps)),
                        'unattached': np.zeros((nsimu, num_steps))
                       }

        if verbose:
            logger.info("Loading data from simulations files")
        for i, (simu_id, meta) in enumerate(self.iter_simulations(raw_path,
                                                        nsimu=nsimu,
                                                        print_progress=verbose)):
            results = meta.evaluate(name="Kineto Attachment", verbose=False)

            attach_state['correct_attached'][i] = results['correct_attached']
            attach_state['incorrect_attached'][i] = results['incorrect_attached']
            attach_state['unattached'][i] = results['unattached']

        logger.disabled = False

        # Mean data
        logger.info("Processing data")
        correct_attached = attach_state['correct_attached'].mean(axis=0)
        correct_attached_std = attach_state['correct_attached'].std(axis=0)
        incorrect_attached = attach_state['incorrect_attached'].mean(axis=0)
        incorrect_attached_std = attach_state['incorrect_attached'].std(axis=0)
        unattached = attach_state['unattached'].mean(axis=0)
        unattached_std = attach_state['unattached'].std(axis=0)

        if not draw:
            return attach_state

        # Configure and plot the graph
        if verbose:
            logger.info("Plotting results")
        timelapse = meta_infos.timelapse

        plot_data = {}
        plot_data['title'] = """Kinetochores attachment rate"""
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time'}
        plot_data['yaxis'] = {'label': 'Attachment rate', 'axis': []}

        # Draw parameters box
        plot_data["params_box"] = [{'name': "Name", 'data': name},
                                   {'name': "Simulations number", 'data': nsimu},
                                   {'name': "Chromosome number", 'data': '%i' % N},
                                   {'name': "Plug site per kinetochore", 'data': '%i' % Mk}
                             ]

        # Add annotation about anaphase onset
        plot_data["annotations"] = []
        plot_data["annotations"].append({'s': 'Anaphase onset: %i' % ana_onset,
                                         'xy': (ana_onset, 0),
                                         'xytext': (0,50),
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

        correct_attached_axis = {'data': correct_attached,
                                 'color': 'green',
                                 'legend': 'Correct attached kinetochore',
                                 'error': correct_attached_std,
                                 }

        incorrect_attached_axis = {'data': incorrect_attached,
                                   'color': 'red',
                                   'legend': 'Incorrect attached kinetochore',
                                   'error': incorrect_attached_std,
                                   }

        unattached_axis = {'data': unattached,
                           'color': 'blue',
                           'legend': 'Unattached kinetochore',
                           'error': unattached_std,
                           }

        plot_data['yaxis']['axis'].append(correct_attached_axis)
        plot_data['yaxis']['axis'].append(incorrect_attached_axis)
        plot_data['yaxis']['axis'].append(unattached_axis)

        dic_plot(plot_data, os.path.join(eval_results_path, "Kineto_Attachment.svg"))

        # Now plot with log xaxis
        plot_data['xaxis']['label'] = 'Time (Logarithmic scale)'
        plot_data['logx'] = True
        dic_plot(plot_data, os.path.join(eval_results_path, "Kineto_Attachment_logarithmic.svg"))

        return attach_state
