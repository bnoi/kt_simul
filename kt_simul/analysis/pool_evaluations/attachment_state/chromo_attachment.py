import numpy as np
import os
import logging

from kt_simul.analysis.pool_evaluations import PoolEvaluation
from kt_simul.draw.plot import dic_plot

logger = logging.getLogger(__name__)


class ChromoAttachment(PoolEvaluation):
    """
    """

    name = "Chromosome Attachment"
    description = """Process and plot mean and std of the returned values of
    Evaluation called 'Chromosome Attachment'"""
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

        chromo_attach = {'monotelic': np.zeros((nsimu, num_steps)),
                         'syntelic': np.zeros((nsimu, num_steps)),
                         'unattached': np.zeros((nsimu, num_steps)),
                         'amphitelic': np.zeros((nsimu, num_steps)),
                         'merotelic-cut': np.zeros((nsimu, num_steps)),
                         'merotelic-amphi': np.zeros((nsimu, num_steps)),
                         'merotelic-synte': np.zeros((nsimu, num_steps)),
                        }

        if verbose:
            logger.info("Loading data from simulations files")
        for i, (simu_id, meta) in enumerate(self.iter_simulations(raw_path,
                                                        nsimu=nsimu,
                                                        print_progress=verbose)):
            results = meta.evaluate(name="Chromosome Attachment", verbose=False)

            for attach_name, value in results.items():
                chromo_attach[attach_name][i] = value

        logger.disabled = False

        if not draw:
            return chromo_attach

        # Configure and plot the graph
        if verbose:
            logger.info("Plotting results")
        timelapse = meta_infos.timelapse

        plot_data = {}
        plot_data['title'] = "Chromosome attachment rate"
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time'}
        plot_data['yaxis'] = {'label': 'Chromosome rate', 'axis': []}
        plot_data['error'] = False

        # Draw parameters box
        plot_data["params_box"] = [{'name': "Name", 'data': name},
                                   {'name': "Simulations number", 'data': nsimu},
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

        plot_data['yaxis']['axis'].append({'data': chromo_attach['monotelic'].mean(axis=0),
                                            'color': 'cyan',
                                            'legend': 'Monotelic',
                                            'error': chromo_attach['monotelic'].std(axis=0)
                                            })

        plot_data['yaxis']['axis'].append({'data': chromo_attach['unattached'].mean(axis=0),
                                            'color': 'blue',
                                            'legend': 'Unattached',
                                            'error': chromo_attach['unattached'].std(axis=0)
                                          })

        # plot_data['yaxis']['axis'].append({'data': chromo_attach['amphitelic'].mean(axis=0),
        #                                     'color': 'green',
        #                                     'legend': 'Amphitelic',
        #                                     'error': chromo_attach['amphitelic'].std(axis=0)
        #                                   })

        plot_data['yaxis']['axis'].append({'data': chromo_attach['syntelic'].mean(axis=0),
                                            'color': 'yellow',
                                            'legend': 'Syntelic',
                                            'error': chromo_attach['syntelic'].std(axis=0)
                                          })

        plot_data['yaxis']['axis'].append({'data': chromo_attach['merotelic-cut'].mean(axis=0),
                                            'color': 'orange',
                                            'legend': 'Balanced Merotelic',
                                            'error': chromo_attach['merotelic-cut'].std(axis=0)
                                          })

        plot_data['yaxis']['axis'].append({'data': chromo_attach['merotelic-amphi'].mean(axis=0),
                                            'color': 'red',
                                            'legend': 'Mero-Amphitelic',
                                            'error': chromo_attach['merotelic-amphi'].std(axis=0)
                                          })

        plot_data['yaxis']['axis'].append({'data': chromo_attach['merotelic-synte'].mean(axis=0),
                                           'color': 'purple',
                                           'legend': 'Mero-Syntelic',
                                            'error': chromo_attach['merotelic-synte'].std(axis=0)
                                          })

        dic_plot(plot_data, os.path.join(eval_results_path, "Chromosome_Attachment.svg"))

        # Now plot with logarithmic xaxis
        plot_data['xaxis']['label'] = 'Time (Logarithmic scale)'
        plot_data['logx'] = True
        dic_plot(plot_data, os.path.join(eval_results_path, "Chromosome_Attachment_logarithmic.svg"))

        # Draw same plot but without all merotelic type : we merge them
        plot_data['logx'] = False
        plot_data['yaxis']['axis'] = plot_data['yaxis']['axis'][0:4]

        mero = chromo_attach['merotelic-cut'] + chromo_attach['merotelic-synte'] + chromo_attach['merotelic-amphi']

        plot_data['yaxis']['axis'].append({'data': mero.mean(axis=0),
                                   'color': 'red',
                                   'legend': 'Merotelic',
                                    'error': mero.std(axis=0)
                                  })
        dic_plot(plot_data, os.path.join(eval_results_path, "Chromosome_Attachment_Merge_Mero.svg"))

        # Now plot with logarithmic xaxis
        plot_data['xaxis']['label'] = 'Time (Logarithmic scale)'
        plot_data['logx'] = True
        dic_plot(plot_data, os.path.join(eval_results_path, "Chromosome_Attachment_Merge_Mero_logarithmic.svg"))


        return chromo_attach
