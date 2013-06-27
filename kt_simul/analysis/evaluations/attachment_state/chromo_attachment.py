# -*- coding: utf-8 -*-
import logging

import numpy as np

from kt_simul.analysis.evaluations import Evaluation
from kt_simul.draw.plot import dic_plot

logger = logging.getLogger(__name__)


class ChromoAttachment(Evaluation):
    """
    """

    name = "Chromosome Attachment"
    description = """Return the chromosomes attachment state"""
    group = "attachment_state"
    enable = True

    def __init__(self,):
        pass

    def run(self, KD, draw=False):
        """
        """

        num_steps = KD.spbR.traj.size
        N = int(KD.params['N'])
        Mk = int(KD.params['Mk'])
        ana_onset = int(KD.params["t_A"])

        value_max = N * Mk * 2

        chromo_attach = {'monotelic': np.zeros(num_steps),
                         'syntelic': np.zeros(num_steps),
                         'unattached': np.zeros(num_steps),
                         'amphitelic': np.zeros(num_steps),
                         'merotelic-cut': np.zeros(num_steps),
                         'merotelic-amphi': np.zeros(num_steps),
                         'merotelic-synte': np.zeros(num_steps),
                        }

        for ch in KD.chromosomes:

            iterator = enumerate(zip(ch.correct_history, ch.erroneous_history))
            for i, (correct, incorrect) in iterator:

                # Kt A n° correct attachments
                ac = correct[0]
                # Kt B n° correct attachments
                bc = correct[1]
                # Kt A n° incorrect attachments
                ai = incorrect[0]
                # Kt B n° incorrect attachments
                bi = incorrect[1]

                if ac + bc + ai + bi == 0:
                    # The actual chromosome is unattached (so easy...)
                    chromo_attach['unattached'][i] += 1

                elif (ac + ai + bi == 0 and bc > 0) or \
                     (bc + bi + ai == 0 and ac > 0):
                     # The actual chromosome is monotelic
                    chromo_attach['monotelic'][i] += 1

                elif ai + bi == 0 and ac > 0 and bc > 0:
                    # The actual chromosome is amphitelic
                    chromo_attach['amphitelic'][i] += 1

                elif (ac == 0 and ai > 0) or (bc == 0 and bi > 0):
                    # The actual chromosome is syntelic
                    chromo_attach['syntelic'][i] += 1

                # It should remain here, only merotelic chromosomes.
                # Let's check it
                elif ai > 0 or bi > 0:
                    # Welcome to the merotelic world !

                    if ac == ai or bc == bi:
                        chromo_attach['merotelic-cut'][i] += 1
                    elif ac < ai or bc < bi:
                        chromo_attach['merotelic-synte'][i] += 1
                    elif ac > ai or bc > bi:
                        chromo_attach['merotelic-amphi'][i] += 1
                    else:
                        # The code should never run here unless we got errors !
                        logger.error("Problem in attachment rate algorithm (second level condition")
                        return False

                else:
                    # The code should never run here unless we got errors !
                    logger.error("Problem in attachment rate algorithm (first level condition")
                    return False

        # Normalize values on 1
        for state in chromo_attach:
            chromo_attach[state] = chromo_attach[state] / N

        if draw:
            # Draw attachment state with matplotlib
            timelapse = np.arange(0, KD.duration, KD.dt)

            plot_data = {}
            plot_data['title'] = "Chromosome attachment rate"
            plot_data['xaxis'] = {'data': timelapse, 'label': 'Time'}
            plot_data['yaxis'] = {'label': 'Attachment rate', 'axis': []}
            plot_data['logx'] = True

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

            plot_data['yaxis']['axis'].append({'data': chromo_attach['monotelic'],
                                                'color': 'cyan',
                                                'legend': 'Monotelic'
                                                })

            plot_data['yaxis']['axis'].append({'data': chromo_attach['unattached'],
                                                'color': 'blue',
                                                'legend': 'Unattached'
                                              })

            plot_data['yaxis']['axis'].append({'data': chromo_attach['amphitelic'],
                                                'color': 'green',
                                                'legend': 'Amphitelic'
                                              })

            plot_data['yaxis']['axis'].append({'data': chromo_attach['syntelic'],
                                                'color': 'yellow',
                                                'legend': 'Syntelic'
                                              })

            plot_data['yaxis']['axis'].append({'data': chromo_attach['merotelic-cut'],
                                                'color': 'orange',
                                                'legend': 'Balanced Merotelic'
                                              })

            plot_data['yaxis']['axis'].append({'data': chromo_attach['merotelic-amphi'],
                                                'color': 'red',
                                                'legend': 'Mero-Amphitelic'
                                              })

            plot_data['yaxis']['axis'].append({'data': chromo_attach['merotelic-synte'],
                                               'color': 'purple',
                                               'legend': 'Mero-Syntelic'
                                              })


            dic_plot(plot_data)

        return chromo_attach
