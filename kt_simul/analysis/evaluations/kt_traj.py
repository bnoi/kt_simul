# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from kt_simul.analysis.evaluations import Evaluation
from kt_simul.draw.plot import dic_plot


class KtTrajectories(Evaluation):
    """
    """

    name = "Kinetochores Trajectories"
    description = """Make profile of the simulation drawing kinetochores trajectories"""
    group = "trajectories"
    enable = True

    def __init__(self,):
        pass

    def run(self, KD, draw=False):
        """
        """

        N = int(KD.params["N"])
        ana_onset = int(KD.params["t_A"])

        timelapse = np.arange(0, KD.duration, KD.dt)

        trajectories = {'spbR': KD.spbR.traj,
                        'spbL': KD.spbL.traj,
                        'chromosomes': [],
                       }

        # Retrieve kinetochores trajectories
        for ch in KD.chromosomes:

            traj = [ch.cen_A.traj,
                    ch.cen_B.traj]

            trajectories['chromosomes'].append(traj)

        if not draw:
            return trajectories

        # Draw attachment state with matplotlib
        timelapse = np.arange(0, KD.duration, KD.dt)

        plot_data = {}
        plot_data['title'] = "Chromosome kinetochores trajectories"
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time'}
        plot_data['yaxis'] = {'label': 'Relative position to SPB axes center', 'axis': []}
        plot_data['logx'] = False
        plot_data['legend'] = False

        # Add annotation about anaphase onset
        plot_data["annotations"] = []
        plot_data["annotations"].append({'s': 'Anaphase onset: %i' % ana_onset,
                                       'xy': (ana_onset, -2),
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

        plot_data['yaxis']['axis'].append({'data': trajectories['spbR'],
                                          'color': 'black',
                                          'plot_args': {'ls': '-',
                                                        'lw': 3}
                                          })

        plot_data['yaxis']['axis'].append({'data': trajectories['spbL'],
                                          'color': 'black',
                                          'plot_args': {'ls': '-',
                                                        'lw': 3}
                                          })

        # User matplotlib to get color gradient
        cm = plt.get_cmap('gist_rainbow')

        for i, ch in enumerate(trajectories['chromosomes']):

            color = cm(1. * i / N)

            plot_data['yaxis']['axis'].append({'data': ch[0],
                                   'color': color,
                                   })
            plot_data['yaxis']['axis'].append({'data': ch[1],
                                   'color': color,
                                   })

        dic_plot(plot_data)

        return trajectories
