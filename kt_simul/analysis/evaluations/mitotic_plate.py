# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from kt_simul.analysis.evaluations import Evaluation
from kt_simul.draw.plot import dic_plot


class MitoticPlate(Evaluation):
    """
    """

    name = "Mitotic Plate"
    description = "Retrieve mitotic plate (Kt dispersion) over time"
    group = "plate"
    enable = True

    def __init__(self,):
        pass

    def run(self, KD, draw=False):
        """
        """

        N = int(KD.params["N"])
        ana_onset = int(KD.params["t_A"])

        timelapse = np.arange(0, KD.duration, KD.dt)
        num_steps = KD.duration / KD.dt

        kt_traj_step = num_steps / 20

        spindle_length = abs(KD.spbR.traj - KD.spbL.traj)

        # Dispersion is size of the segment made by extrem Kt scaled
        # by spindle size
        kt_plate = {'dispersion': None,
                    'kt_trajs': [],
                    'plate_shape_right': np.zeros((N , num_steps), dtype="float"),
                    'plate_shape_left': np.zeros((N , num_steps), dtype="float")
                    }

        max_kt = np.zeros((N * 2, num_steps), dtype="float")
        min_kt = np.zeros((N * 2, num_steps), dtype="float")

        kt_trajs = []

        for i, ch in enumerate(KD.chromosomes):
            id1 = (i * 2) - 1
            id2 = (i * 2)

            max_kt[id1] = ch.cen_A.traj
            max_kt[id2] = ch.cen_B.traj

            min_kt[id1] = ch.cen_A.traj
            min_kt[id2] = ch.cen_B.traj

            kt_trajs.append(ch.cen_A.traj / spindle_length)
            kt_trajs.append(ch.cen_B.traj / spindle_length)

            if ch.cen_A.traj.mean() < ch.cen_B.traj.mean():
                kt_plate['plate_shape_left'][i] = ch.cen_A.traj
                kt_plate['plate_shape_right'][i] = ch.cen_B.traj
            else:
                kt_plate['plate_shape_right'][i] = ch.cen_A.traj
                kt_plate['plate_shape_left'][i] = ch.cen_B.traj

        kt_plate['plate_shape_left'] = kt_plate['plate_shape_left'].std(axis=0)
        kt_plate['plate_shape_right'] = kt_plate['plate_shape_right'].std(axis=0)

        max_kt_distance = max_kt.max(axis=0)
        min_kt_distance = min_kt.min(axis=0)

        kt_plate['dispersion'] = (abs(max_kt_distance - min_kt_distance) /
                                    spindle_length)

        # Scale kt trajectories
        for kt in kt_trajs:
            kt = kt + kt_plate['dispersion']

            # Lighten array
            light_traj = np.zeros(num_steps, dtype='float')
            light_traj[::kt_traj_step] = 1
            light_traj *= kt
            light_traj[light_traj == 0] = np.nan

            kt_plate['kt_trajs'].append(light_traj)

        if not draw:
            return kt_plate

        # Draw attachment state with matplotlib
        timelapse = np.arange(0, KD.duration, KD.dt)

        plot_data = {}
        plot_data['title'] = "Centromeres dispersion variation"
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time (second)'}
        plot_data['yaxis'] = {'label': "Dispersion measure (relative to the spindle length)", 'axis': []}
        plot_data['logx'] = False
        plot_data['legend'] = False
        plot_data['limit_y_min'] = 0

        # Draw parameters box
        plot_data["params_box"] = [{'name': "Lenght Dependance factor", 'data': KD.params["ldep"]},
                                  ]

        # Add annotation about anaphase onset
        plot_data["annotations"] = []
        plot_data["annotations"].append({'s': 'Anaphase onset: %i' % ana_onset,
                                       'xy': (ana_onset, 0),
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

        plot_data['yaxis']['axis'].append({'data': kt_plate['dispersion'],
                                          'color': 'blue',
                                          'plot_args' : {'lw': 2}
                                          })

        dic_plot(plot_data)

        # User matplotlib to get color gradient
        cm = plt.get_cmap('gist_rainbow')
        for j, kt in enumerate(kt_plate['kt_trajs']):

            if j%2 == 0:
                color = cm(1. * j / (N))

            plot_data['yaxis']['axis'].append({'data': kt,
                                               'color': color,
                                               'plot_args' : {'marker': 'o'}
                                               })

        dic_plot(plot_data)

        # Plot the plate shape
        plot_data['title'] = "Mitotic plate formation"
        plot_data['yaxis'] = {'label': "Centromeres position variance", 'axis': []}
        del plot_data['limit_y_min']

        plot_data['yaxis']['axis'].append({'data': kt_plate['plate_shape_right'],
                                          'color': 'blue',
                                          'legend': 'Right kinetochores'
                                          })
        plot_data['yaxis']['axis'].append({'data': kt_plate['plate_shape_left'],
                                          'color': 'red',
                                          'legend': 'Left kinetochores'
                                          })

        dic_plot(plot_data)

        return kt_plate
