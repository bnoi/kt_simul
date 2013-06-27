import numpy as np
import os
import logging

from kt_simul.analysis.pool_evaluations import PoolEvaluation
from kt_simul.draw.plot import dic_plot

logger = logging.getLogger(__name__)


class MitoticPlate(PoolEvaluation):
    """
    """

    name = "Mitotic Plate"
    description = "Retrieve mitotic plate (Kt dispersion) over time"
    group = "plate"
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

        kt_plate = {'dispersion': np.zeros((nsimu, num_steps)),
                    'plate_shape_right': np.zeros((nsimu, num_steps)),
                    'plate_shape_left': np.zeros((nsimu, num_steps)),
                    'params': params
                   }

        if verbose:
            logger.info("Loading data from simulations files")
        for i, (simu_id, meta) in enumerate(self.iter_simulations(raw_path,
                                                        nsimu=nsimu,
                                                        print_progress=verbose)):
            results = meta.evaluate(name="Mitotic Plate", verbose=False)

            kt_plate['dispersion'][i] = results['dispersion']
            kt_plate['plate_shape_left'][i] = results['plate_shape_left']
            kt_plate['plate_shape_right'][i] = results['plate_shape_right']

        logger.disabled = False

        kt_plate['dispersion_std'] = kt_plate['dispersion'].std(axis=0)
        kt_plate['dispersion'] = kt_plate['dispersion'].mean(axis=0)

        kt_plate['plate_shape_right'] = kt_plate['plate_shape_right'].mean(axis=0)
        kt_plate['plate_shape_right_std'] = kt_plate['plate_shape_right'].std(axis=0)

        kt_plate['plate_shape_left'] = kt_plate['plate_shape_left'].mean(axis=0)
        kt_plate['plate_shape_left_std'] = kt_plate['plate_shape_left'].std(axis=0)

        if not draw:
            return kt_plate

        # Configure and plot the graph
        if verbose:
            logger.info("Plotting results")
        timelapse = meta_infos.timelapse

        plot_data = {}
        plot_data['title'] = "Centromeres dispersion variation"
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time (second)'}
        plot_data['yaxis'] = {'label': 'Dispersion measure (relative to the spindle length)', 'axis': []}
        plot_data['error'] = True
        plot_data['legend'] = False
        plot_data['limit_y_min'] = 0

        # Draw parameters box
        plot_data["params_box"] = [{'name': "Name", 'data': name},
                                   {'name': "Simulations number", 'data': nsimu},
                                   {'name': "Lenght Dependance factor", 'data': params["ldep"]}
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

        plot_data['yaxis']['axis'].append({'data': kt_plate['dispersion'],
                                            'color': 'blue',
                                            'error': kt_plate['dispersion_std']
                                            })

        dic_plot(plot_data, os.path.join(eval_results_path, "Mitotic_Plate_Formation.svg"))

        # Plot the plate shape
        plot_data['title'] = "Mitotic plate formation"
        plot_data['yaxis'] = {'label': "Centromeres position variance", 'axis': []}
        del plot_data['limit_y_min']
        plot_data['legend'] = True

        plot_data['yaxis']['axis'].append({'data': kt_plate['plate_shape_right'],
                                          'color': 'blue',
                                          'legend': 'Right kinetochores',
                                          'error': kt_plate['plate_shape_right_std']
                                          })
        plot_data['yaxis']['axis'].append({'data': kt_plate['plate_shape_left'],
                                          'color': 'red',
                                          'legend': 'Left kinetochores',
                                          'error': kt_plate['plate_shape_left_std']
                                          })

        dic_plot(plot_data, os.path.join(eval_results_path, "Mitotic_Plate_Formation2.svg"))

        return kt_plate
