import numpy as np

from kt_simul.analysis.evaluations import Evaluation
from kt_simul.draw.plot import dic_plot


class KinetoAttachment(Evaluation):
    """
    """

    name = "Kineto Attachment"
    description = """Return the naive attachment state of kinetochore
    (correct, erroneous or unattached)"""
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

        value_max = N * Mk * 2

        kt_attach = {'correct_attached': np.zeros(num_steps, dtype='float'),
                        'incorrect_attached': np.zeros(num_steps, dtype='float'),
                        'unattached': np.zeros(num_steps, dtype='float')
                       }

        for ch in KD.chromosomes:

            ch_incorrect = np.array(ch.erroneous_history).sum(1)
            ch_correct = np.array(ch.correct_history).sum(1)

            kt_attach['correct_attached'] += ch_correct
            kt_attach['incorrect_attached'] += ch_incorrect

        kt_attach['unattached'] = np.ones(num_steps, dtype='float') * value_max
        kt_attach['unattached'] -= kt_attach['correct_attached']
        kt_attach['unattached'] -= kt_attach['incorrect_attached']

        # Scale values on 1
        kt_attach['correct_attached'] /= value_max
        kt_attach['incorrect_attached'] /= value_max
        kt_attach['unattached'] /= value_max

        if draw:
            # Draw attachment state with matplotlib
            timelapse = np.arange(0, KD.duration, KD.dt)

            plot_data = {}
            plot_data['title'] = "Kinetochores attachment rate"
            plot_data['xaxis'] = {'data': timelapse, 'label': 'Time'}
            plot_data['yaxis'] = {'label': 'Attachment rate', 'axis': []}

            # Draw parameters box

            correct_attached_axis = {'data': kt_attach['correct_attached'],
                                     'color': 'green',
                                     'legend': 'Correct attached kinetochore'}

            incorrect_attached_axis = {'data': kt_attach['incorrect_attached'],
                                       'color': 'red',
                                       'legend': 'Incorrect attached kinetochore'}

            unattached_axis = {'data': kt_attach['unattached'],
                               'color': 'blue',
                               'legend': 'Unattached kinetochore'}

            plot_data['yaxis']['axis'].append(correct_attached_axis)
            plot_data['yaxis']['axis'].append(incorrect_attached_axis)
            plot_data['yaxis']['axis'].append(unattached_axis)

            dic_plot(plot_data)

        return kt_attach
