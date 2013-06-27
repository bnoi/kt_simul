from  kt_simul.analysis.evaluations import Evaluation

import numpy as np

class PlugedStats(Evaluation):

    name = "PlugedStats"
    description = "PlugedStats"
    group = None
    enable = True

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """

        N = int(KD.params['N'])
        trans_MA = KD.params['t_A']
        dt = KD.params['dt']
        trans_MA = int(trans_MA/dt)
        tot_avg = 0
        delta_avg = 0
        for ch in KD.chromosomes:
            tot = ch.correct_history[:trans_MA].sum(axis=1)
            delta = np.diff(ch.correct_history[:trans_MA], axis=1)
            delta = abs(delta[:trans_MA])
            tot_avg += np.mean(tot)/2
            delta_avg += np.mean(delta)
        tot_avg /= N
        delta_avg /= N
        return (tot_avg, delta_avg)