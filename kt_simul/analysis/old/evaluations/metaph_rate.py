from  kt_simul.analysis.evaluations import Evaluation

import numpy as np

class MetaphaseRate(Evaluation):

    name = "Metaphase rate"
    description = "Calculates metaphase elongation rate."
    group = None
    enable = True

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """

        N = int(KD.params['N'])
        t_A = KD.params['t_A']
        dt = KD.params['dt']
        spindle_length = KD.spbR.traj - KD.spbL.traj
        if spindle_length.size >= int(t_A/dt):
            elapsed = np.arange(t_A, step = dt)
            (a, b)=np.polyfit(elapsed, spindle_length[:elapsed.shape[0]], 1)
        else:
            elapsed = np.arange(spindle_length.size*dt, step = dt)
            (a, b)=np.polyfit(elapsed, spindle_length, 1)

        return a

