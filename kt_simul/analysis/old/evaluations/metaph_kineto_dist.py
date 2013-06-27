from  kt_simul.analysis.evaluations import Evaluation

import numpy as np

class MetaphaseKinetoDistance(Evaluation):

    name = "Metaphase Kineto Distance"
    description = "Calculate the mean and sdev of the distance between the kinetochores of each chromosome during metaphase. Returns a tuple (mean, sdev)"
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
        average = 0
        stdev = 0
        for n in range(N):
            kinetoR = KD.chromosomes[n].cen_A.traj
            kinetoL = KD.chromosomes[n].cen_B.traj
            dist = kinetoR - kinetoL
            dist = dist[:trans_MA]
            average += np.mean(dist)
            stdev += np.std(dist)
        average /= N
        stdev /= N
        return average, stdev