from  kt_simul.analysis.evaluations import Evaluation

import numpy as np

class AnaphaseRate(Evaluation):
    """
    Calculates anaphase B elongation rate.
    """

    name = "Anaphase rate"
    description = "Calculates anaphase B elongation rate."
    group = None
    enable = True

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """

        t_A = KD.params['t_A']
        dt = KD.params['dt']
        trans_AB = self.anaph_transAB(KD)
        stop = len(KD.spbR.traj)*dt
        if not trans_AB:
            return 0
        else:
            spindle_length = np.array(KD.spbR.traj) - np.array(KD.spbL.traj)
            elapsed = np.arange(trans_AB, stop, dt)
            (a, b)=np.polyfit(elapsed, spindle_length[int(trans_AB/dt):], 1)
        return a

    def anaph_transAB(self, KD):
        """
        As anaphase A -> B transition is defined (at least in the lab)
        by the date at which the first kinetochore reaches and stay at the spb
        we have to retrieve this moment. Fortunately, we already know
        times of arrival
        """

        toas = []
        for ch in KD.chromosomes:
            toas.append(ch.cen_A.toa)
            toas.append(ch.cen_B.toa)
        trans_AB = min(toas)
        return trans_AB
