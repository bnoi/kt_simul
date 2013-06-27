from  kt_simul.analysis.evaluations import Evaluation

import numpy as np

class PolewardSpeed(Evaluation):

    name = "Poleward Speed"
    description = "Poleward Speed"
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
        start = int(trans_MA/dt)

        trans_MA = KD.params['t_A']
        dt = KD.params['dt']
        trans_AB = self.anaph_transAB(KD)

        pole_speeds = []
        for ch in KD.chromosomes:
            r_stop = int(ch.cen_A.toa/dt)
            if r_stop - start > 2:
                r_dist = - KD.spbR.traj[start:r_stop] + ch.cen_A.traj[start:r_stop]
                elapsed = np.r_[trans_MA:ch.cen_A.toa:dt]
                (ra,rb) = np.polyfit(elapsed, r_dist, 1)
                pole_speeds.append(ra)
            l_stop = int(ch.cen_B.toa/dt)
            if l_stop - start > 2:
                l_dist = KD.spbL.traj[start:l_stop] - ch.cen_B.traj[start:l_stop]
                elapsed = np.r_[trans_MA:ch.cen_B.toa:dt]
                (la,lb) = np.polyfit(elapsed, l_dist, 1)
                pole_speeds.append(la)

        pole_speeds = np.array(pole_speeds)

        return pole_speeds.mean(), pole_speeds.std()

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