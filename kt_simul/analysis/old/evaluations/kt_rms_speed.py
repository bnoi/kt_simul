from  kt_simul.analysis.evaluations import Evaluation

import numpy as np

class KtRmsSpeed(Evaluation):
    """
    KtRmsSpeed
    """

    name = "KtRmsSpeed"
    description = "KtRmsSpeed"
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
        stop = int(trans_MA/dt)
        rms_speeds = []

        for ch in KD.chromosomes:
            r_speed = np.diff(ch.cen_A.traj[:stop])
            r_rmss = np.sqrt((r_speed**2).mean())
            rms_speeds.append(r_rmss)

            l_speed = np.diff(ch.cen_B.traj[:stop])
            l_rmss = np.sqrt((l_speed**2).mean())
            rms_speeds.append(l_rmss)

        rms_speeds = np.array(rms_speeds)

        return rms_speeds.mean(), rms_speeds.std()

