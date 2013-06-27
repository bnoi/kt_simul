from  kt_simul.analysis.evaluations import Evaluation

import numpy as np

class AutoCorel(Evaluation):
    """
    Calculates the "pitch" of chromosomes"trajectory (kind of a Pitch detection algorithm
    """

    name = "AutoCorel"
    description = 'Calculates the "pitch" of chromosomes"trajectory (kind of a Pitch detection algorithm'
    group = None
    enable = False

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """
        trans_MA = KD.params['t_A']
        dt = KD.params['dt']
        pitches = []
        if len(KD.spbR.traj) <= trans_MA/dt: #no anaphase execution
            elapsed = np.arange(len(KD.spbR.traj)*dt, step = dt)
        else:
            elapsed = np.arange(0, trans_MA, dt)
        smth = int(smooth/dt)

        for ch in KD.chromosomes:

            cen_A = ch.cen_A.traj[:elapsed.shape[0]]
            # In order to compare with the 'real world' tracked kinetochores
            # we smooth by a factor smooth/dt
            cen_A_tck = splrep(elapsed, cen_A, t = elapsed[smth:-smth:smth])
            cen_A_s = splev(elapsed, cen_A_tck, der = 1)
            m_speed, st_speed = cen_A_s.mean(), cen_A_s.std()
            cen_A_sc = (cen_A_s - m_speed)/st_speed
            co_cen_A = np.correlate(cen_A_sc, cen_A_sc, 'full') / cen_A_sc.size
            pitches.append(first_min(co_cen_A[-co_cen_A.size//2:]))

            cen_B = ch.cen_B.traj[:elapsed.shape[0]]
            cen_B_tck = splrep(elapsed, cen_B, t = elapsed[smth:-smth:smth])
            cen_B_s = splev(elapsed, cen_B_tck, der = 1)
            m_speed, st_speed = cen_B_s.mean(), cen_B_s.std()
            cen_B_sc = (cen_B_s - m_speed)/st_speed
            co_cen_B = np.correlate(cen_B_sc, cen_B_sc, 'full') / cen_B_sc.size
            pitches.append(first_min(co_ktL[-co_ktL.size//2:]))
        try:
            pitches = 1/(np.array(pitches)*dt)
        except:
            return 0
        return pitches.mean(), pitches.std()