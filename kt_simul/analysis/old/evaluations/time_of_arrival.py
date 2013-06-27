from  kt_simul.analysis.evaluations import Evaluation

import numpy as np

class TimeOfArrival(Evaluation):

    name = "TimeOfArrival"
    description = "TimeOfArrival"
    group = None
    enable = True

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """

        toas = []
        for ch in KD.chromosomes:
            toas.append(ch.cen_A.toa)
            toas.append(ch.cen_B.toa)

        toas = np.array(toas)
        if any(toas):
            toas -= toas[toas>0].min()
        else :
            toas -= 1.
        return toas


