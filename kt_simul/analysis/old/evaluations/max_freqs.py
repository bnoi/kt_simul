from  kt_simul.analysis.evaluations import Evaluation

import numpy as np
from scipy.interpolate import splrep, splev
from numpy.fft import rfft, fftfreq

class MaxFreqs(Evaluation):

    name = "MaxFreqs"
    description = "MaxFreqs"
    group = None
    enable = True

    def __init__(self,):
        pass

    def run(self, KD, show_fig = False):
        """
        """

        trans_MA = KD.params['t_A']
        dt = KD.params['dt']
        smth = int(10./dt)
        elapsed = np.arange(0, trans_MA , dt)
        if show_fig : figure(100)
        max_fs = []
        for ch in KD.chromosomes:
            center = (ch.cen_A.traj + ch.cen_B.traj)/2
            center_metaphase = center[:elapsed.shape[0]]
            knots = elapsed[smth:-smth:smth]
            tck = splrep(elapsed, center_metaphase, t = knots)
            sp_speed = splev(elapsed, tck, der = 1)
            sp_fft = rfft(sp_speed)
            freqs = fftfreq(elapsed.shape[0], dt)
            m_freq = freqs[np.argmax(sp_fft)]
            max_fs.append(m_freq)
            if show_fig:
                plot(freqs[:len(sp_fft)], abs(sp_fft), 'o', alpha = 0.6)

        return max_fs