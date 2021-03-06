import numpy as np
from scipy import signal

from ..signal.fft import get_fft


def simu_analyzer(simu, params):

    results = []

    times = simu.time
    anaphase = simu.time_anaphase
    index_anaphase = simu.index_anaphase

    spindle_size = np.abs(simu.KD.spbL.traj - simu.KD.spbR.traj)

    for j, ch in enumerate(simu.KD.chromosomes):

        d = {}
        d.update(params.to_dict())
        d['ch_id'] = j

        d['spindle_size_metaphase'] = spindle_size[:index_anaphase].mean()
        d['spindle_size_ana_onset'] = spindle_size[index_anaphase]
        d['spindle_elongation'] = np.abs(spindle_size[index_anaphase] - spindle_size[0]) / anaphase

        d['stretch_mean_metaphase'] = np.abs(ch.cen_A.traj - ch.cen_B.traj)[:index_anaphase].mean()
        d['stretch_mean_ana_onset'] = np.abs(ch.cen_A.traj - ch.cen_B.traj)[index_anaphase]

        kt_traj = np.abs((ch.cen_A.traj + ch.cen_B.traj) / 2)

        d['kt_pos_mean_ana_onset'] = kt_traj[index_anaphase]
        d['kt_pos_mean_metaphase'] = kt_traj[:index_anaphase].mean()
        d['kt_pos_mean_ana_onset_relative'] = (kt_traj / spindle_size)[index_anaphase]
        d['kt_pos_mean_metaphase_relative'] = (kt_traj / spindle_size)[:index_anaphase].mean()

        fft_mag, rfreqs = get_fft(kt_traj[400:], simu.KD.dt)

        # Get peaks
        max_peaks_idxs = signal.argrelmax(fft_mag)[0]

        # Remove peaks which are greater than half the entire duration
        if len(max_peaks_idxs) > 0:
            first_max_peak_id = max_peaks_idxs[0]
            peak_duration = 1 / rfreqs[first_max_peak_id]
            event_half_duration = times.max() / 2
            if peak_duration > event_half_duration:
                fft_mag = np.delete(fft_mag, first_max_peak_id)
                rfreqs = np.delete(rfreqs, first_max_peak_id)

            # Set freq 0hZ to -np.inf
            fft_mag[0] = -np.inf

            # Sort peaks by magnitude
            sorted_idxs = np.argsort(fft_mag[max_peaks_idxs])[::-1]
            max_magn_idxs = max_peaks_idxs[sorted_idxs]

            for i, peak_idx in zip(range(1, 3), max_magn_idxs):
                d["magn_{}".format(i)] = fft_mag[peak_idx] * 2
                d["period_{}".format(i)] = 1 / rfreqs[peak_idx]

        # Attachment defect at anaphase
        atts = simu.get_attachment_vector()
        atts[atts != 0] = 1
        d['attachment_defect_ana_onset'] = atts[index_anaphase].mean()

        results.append(d)

    return results


def attach_analyzer(simu, params):
    """
    """

    results = []

    d = {}
    d.update(params.to_dict())

    d['attachments'] = simu.get_attachment_vector()
    d['times'] = simu.time

    results.append(d)

    return results
