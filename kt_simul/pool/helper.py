import os
import logging
import gc

import numpy as np
from scipy import signal

from ..signal.fft import get_fft
from . import Pool


def simu_analyzer(simu, params):

    d = {}
    d.update(params.to_dict())

    # Load trajs
    times = simu.time

    anaphase = simu.time_anaphase
    index_anaphase = simu.index_anaphase

    # Get mean spindle length in metaphase
    spindle_length = np.abs(simu.KD.spbL.traj - simu.KD.spbR.traj)
    d['spindle_size_metaphase'] = spindle_length[:index_anaphase].mean()

    # Get spindle length at anaphase onset
    d['spindle_size_ana_onset'] = spindle_length[index_anaphase]

    # Get spindle elongation rate
    d['spindle_elongation'] = np.abs(spindle_length[index_anaphase] - spindle_length[0]) / anaphase

    d['stretch_mean_metaphase'] = []
    d['stretch_mean_ana_onset'] = []
    d['kt_pos_mean_ana_onset'] = []
    d['kt_pos_mean_metaphase'] = []
    d['kt_pos_mean_ana_onset_relative'] = []
    d['kt_pos_mean_metaphase_relative'] = []
    for i in range(1, 3):
        d["magn_{}".format(i)] = []
        d["period_{}".format(i)] = []

    spindle_size = np.abs(simu.KD.spbL.traj - simu.KD.spbR.traj)

    for ch in simu.KD.chromosomes:

        # Get mean stretch in metaphase
        d['stretch_mean_metaphase'].append(np.abs(ch.cen_A.traj - ch.cen_B.traj)[:index_anaphase].mean())

        # Get stretch at anaphase onset
        d['stretch_mean_ana_onset'].append(np.abs(ch.cen_A.traj - ch.cen_B.traj)[index_anaphase])

        # Get mean kt psoition during metaphase
        kt_traj = np.abs((ch.cen_A.traj + ch.cen_B.traj) / 2)

        d['kt_pos_mean_ana_onset'].append(kt_traj[index_anaphase])
        d['kt_pos_mean_metaphase'].append(kt_traj[:index_anaphase].mean())
        d['kt_pos_mean_ana_onset_relative'].append((kt_traj / spindle_size)[index_anaphase])
        d['kt_pos_mean_metaphase_relative'].append((kt_traj / spindle_size)[:index_anaphase].mean())

        # Amplitude and Period of oscillations in metaphase
        ch_traj = (ch.cen_A.traj + ch.cen_B.traj) / 2
        fft_mag, rfreqs = get_fft(ch_traj[400:], simu.KD.dt)

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
                d["magn_{}".format(i)].append(fft_mag[peak_idx] * 2)
                d["period_{}".format(i)].append(1 / rfreqs[peak_idx])

    d['stretch_mean_metaphase'] = np.mean(d['stretch_mean_metaphase'])
    d['stretch_mean_ana_onset'] = np.mean(d['stretch_mean_ana_onset'])

    d['kt_pos_mean_ana_onset'] = np.mean(d['kt_pos_mean_ana_onset'])
    d['kt_pos_mean_metaphase'] = np.mean(d['kt_pos_mean_metaphase'])
    d['kt_pos_mean_ana_onset_relative'] = np.mean(d['kt_pos_mean_ana_onset_relative'])
    d['kt_pos_mean_metaphase_relative'] = np.mean(d['kt_pos_mean_metaphase_relative'])

    for i in range(1, 3):
        d["magn_{}".format(i)] = np.mean(d["magn_{}".format(i)])
        d["period_{}".format(i)] = np.mean(d["period_{}".format(i)])

    # Attachment defect at anaphase
    atts = simu.get_attachment_vector()
    d['attachment_defect_ana_onset'] = np.sum(atts[index_anaphase] != 0) / atts[index_anaphase].shape[0]

    d['time_anaphase'] = simu.time_anaphase

    # Get aneuploidy rate
    # TODO

    return d


def run_simu(i, params, simu_path, paramtree, measuretree, pool_parameters):
    """
    """

    simu_name = os.path.join(simu_path, 'simu_{}.h5'.format(i))

    current_paramtree = paramtree.copy()
    current_measuretree = measuretree.copy()

    for k, v in params.iteritems():
        current_paramtree[k] = v

    pool_params = {'simu_path': simu_name,
                   'load': False,
                   'n_simu': pool_parameters['N'],
                   'paramtree': current_paramtree,
                   'measuretree': current_measuretree,
                   'initial_plug': pool_parameters['initial_plug'],
                   'parallel': False,
                   'verbose': False,
                   'erase': False}

    pool = Pool(**pool_params)
    pool.run()

    return pool


def processor(i, params, simu_path, paramtree, measuretree, pool_parameters, n_params):

    # Log
    real_i = (i + 1) * pool_parameters["N"]
    logging.info("Processing simulations: {}/{}".format(real_i, n_params))

    # Run pool of simus
    pool = run_simu(i, params, simu_path, paramtree, measuretree, pool_parameters)

    data = []

    # Load and analyze simus
    for j, simu in enumerate(pool.load_metaphases()):

        d = simu_analyzer(simu, params)
        del simu
        gc.collect()

        data.append(d)

    os.remove(pool.simu_path)

    del pool
    gc.collect()

    return data
