import logging
log = logging.getLogger(__name__)

import numpy as np
import pandas as pd

__all__ = ["get_fft", "get_fft_downsampled", "filter_signal"]


def get_fft(x, dt):
    """
    """
    n = len(x)

    fft_output = np.fft.rfft(x)     # Perform real fft
    rfreqs = np.fft.rfftfreq(n, dt) # Calculatel frequency bins
    fft_mag = np.abs(fft_output)    # Take only magnitude of spectrum

    # Normalize the amplitude by number of bins and multiply by 2
    # because we removed second half of spectrum above the Nyqist frequency
    # and energy must be preserved
    fft_mag = fft_mag * 2 / n

    return np.array(fft_mag), np.array(rfreqs)


def get_fft_downsampled(x, dt, new_dt, verbose=False):
    """
    See http://www.dspguide.com/ch9/1.htm
    """

    n = x.size
    n_jump = new_dt / dt
    n_segments = n / n_jump

    if verbose:
        log.info("Numver of segments : {}".format(n_segments))
        log.info("Space between each element in one segment : {}".format(n_jump))
        log.info("New dt : {} s".format(new_dt))
        log.info("Max frequency : {} Hz".format(1 / new_dt))

    all_fft_mag = []
    all_rfreqs = []

    for i in np.arange(n_segments):
        new_x = x[i::n_jump]
        _fft_mag, _rfreqs = get_fft(new_x, new_dt)
        all_fft_mag.append(_fft_mag)
        all_rfreqs.append(_rfreqs)

    all_fft_mag = pd.DataFrame(all_fft_mag).values
    all_rfreqs = pd.DataFrame(all_rfreqs).values

    # Remove columns with Nan values
    all_fft_mag = np.ma.compress_rows(np.ma.fix_invalid(all_fft_mag.T)).T
    all_rfreqs = np.ma.compress_rows(np.ma.fix_invalid(all_rfreqs.T)).T

    # Get mean values along columns
    mean_fft_mag = all_fft_mag.mean(axis=0)
    mean_freqs = all_rfreqs.mean(axis=0)

    return mean_fft_mag, mean_freqs


def filter_signal(x, dt, filter_value):
    n = len(x)

    fft_output = np.fft.rfft(x)

    rfreqs = np.fft.rfftfreq(n, d=dt)

    # Filtering
    fft_filtered = fft_output.copy()
    fft_filtered[(rfreqs > filter_value)] = 0

    # Inversed TF
    x_filtered = np.fft.irfft(fft_filtered)

    times = np.arange(0, n, dt)
    if x_filtered.size != times.size:
        times_filtered = times[:-1]
    else:
        times_filtered = times

    return x_filtered, times_filtered
