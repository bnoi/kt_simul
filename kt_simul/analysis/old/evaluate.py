# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
from numpy.fft import rfft, fftfreq
from numpy.random import normal
from scipy.interpolate import splrep, splev

import pyximport
pyximport.install()

from kt_simul.core.xml_handler import *
from kt_simul.core.spindle_dynamics import *

__all__ = ["evaluations", "anaphase_rate", "metaph_rate",
           "metaph_kineto_dist", "auto_corel", "poleward_speed",
           "kt_rms_speed", "time_of_arrival", "pluged_stats"]

# Launch a bunch of function to evaluate
# and analyse the simulation

def evaluations():
    function_dic = {'anaphase_rate': anaphase_rate,
                    'metaph_rate': metaph_rate,
                    'mean_metaph_k_dist': metaph_kineto_dist,
                   # 'pitch' : auto_corel,
                    'poleward_speed': poleward_speed,
                    'kt_rms_speed': kt_rms_speed,
                    'times_of_arrival': time_of_arrival,
                    'pluged_stats': pluged_stats}
    return function_dic

def get_kt_speeds(KD, step = 2):

    speeds = []
    for ch in KD.chromosomes.values():
        speeds.append(np.diff(ch.cen_A.traj[::step])/step)
        speeds.append(np.diff(ch.cen_B.traj[::step])/step)

    return speeds

def get_spb_speed(KD, step = 2):

    speed = np.diff(KD.spbR.traj[::step])/step
    return speed


def anaph_transAB(KD):
    '''
    As anaphase A -> B transition is defined (at least in the lab)
    by the date at which the first kinetochore reaches and stay at the spb
    we have to retrieve this moment. Fortunately, we already know
    times of arrival
    '''
    toas = []
    for ch in KD.chromosomes:
        toas.append(ch.cen_A.toa)
        toas.append(ch.cen_B.toa)
    trans_AB = min(toas)
    return trans_AB

def anaphase_rate(KD):
    '''
    Calculates anaphase B elongation rate.
    '''
    t_A = KD.params['t_A']
    dt = KD.params['dt']
    trans_AB = anaph_transAB(KD)
    stop = len(KD.spbR.traj)*dt
    if not trans_AB:
        return 0
    else:
        spindle_length = np.array(KD.spbR.traj) - np.array(KD.spbL.traj)
        elapsed = np.arange(trans_AB, stop, dt)
        (a, b)=np.polyfit(elapsed, spindle_length[int(trans_AB/dt):], 1)
    return a

def metaph_rate(KD):
    '''Calculates metaphase elongation rate.
    '''
    N = int(KD.params['N'])
    t_A = KD.params['t_A']
    dt = KD.params['dt']
    spindle_length = KD.spbR.traj - KD.spbL.traj
    if spindle_length.size >= int(t_A/dt):
        elapsed = np.arange(t_A, step = dt)
        (a, b)=np.polyfit(elapsed, spindle_length[:elapsed.shape[0]], 1)
    else:
        elapsed = np.arange(spindle_length.size*dt, step = dt)
        (a, b)=np.polyfit(elapsed, spindle_length, 1)

    return a


def metaph_kineto_dist(KD):
    '''calculates the mean and sdev of the distance between the kinetochores
    of each chromosome during metaphase. Returns a tuple (mean, sdev)
    '''
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

def poleward_speed(KD):

    N = int(KD.params['N'])
    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    start = int(trans_MA/dt)

    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    trans_AB = anaph_transAB(KD)

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

def kt_rms_speed(KD):

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

def time_of_arrival(KD):

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

def pluged_stats(KD):

    N = int(KD.params['N'])
    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    trans_MA = int(trans_MA/dt)
    tot_avg = 0
    delta_avg = 0
    for ch in KD.chromosomes:
        tot = ch.correct_history[:trans_MA].sum(axis=1)
        delta = np.diff(ch.correct_history[:trans_MA], axis=1)
        delta = abs(delta[:trans_MA])
        tot_avg += np.mean(tot)/2
        delta_avg += np.mean(delta)
    tot_avg /= N
    delta_avg /= N
    return (tot_avg, delta_avg)

def auto_corel(KD, smooth = 10.):
    '''
    Calculates the "pitch" of chromosomes"trajectory (kind of a
    Pitch detection algorithm
    '''

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


def max_freqs(KD, show_fig = True):

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


def fwhm(wave):
    '''
    Returns full width at half maximum
    '''

    m = wave.max()
    x = wave.argmax()

    right = left = nan

    wave /= m
    for i in range(min(wave[x:-1].size, wave[1:x].size)):
        if wave[x+i] <= 0.5 < wave[x+i+1]:
            right = i
        if wave[x-i] <= 0.5 < wave[x-i-1]:
            left = i

        if isfinite(right) and isfinite(left):
            return right + left


def first_min(wave):
    '''
    returns the first minimum of an array, starting from the first element
    if interp = True, evaluates the minimum position via linear interpolation
    of np.diff(wave)

    '''

    for i in range(1,wave[1:-1].size):
        if wave[i-1] > wave[i] < wave[i+1]:
            return i
    return None
