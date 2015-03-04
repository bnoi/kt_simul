# -*- coding: utf-8 -*-
"""
Module dealing with simulation parameters
"""

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import logging
import os

import numpy as np
from scipy.stats import cauchy

from ..io import ParamTree

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
PARAMFILE = os.path.join(ROOT_DIR, 'io', 'params.json')
MEASUREFILE = os.path.join(ROOT_DIR, 'io', 'measures.json')

log = logging.getLogger(__name__)


def get_default_paramtree():
    """
    """
    return ParamTree(PARAMFILE)


def get_default_measuretree():
    """
    """
    return ParamTree(MEASUREFILE, adimentionalized=False)


def reduce_params(paramtree, measuretree, force_parameters=[]):
    """
    This functions changes the parameters so that
    the dynamical characteristics complies with the measures [1]_.

    Parameters
    ----------

    paramtree : :class:`~kt_simul.io.ParamTree` instance
        Modified in place.
    measuretree : :class:`~kt_simul.io.MeasureTree` instance

    References
    ----------
    .. [1] G. Gay, T.Courthéoux, C. Reyes, S. Tournier, Y. Gachet.
           J. Cell Biol 2012 http://dx.doi.org/10.1083/jcb.201107124

    """
    params = paramtree.absolute_dic
    measures = measuretree.absolute_dic

    try:
        poleward_speed = measures['poleward_speed']
        metaph_rate = measures['metaph_rate']
        anaph_rate = measures['anaph_rate']
        mean_metaph_k_dist = measures['mean_metaph_k_dist']
        max_metaph_k_dist = measures['max_metaph_k_dist']
        outer_inner_dist = measures['oi_dist']
        tau_k = measures['tau_k']
        tau_c = measures['tau_c']
        obs_d0 = measures['obs_d0']
        mean_kt_spb_dist = measures["mean_kt_spb_dist"]
        mean_spb_size = measures["mean_spb_size"]

        k_a = params['k_a']
        d_alpha = params['d_alpha']
        N = int(params['N'])
        Mk = int(params['Mk'])
        mus = params['mus']
        Fk = params['Fk']

        ldep_mu = params['ldep_for_attachment_mu']
        ldep_gamma = params['ldep_for_attachment_gamma']
        ldep_N_mt = params['ldep_for_attachment_N_mt']

    except KeyError:
        log.warning("Some parametree / measuretree values are missing")
        return False

    # Set various parameters from measures
    if 'd0' not in force_parameters:
        params['d0'] = obs_d0
    d0 = params['d0']

    if 'ldep_balance' not in force_parameters:
        params['ldep_balance'] = mean_kt_spb_dist

    # Set max speed from measures
    if 'Vk' not in force_parameters:
        params['Vk'] = poleward_speed
    Vk = params['Vk']

    if 'Vmz' not in force_parameters:
        params['Vmz'] = anaph_rate
    Vmz = params['Vmz']

    # Set spring constant from measures
    if 'kappa_c' not in force_parameters:
        params['kappa_c'] = Fk * Mk * 2 / (max_metaph_k_dist - d0)
    kappa_c = params['kappa_c']

    if 'kappa_k' not in force_parameters:
        params['kappa_k'] = Fk * Mk * 2 / (2 * outer_inner_dist)
    kappa_k = params['kappa_k']

    # Set drag coefficient from measure
    if 'muc' not in force_parameters:
        muc = (tau_c * kappa_c)
        params['muc'] = muc

    if 'muk' not in force_parameters:
        muk = (tau_k * kappa_k)
        params['muk'] = muk

    # Get "average" constant of simulation

    # Get mean detachment rate (depend on Aurora parameter)
    if d_alpha != 0:
        # Cauchy distribution
        # x0 = 0
        # prefactor = (2 / (np.pi * d_alpha)) / (1 + ((mean_metaph_k_dist - x0) / (d_alpha / 2)) ** 2)

        # Inverse distribution
        prefactor = d_alpha / (mean_metaph_k_dist / 2)

        k_d_eff = k_a * prefactor
    else:
        k_d_eff = k_a

    # Get mean attachment rate
    if ldep_N_mt > 0:
        x = np.linspace(0, mean_spb_size, 100)
        cauchy_cdf = cauchy.pdf(x, loc=ldep_mu, scale=ldep_gamma)
        total_area = np.trapz(cauchy_cdf, x=x)
        prefactor = cauchy.pdf(mean_kt_spb_dist, loc=ldep_mu, scale=ldep_gamma) / total_area
        prefactor *= ldep_N_mt

        k_a_eff = k_a * prefactor
    else:
        k_a_eff = k_a

    # Compute mean attachment state during a simulation
    alpha_mean = 1 / (1 + k_d_eff / k_a_eff)

    # Compute max force generated by spindle midzone motor
    if 'Fmz' not in force_parameters:
        Fmz = (Fk * N * Mk * alpha_mean * (1 + metaph_rate / (2 * Vk)) + mus * metaph_rate / 2.)
        Fmz /= (1 - metaph_rate / Vmz)
        params['Fmz'] = Fmz

    for key, val in list(params.items()):
        if key not in force_parameters:
            paramtree[key] = val
