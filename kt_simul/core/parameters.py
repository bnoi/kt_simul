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
from scipy import linalg

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

        k_d0 = params['k_d0']
        k_a = params['k_a']
        d_alpha = params['d_alpha']
        N = int(params['N'])
        Mk = int(params['Mk'])
        mus = params['mus']
        Fk = params['Fk']

    except KeyError:
        log.warning("Some parametree / measuretree values are missing")
        return False

    # Three states
    k_res_max = params['k_res']
    k_cat_max = params['k_cat']
    vg = params['vg']
    vs = params['vs']
    mean_kt_v = measures['mean_kt_v']

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

    if 'k_d0' not in force_parameters:
        params['k_d0'] = k_a

    # Get "average" constant of simulation
    # Compute average detachment rate
    if d_alpha != 0:
        k_d_eff = k_a * d_alpha / (mean_metaph_k_dist / 2)
    else:
        k_d_eff = k_d0

    ### Compute average probability to be in attached state

    ## Two states version
    alpha_mean = k_a / (k_a + k_d_eff)

    ## Three states naive version
    alpha_mean = (k_a / (k_a + k_d_eff)) #* (k_cat_max / (k_cat_max + k_res_max))

    # Three states markov chain version
    mean_k_cat = k_cat_max #/ (1 + np.exp((mean_kt_v - vg) / vg))
    mean_k_res = k_res_max #/ (1 + np.exp((mean_kt_v - vs) / vs))

    # We define the transition matrix P as follow (see doc for more details)
    # q = [p_GA, p_SA, p_D]
    P = np.array([[1 - mean_k_cat, mean_k_cat,                 0],
                  [mean_k_res,     1 - (mean_k_res + k_d_eff), k_d_eff],
                  [k_a,            0,                          1 - k_a]])

    # We want to find the steady state of the transition matrix P
    # We make the hypothesis that at steady state q does not depend on P anymore
    # q * P = q
    w, v = linalg.eig(P, left=True, right=False)

    scale_unit_index = np.where(np.isclose(w, 1))[0]
    if len(scale_unit_index) < 1:
        raise Exception("Can't find correct eigen vector, fallback to naive three state Markov chain version")

    v = v[:, scale_unit_index].real
    v /= v.sum()

    # Average probability of being attached and producing a force is GA
    #alpha_mean = v[0, 0]

    log.debug("Probability of being attached and producing a force is {}".format(alpha_mean))

    # Compute max force generated by spindle midzone motor
    if 'Fmz' not in force_parameters:
        Fmz = (Fk * N * Mk * alpha_mean * (1 + metaph_rate / (2 * Vk)) + mus * metaph_rate / 2.)
        Fmz /= (1 - metaph_rate / Vmz)
        params['Fmz'] = Fmz

    for key, val in list(params.items()):
        if key not in force_parameters:
            paramtree[key] = val
