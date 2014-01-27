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

from kt_simul.io.xml_handler import ParamTree

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
PARAMFILE = os.path.join(ROOT_DIR, 'data', 'params.xml')
MEASUREFILE = os.path.join(ROOT_DIR, 'data', 'measures.xml')
MEASURETREE = ParamTree(MEASUREFILE, adimentionalized=False)
MEASURES = MEASURETREE.absolute_dic

log = logging.getLogger(__name__)


def reduce_params(paramtree, measuretree, force_parameters=[]):
    """
    This functions changes the parameters so that
    the dynamical characteristics complies with the measures [1]_.

    Parameters
    ----------

    paramtree : :class:`~kt_simul.io.xml_handler.ParamTree` instance
        Modified in place.
    measuretree : :class:`~kt_simul.io.xml_handler.MeasureTree` instance

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
    except KeyError:
        log.warning("The measures dictionary should contain"
                       "at least the following keys: ")
        log.warning(list(MEASURES.keys()))
        return False

    k_a = params['k_a']  # 'free' attachement event frequency
    k_d0 = params['k_d0']  # 'free' detachement event frequency
    d_alpha = params['d_alpha']
    N = int(params['N'])
    Mk = int(params['Mk'])
    kappa_k = params['kappa_k']
    Fk = params['Fk']

    #Let's go for the direct relations
    d0 = obs_d0
    if 'd0' not in force_parameters:
        params['d0'] = obs_d0

    Vk = poleward_speed
    if 'Vk' not in force_parameters:
        params['Vk'] = poleward_speed

    Vmz = anaph_rate
    if 'Vmz' not in force_parameters:
        params['Vmz'] = anaph_rate

    if 'ldep_balance' not in force_parameters:
        params['ldep_balance'] = mean_kt_spb_dist

    #Aurora modifies fd
    if d_alpha != 0:
        k_d_eff = k_a * d_alpha / (mean_metaph_k_dist / 2)
    else:
        #log.warning("Things don't go well without Aurora ")
        k_d_eff = k_d0

    # alpha_mean = float(mean_attachment(k_a/fd_eff) / Mk)
    alpha_mean = 1 / (1 + k_d_eff / k_a)
    # Take metaphase kt pair distance as the maximum one
    # TODO : kc = Fk * Mt * alpha_mean / (max_metaph_k_dist - d0)
    kappa_c = Fk * Mk * 2 / (max_metaph_k_dist - d0)
    if 'kappa_c' not in force_parameters:
        params['kappa_c'] = kappa_c

    #kop = alpha_mean * ( 1 + metaph_rate/2 ) / ( outer_inner_dist )
    kappa_k = Fk * Mk * 2 / (2 * outer_inner_dist)
    if 'kappa_k' not in force_parameters:
        params['kappa_k'] = kappa_k
    #Ensure we have sufficientely small time steps
    dt = params['dt']
    # params['dt'] = min(tau_c / 4., tau_k / 4., params['dt'])
    # if params['dt'] != dt:
    #     log.info('Time step changed')

    mus = params['mus']
    Fmz = (Fk * N * Mk * alpha_mean * (1 + metaph_rate / (2 * Vk))
           + mus * metaph_rate / 2.) / (1 - metaph_rate / Vmz)
    if 'Fmz' not in force_parameters:
        params['Fmz'] = Fmz
    muc = (tau_c * kappa_c)
    if 'muc' not in force_parameters:
        params['muc'] = muc
    muk = (tau_k * kappa_k)
    if 'muk' not in force_parameters:
        params['muk'] = muk

    for key, val in list(params.items()):
        if key not in force_parameters:
            paramtree.change_dic(key, val, verbose=False)
