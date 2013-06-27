#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Title: attachment_state
## Description: Here are several calculi of the attachment state evolution
## and the various errors
## Author:uillaume Gay<elagachado AT  gmail DOT com>
## Commentary:

from numpy import *

def get_history(kd):
    """
    """
    meros = hstack([array(ch.mero_history) for ch in kd.chromosomes.values()])
    plugs = hstack([array(ch.pluged_history) for ch in kd.chromosomes.values()])
    return plugs, meros


def balance_histories(kd):

    '''
    input: a kt_simul.spindle_dynamics.KinetochoreDynamics()
    after a simulation instance.

    Returns the difference between the number of correctly
    and erroneously attached plug sites when those two are != 0.
    This is  stored in an 2D array for which each line gives the
    number of kt in each of the possible cases
    (i.e balance = -Mk+2, -Mk+3, .., Mk-2) for each time point.


    '''
    Mk = int(kd.params['Mk'])
    num_steps = len(kd.spbR.traj)

    #The number of cases is (Mk-1) * 2 + 1

    balance = vsplit(zeros((2 * Mk - 1, num_steps)), 2 * Mk - 1)
    num_cuts = zeros(num_steps)

    merotelic_types = {'corrected':zeros(num_steps),
                       'cut':zeros(num_steps),
                       'monotelic':zeros(num_steps),
                       'syntelic':zeros(num_steps)}

    for ch in kd.chromosomes:

        mh = array(ch.mero_history)
        ph = array(ch.pluged_history)

        rbalance = zeros(num_steps)
        lbalance = zeros(num_steps)
        for j, (m, p) in enumerate(zip(mh, ph)):
            if m[0] * p[0] == 0 : #So the kt is not merotelic
                rbalance[j] =  nan
            else:
                rbalance[j] = p[0] - m[0]
            if m[1] * p[1] == 0 :
                lbalance[j] =  nan
            else:
                lbalance[j] = p[1] - m[1]

            if isfinite(lbalance[j]) or isfinite(rbalance[j]):
                if p[0] == m[0]:
                    if p[0] != 0:
                        merotelic_types['cut'][j] += 1
                    else:
                        merotelic_types['monotelic'][j] += 1
                if p[1] == m[1]:
                    if p[1] != 0:
                        merotelic_types['cut'][j] += 1
                    else:
                        merotelic_types['monotelic'][j] += 1

                if (p[1] - m[1]) * (p[0] - m[0]) < 0: #They are pulled the same way
                    merotelic_types['syntelic'][j] += 1

                if (p[1] - m[1]) * (p[0] - m[0]) > 0: #they are pulled appart
                    merotelic_types['corrected'][j] += 1



        for i, bal in enumerate(balance):
            bal += array([rbalance == i - Mk + 2]).flatten().astype(int)
            bal += array([lbalance == i - Mk + 2]).flatten().astype(int)

    return vstack(balance), merotelic_types

def transition_matrix(kd):

    def_tags = {'merotelic':3, 'syntelic':4, 'amphitelic':2,
                    'monotelic':1, 'unattached':0}
    N = int(kd.params['N'])
    num_defects, were_defects = defect_histories(kd)
    trans_mat = zeros((5,5))
    for key, wd in were_defects.items():
        for key2, wd2 in were_defects.items():
            num_trans = 0
            for i in range(N):
                num_trans += (wd.astype(int)[:-1,i]*wd2.astype(int)[1:,i]).sum()
                trans_mat[def_tags[key], def_tags[key2]] += num_trans

    return trans_mat


def defect_histories(kd):
    """
    returns the history of the various attachment states.
    Takes a KinetochoreDynamics instance as unique argument

    output:
    num_defects, were_defects: dectionnaries

    The following outputs are established on a per chromosome
    (i.e. kt pairs) basis

    Those outputs are returned in a dictionnary called num_defects:
    num_defects = {"amphitelic":the number of correctely plugged chromosomes,
               "merotelic":the number of merotelic chromosomes,
               "monotelic":the number of monotelic chromosomes,
               "syntelic":the number of syntelic chromosomes
               "unattached":the number of unpluged chromosomes}
    Normaly num_plugs+num_synt+num_mono+num_mero+num_unplugs = N

    were_plug, were_mero, etc. : 2D arrays of booleans giving the state
    of each chromosome
    """

    Mk = int(kd.params['Mk'])
    N = int(kd.params['N'])

    meros = [array(ch.mero_history) for ch in kd.chromosomes.values()]
    plugs = [array(ch.pluged_history) for ch in kd.chromosomes.values()]

    meros = hstack(meros)
    plugs = hstack(plugs)
    mh_b = meros.astype(bool)
    ph_b = plugs.astype(bool)

    #Put all this in a dictionary
    were_defects = {'amphitelic':were_plug(mh_b, ph_b),
                   'merotelic':were_mero(mh_b, ph_b),
                   'monotelic':were_mono(mh_b, ph_b),
                   'syntelic':were_synt(mh_b, ph_b),
                   'unattached':were_unat(mh_b, ph_b)}

    num_defects = {}

    for key, val in were_defects.items():
        num_defects[key] = val.sum(axis = 1)
    return num_defects, were_defects


def were_plug(mh_b, ph_b):
    """
    Test wether chromosomes were  amphitelic
                Correct  Merotelic
    -------------------------------
         Left  |  True     False
    AND  Right |  True     False
    OR
         Left  |  False    True
    AND  Right |  False    True
    """

    correct_r = ph_b[:,::2] & logical_not(mh_b[:,::2])
    correct_l = ph_b[:,1::2] & logical_not(mh_b[:,1::2])
    and_x = correct_r & correct_l

    mero_r = mh_b[:,::2] & logical_not(ph_b[:,::2])
    mero_l = mh_b[:,1::2] & logical_not(ph_b[:,1::2])
    and_x_mero = mero_r & mero_l

    return (and_x | and_x_mero)


def were_synt(mh_b, ph_b):

    #The number of syntelic chromosomes
    #        Correct  Merotelic OR  Correct  Merotelic
    #------------------------------------------------
    # Left  |  True     False        False    True
    # Right | False     True         True     False

    xor_l = mh_b[:,::2] ^ ph_b[:,::2]   #left
    xor_r = mh_b[:,1::2] ^ ph_b[:,1::2] #right
    xor_x =  ph_b[:,::2] ^ ph_b[:,1::2] #crossed
    xor_t = xor_l & xor_r & xor_x
    return xor_t


def were_mero(mh_b, ph_b):

    #The number of merotelic chromosomes
    and_r = mh_b[:,::2] & ph_b[:,::2]   #          Correct  Merotelic
    and_l = mh_b[:,1::2] & ph_b[:,1::2] #     Left   True     True
    or_x = and_r | and_l                # OR  Right  True     True
    return or_x

#On a phenotype point of view, one should futher distinguish between merotelic chromosomes :
###  * The ones that are amphitelic - merotelic (i.e kts tend to segregate correctly)
###  * The ones that are syntelic - merotelic (both kts go to the same poles)
###  * The ones that are monotelic - merotelic (one kt is merotelic and the other is unattached)
###  * The ones that are cut (one or both kts are equaly attached to both poles)


#So let's go
def were_amphi_mero(mh_b, ph_b):

    and_r = mh_b[:,::2] & ph_b[:,::2]
    and_l = mh_b[:,1::2] & ph_b[:,1::2]

    nor_r = logical_not(mh_b[:,::2] | ph_b[:,::2])
    nor_l = logical_not(mh_b[:,1::2] | ph_b[:,1::2])

    return ( and_r & nor_l ) ^ ( and_l & nor_r )


def were_mono(mh_b, ph_b):

    xor_r = mh_b[:,1::2] ^ ph_b[:,1::2] #right
    xor_l = mh_b[:,::2] ^ ph_b[:,::2]   #left
    #The number of monotelic chromosomes
    nand_r = logical_not(mh_b[:,::2] | ph_b[:,::2])    #            Correct  Merotelic
    nand_l = logical_not(mh_b[:,1::2] | ph_b[:,1::2])   #      Left   False     False
    xor_x =  (nand_r & xor_r)  ^  (nand_l & xor_l)  # XOR  Right  False     False
    return xor_x


def were_unat(mh_b, ph_b):
    #The number of unattached chromosomes
    #             Correct  Merotelic
    # -------------------------------
    #      Left  |  False     False
    # AND  Right |  False     False
    xor_r = mh_b[:,1::2] ^ ph_b[:,1::2] #right
    xor_l = mh_b[:,::2] ^ ph_b[:,::2]   #left
    nand_r = logical_not(mh_b[:,::2] | ph_b[:,::2])
    nand_l = logical_not(mh_b[:,1::2] | ph_b[:,1::2])
    and_x = nand_r & nand_l
    return and_x
