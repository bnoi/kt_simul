#!/usr/bin/python
# -*- coding: utf-8 -*-

from numpy import *
import sys, os
import time

import pyximport
pyximport.install()

from kt_simul.core.simul_spindle import *
from kt_simul.core.xml_handler import *
from kt_simul.analysis.attachment_state import *

#You can change the path to the default parameter and measure files here
try:
    paramfile = os.path.join(os.path.dirname(__file__),
                             'default', 'params.xml')
    measurefile = os.path.join(os.path.dirname(__file__),
                               'default', 'measures.xml')

except NameError: #Not importing from a module
    paramfile = 'default/params.xml'
    measurefile = 'default/measures.xml'


def full_simul(new_params, plug = 'monotelic'):

    '''
    this functions runs a kinetochore dynamics simulation after having changed
    the parameters specified in the dictionnary new_params. plug is the initial
    attachment state of the kinetochores.

    accepted values for plug are:
    'monotelic', 'merotelic', 'null' or 'random'
    
    '''
    m = Metaphase()
    for key, new_value in new_params.items(): 
        m.paramtree.change_dic(key, new_value, write = False,
                               back_up = False, verbose = False)
    #we don't want to execute anaphase during this kind of simulation
    m.paramtree.change_dic('t_A', num_steps, write = False,
                           back_up = False, verbose = False)
    m.__init__(paramtree = m.paramtree,  plug = plug)
    m.simul()
    return m
    

def explore_2D(pcs1s, auroras, num_steps, num_ech, plug, dt = 1):

    new_params = {}
    logfile = file('%s.log' %base_name, 'w+')

    for n, pcs1 in enumerate(pcs1s):
        for m, d_alpha in enumerate(auroras):
            new_params['d_alpha'] = d_alpha
            new_params['orientation'] = pcs1
            new_params['span'] = num_steps * dt
            new_params['dt'] = dt

            for i in range(num_ech):
                #run the simulation
                mp = full_simul(new_params, plug = plug)
                #save results
                xmlfname = '%s_res_pcs1-%03i_auroras-%03i_%03i.xml' %(base_name, m, i)
                datafname = '%s_data_pcs1-%03i_auroras-%03i_%03i.npy' %(base_name, m, i)
                mp.write_results(xmlfname, datafname)
                del mp

            logfile.write('ran plug = %s, pcs1 = %03f, d_alpha = %03f\n' %(plug, pcs1, d_alpha))
    return 0
 
if __name__ == "__main__":

    if len(sys.argv) != 2 :
        print 'Usage: python batch_explore <base_directory/simulations_prefix>'
    base_name = sys.argv[1]
    t0 = time.time()

    
    #Generate an arrray of values for the parameters we want
    # to explore
    #Here the orientation and aurora parameters
    pcs1s = linspace(0., 1., 31)
    auroras = logspace(-3., .2, 51)

    #Save the explored values in separate files for later
    pcs1sname = '%s_pcs1s.npy' %base_name
    save(pcs1sname, pcs1s)

    auroraname = '%s_auroras.npy' %base_name
    save(auroraname, auroras)

    #Set some conditions
    plug = 'random'
    num_steps = 800
    dt = 1.
    num_ech = 5
    #run the simulations and save the results (see above)
    explore_2D(pcs1s, auroras, num_steps, num_ech, plug, dt)
    print 'time: %.3f' %(time.time() - t0)

