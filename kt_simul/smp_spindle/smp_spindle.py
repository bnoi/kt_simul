# -*- coding: utf-8 -*-

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

from sympy import symbols, Matrix
from sympy import sympify
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics import ReferenceFrame, Point, Particle



from sympy.physics.mechanics import mechanics_printing
mechanics_printing(pretty_print=True) # Shinny

from .utils import indexer


class SympySpindle:
    '''
    container class for the implementation of a formal
    or symbolic model in sympy.

    We collect all we instanciate at the sympy level

    Attribues:
    ----------
    TODO


    '''

    def __init__(self, name):
        '''
        A mitotic spindle in sympy.

        Contains the generalized coordinates and speeds,
        kinematic diff. equations and a `forcebalance`
        dictionnary of lists containing the dynamical
        aspects (e.g. particles & forces). Please refer to
        the documentation (for now, notebooks here:
        http://nbviewer.ipython.org/github/glyg/Notebooks/blob/master/kt_sim)
        for further details.

        Once all elements are collected, the spindle is fed to
        a :class:`sympy.physics.mechanics.KanesMethod` object.
        '''
        self.name = name

        self.forces = []
        ## Coordinates
        self.q_ind = []
        self.qd_ind = []
        self.u_ind = []

        ## points
        self.points = []
        ### Kinematic diff
        self.kd = []
        self._setup_done = False
