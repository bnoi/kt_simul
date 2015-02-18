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
from .forces import ViscousDrag, DampedSpring, LinearFV



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
        '''A mitotic spindle in sympy.

        Contains the generalized coordinates and speeds, kinematic
        diff. equations and lists containing the dynamical aspects
        (e.g. particles & forces). Please refer to the documentation
        (for now, notebooks here:
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
        self.S = ReferenceFrame('S')
        self.center = Point('C')
        self.center.set_vel(self.S, 0)
        self.length = dynamicsymbols('L')
        self.elong_rate = dynamicsymbols('vL')
        self._setup_done = False


class Organite:

    def __init__(self, name, parent,
                 ref, pos, vel):

        self.name = name
        self.ref = ref
        self.point = parent.locatenew(name, pos)
        self.point.set_vel(ref, vel)
        self.forces = []

    def viscous(self, mu):

        viscous = ViscousDrag(self.ref,
                              self.point, mu)
        self.forces.append((self.point, viscous.F))


class Spb(Organite):

    def __init__(self, side, spindle):

        name = 'spbR' if side == 1 else 'spbL'
        pos = side * spindle.length * spindle.S.x / 2
        vel = side * spindle.elong_rate * spindle.S.x / 2
        self.side = side
        self.spindle = spindle
        Organite.__init__(name, spindle.center, spindle.S, pos, vel)
