# -*- coding: utf-8 -*-

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

from sympy import symbols, Matrix
from sympy import sympify
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics import ReferenceFrame, Point, Particle
from sympy.physics.mechanics import LagrangesMethod, Lagrangian

from sympy.physics.mechanics import mechanics_printing
mechanics_printing(pretty_print=True) # Shinny

from collections import defaultdict



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

        ## dictionary of lists
        self.forcebalance = defaultdict(list)

        ## Coordinates
        ### independant
        self.q_ind = []
        self.qd_ind = []

        self.u_ind = [] #Kanes
        self.ud_ind = []

        ### dependant
        self.q_dep = []
        self.qd_dep = []

        self.u_dep = []
        self.ud_dep = []
        self._setup_done = False

    @property
    def lagrangian(self):
        return sympify(sum(self.forcebalance['lagrangian']))

    @property
    def q(self):
        return self.forcebalance['q_ind'] + self.forcebalance['q_dep']
    @property
    def qd(self):
        return self.forcebalance['q_ind'] + self.forcebalance['q_dep']
    @property
    def u(self):
        return self.forcebalance['q_ind'] + self.forcebalance['q_dep']
    @property
    def ud(self):
        return self.forcebalance['q_ind'] + self.forcebalance['q_dep']

class ViscousDrag:
    ''' Simple viscous drag

    '''
    def __init__(self, qd, mu, point, e_F):
        self.point = point
        self.F = - mu * qd * e_F

    @classmethod
    def new(cls, spindle, *args, **kwargs):
        '''
        Parameters
        ----------
        spindle: a :class:`kt_simul.SympySpinde` instance
        qd : symbol
          the generalized speed
        mu : the friction coefficient
        point : point of aplication of the viscous drag
        e_F : force unit vector
        '''
        viscous = ViscousDrag(*args, **kwargs)
        spindle.forcebalance['forces'].append((viscous.point,
                                               viscous.F))
        return viscous

class DampedSpring:
    '''
    A Damped spring
    ===============

    a damped spring is governed by its deformation length - assumed
    corected for equilibrium distance and the time diferentiate of this deformation

    Attributes
    ----------

    TODO
    '''

    def __init__(self, length, lengthd,
                 kappa, mu,
                 point1, point2, e_F):

        self.length = length
        self.lengthd = lengthd
        self.F_k = kappa * self.length * e_F
        self.F_v = - mu * self.lengthd * e_F
        self.point1 = point1
        self.point2 = point2
        self.lagrangian = kappa * self.length**2 / 2

    @classmethod
    def new(cls, spindle, *args, **kwargs):
        '''Adds an over-damped spring in `spindle`

        Parameters
        ----------
        spindle: a :class:`kt_simul.SympySpinde` instance
        length : dynamicsymbol
          the extension of the spring corrected for eq. distance.
        lengthd : dynamicsymbol
          the damped speed
        kappa : Young modulus
        mu : friction coefficient
        poitn1, point2: extremities of the spring
        e_F : unit vector allong the force axis
        '''

        spring = DampedSpring(*args, **kwargs)
        spindle.forcebalance['forces'].extend([(spring.point1, spring.F_k),
                                               (spring.point2, -spring.F_k),
                                               (spring.point1, spring.F_v),
                                               (spring.point2, -spring.F_v)])

        spindle.forcebalance['lagrangian'].append(spring.lagrangian)
        return spring

class LinearFV():
    '''
    Linear Force velocity relationship

    In order to use methods such as Lagrange or Kanes, one must
    get rid of constant terms in force expressions. This is achieved
    through the definition of constrained generalized speeds.
    In this class, both the generalized speed and the force vector
    are computed.

    Attributes
    ----------

    ell : dependant generalized coordinate
    v_g : dependant generalized speed
    coord_constraint: relation between the generalized coordinate
      and the input positions
    speed_constraint : relation between the generalized speed
      and the input parameters
    F : force vector, as a function of the generalized speed and coordinate

    We have
    ..math:
        F = F_{max} ( 1 - \gamma \frac{\dot{l}}{V_{max}})
    and
    ..math:
        F = \frac{F_{max} v_g}{V_{max}} \mathbf{e}_F

    Example
    -------
    >>> ### container
    >>> expl_spindle = SympySpindle('expl')

    >>> ### coordinates
    >>>  q1, q2 = dynamicsymbols("q1, q2")
    >>> q1d, q2d = dynamicsymbols("q1d, q2d")
    >>> N = ReferenceFrame('N')
    >>> ### points
    >>> O = Point('O')
    >>> O.set_vel(N, 0)
    >>> p1 = O.locatenew('p1', q1 * N.x)
    >>> p1.set_vel(N, q1d * N.x)
    >>> p2 = O.locatenew('p2', q2 * N.x)
    >>> p2.set_vel(N, q2d * N.x)
    >>> ## Motor Parameters
    >>> F_m, V_m = symbols("F_m, V_m")
    >>> ### Add to the container
    >>> lfv = LinearFV.new(expl_spindle,
    ...                    p1, p2, F_m, V_m, gamma=1,
    ...                    e_F=N.x, index='o',  expo='ex')
    >>> print(lfv.coord_constraint,  lfv.speed_constraint)
    >>> print(lfv.F)
    >>> assert lfv.ell.diff() == lfv.elld
    '''

    def __init__(self, point1, point2, F_max,
                 V_max, gamma, e_F,
                 index=None, expo=None):

        '''
        point1, point2: points
            between which the force is applied.
        F_max: symbol,
          the stall force
        V_max: symbol,
         the maximum speed (in intended regime)
        gamma: int
         if gamma = 1, force tends to extend the distance ell
         if gamma = 0, force tends to reduce ell
        e_F : unit vector
         the unit vector along which the force is applied
         (should probably be  in the reference frame used to
          compute ld)



        '''

        self.point1 = point1
        self.point2 = point2

        ## generalized coordinate
        ell_name = indexer('ell', index, expo)
        self.ell = dynamicsymbols(ell_name)
        self.elld = dynamicsymbols(ell_name, 1)
        ## coordinate constraint
        self.coord_constraint  = self.ell - (point1.pos_from(point2) & e_F)

        ## generalized speed
        self.v_g = dynamicsymbols(indexer('v_g', index, expo))
        ## speed constraint (v_g = V_max - gamma * ld)
        self.speed_constraint = self.v_g + gamma * self.elld - V_max # == 0

        ## Force
        self.F = (F_max / V_max) * self.v_g * e_F

        self.F_ = self.F.subs(self.v_g, V_max - gamma * self.elld)

    @classmethod
    def new(cls, spindle, *args, **kwargs):


        lfv = cls(*args, **kwargs)
        ## Register
        spindle.forcebalance['forces'].extend([(lfv.point1, lfv.F),
                                               (lfv.point2, - lfv.F)])

        spindle.q_dep.append(lfv.ell)
        spindle.qd_dep.append(lfv.elld)
        spindle.u_dep.append(lfv.v_g)
        spindle.forcebalance['kd'].append(lfv.elld - lfv.ell)
        spindle.forcebalance['configuration_constraints'].append(lfv.coord_constraint)
        spindle.forcebalance['speed_constraints'].append(lfv.speed_constraint)
        return lfv

def indexer(name, index=None, expo=None):
    '''Small utility to append an index  and
    an exponent to the string `name`
    '''
    if index is not None:
        name = name + '_' + index
    if expo is not None:
        name = name + '^' + expo
    return name
