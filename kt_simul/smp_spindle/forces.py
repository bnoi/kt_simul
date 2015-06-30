# -*- coding: utf-8 -*-

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

class ViscousDrag:
    ''' Simple viscous drag

    '''
    def __init__(self, N, point, mu):
        '''
        '''
        self.point = point
        self.F = - mu * point.vel(N)

    @classmethod
    def new(cls, spindle, *args, **kwargs):
        '''
        Parameters
        ----------
        spindle: a :class:`kt_simul.SympySpinde` instance
        u : symbol
          the generalized speed
        mu : the friction coefficient
        point : point of aplication of the viscous drag
        e_F : force unit vector
        '''
        viscous = ViscousDrag(*args, **kwargs)
        spindle.forces.append((viscous.point,
                               viscous.F))
        return viscous

class DampedSpring:
    '''
    A Damped spring
    ===============

    a damped spring is governed by its deformation length - assumed
    corrected for equilibrium distance and the time diferentiate of this
    deformation

    Attributes
    ----------

    TODO
    '''

    def __init__(self, N, point1, point2,
                 kappa, mu, d_eq=0, e_F=None):

        self.r = point2.pos_from(point1)
        if e_F is None:
            self.e_F = self.r.normalize()
        else:
            self.e_F = e_F
        self.F_k = - kappa *  (self.r - d_eq * self.e_F)
        self.F_v = - mu * (point2.vel(N) - point1.vel(N))
        self.point1 = point1
        self.point2 = point2

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
        spindle.forces.extend([(spring.point1, -spring.F_k),
                               (spring.point2, spring.F_k),
                               (spring.point1, -spring.F_v),
                               (spring.point2, spring.F_v)])
        return spring

class LinearFV:
    '''
    Linear Force velocity relationship
    F = F_max(1 - gamma * (u2 - u1)/V_max)e_F

        p1      p2
          .____.
       <--     -->
        -F      F gamma = 1

        -->   <--
        -F      F gamma = -1

    Attributes
    ----------

    F : force vector
    v: expression
         the relative speed of point2 wrt to point 1

    We have
    ..math:
        F = \gamma F_{max} ( 1 - \frac{v}{V_{max}})

    '''

    def __init__(self, N, point1, point2,
                 F_max, V_max, gamma, e_F=None):

        '''
        point1, point2: points
            between which the force is applied.
        F_max: symbol,
          the stall force
        V_max: symbol,
         the maximum speed (in intended regime)
        gamma: int
         if gamma = 1, force tends to extend the distance ell
         if gamma = -1, force tends to reduce ell
        '''
        self.point1 = point1
        self.point2 = point2
        self.v = point2.vel(N) - point1.vel(N)
        self.r = point2.pos_from(point1)
        if e_F is None:
            self.e_F = self.r.normalize()
        else:
            self.e_F = e_F
        ### Force
        self.F = gamma * F_max * (self.e_F - gamma * self.v / V_max)

    @classmethod
    def new(cls, spindle, *args, **kwargs):
        '''
        point1, point2: points
            between which the force is applied.
        F_max: symbol,
          the stall force
        V_max: symbol,
         the maximum speed (in intended regime)
        gamma: int
         if gamma = 1, force tends to extend the distance ell
         if gamma = -1, force tends to reduce ell
        e_F : unit vector
         the unit vector along which the force is applied
         if e_F is None, it will be computed as the unit vector of
         the relative speed
        '''

        lfv = cls(*args, **kwargs)
        ## Register
        # if extending, force goes outward the segment,
        # thus is negative for point1
        spindle.forces.extend([(lfv.point1, -lfv.F),
                               (lfv.point2, lfv.F)])
        return lfv
