# -*- coding: utf-8 -*-

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

from sympy import symbols, Matrix
from sympy import sympify, lambdify
from sympy.physics.vector import cross
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics import ReferenceFrame, Point, Particle
from sympy.physics.mechanics import KanesMethod

from sympy.physics.mechanics import mechanics_printing

from .utils import indexer
from .forces import ViscousDrag, DampedSpring, LinearFV

mechanics_printing(pretty_print=True)  # Shinny

param_names = ['mu_s', 'mu_ch',
               'kappa_c', 'd_0', 'mu_c',
               'kappa_k', 'mu_k',
               'F_mz', 'V_mz', 'F_k', 'V_k']

parameters = {name: symbols(name)
              for name in param_names}


class SympySpindle:
    '''
    container class for the implementation of a formal
    or symbolic model in sympy.

    We collect all we instanciate at the sympy level

    Attribues
    ---------

    TODO


    '''

    def __init__(self, name, N=1, Mk=1):
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
        self.N = N
        self.Mk = Mk

        self.points = []
        self.S = ReferenceFrame('S')
        self.forces = []
        self.torques = []

        # ## Coordinates
        self.q_ind = []
        self.qd_ind = []
        self.u_ind = []

        # ## points
        self.center = Point('C')
        self.center.set_vel(self.S, 0)

        self.length = dynamicsymbols('L')
        self.elong_rate = dynamicsymbols('vL')
        self.q_ind.append(self.length)
        self.u_ind.append(self.elong_rate)
        self._setup_done = False
        self.add_spbs()
        self.chromosomes = []
        for n in range(self.N):
            self.add_chromosome(n)
        self.qd_ind = [coord.diff()
                       for coord in self.q_ind]

    def add_spbs(self):
        self.spbR = Spb(1, self)
        self.spbL = Spb(-1, self)
        F_mz = parameters['F_mz']
        V_mz = parameters['V_mz']
        LinearFV.new(self, self.S,
                     self.spbL.point, self.spbR.point,
                     F_mz, V_mz, 1, e_F=self.S.x)

    def add_chromosome(self, n):
        ch = Chromosome(self, n)
        self.chromosomes.append(ch)

    @property
    def kd(self):
        return [-u + qd for (u, qd) in zip(self.u_ind, self.qd_ind)]

    def kanesmethod(self):
        self.KM = KanesMethod(self.S, self.q_ind, self.u_ind, self.kd)
        particles = [Particle('p{}'.format(point.name), point, 0)
                     for point in self.points]
        self.eoms, frstar = self.KM.kanes_equations(self.forces+self.torques,
                                                    particles)

    @property
    def A_uu(self):
        return self.eoms.jacobian(Matrix(self.u_ind))

    @property
    def B_qq(self):
        return self.eoms.jacobian(Matrix(self.q_ind))

    @property
    def C(self):
        u_zero = {u: 0 for u in self.u_ind}
        q_zero = {q: 0 for q in self.q_ind}
        return self.eoms.subs(u_zero).subs(q_zero)

    @property
    def attach_state(self):
        astate = Matrix([[site_A.lbda, site_A.rho, site_B.lbda, site_B.rho]
                         for ch in self.chromosomes for (site_A, site_B)
                         in zip(ch.cen_A.plugsites, ch.cen_B.plugsites)
                         ]).reshape(2 * self.N * self.Mk, 2)
        return astate

    def get_np_A0(self):
        adim = {parameters['F_k']: 1, parameters['V_k']: 1}
        passive_params = ['mu_s', 'mu_ch', 'mu_c', 'mu_k',
                          'kappa_c', 'kappa_k', 'd_0']
        passive_zero = {parameters[name]: 0 for name in passive_params}
        A0 = self.A_uu.subs(passive_zero).subs(adim)
        A0_args = [parameters[name] for name in passive_params]
        A0_np = lambdify(A0_args, A0, 'numpy', dummify=False)
        return A0_args, A0_np

    def get_np_At(self):
        adim = {parameters['F_k']: 1, parameters['V_k']: 1}
        active_params = ['F_k', 'V_k', 'F_mz', 'V_mz']
        active_zero = {parameters[name]: 0 for name in active_params}
        At = self.A_uu.subs(active_zero).subs(adim)
        At_args = [parameters[name]
                   for name in active_params] + list(self.attach_state)
        At_np = lambdify(At_args, At, 'numpy', dummify=False)
        return At_np, At_args


class Organite:
    ''' Creates an organite

    Attributes
    ----------
    name: str,
        the organite's name
    point: `:class:sympy.mechanics.Point`
    ref: `:class:sympy.mechanics.ReferenceFrame`
        in which the point's velocity is defined
    '''

    def __init__(self, name, parent,
                 ref, pos, vel):

        self.name = name
        self.ref = ref
        self.point = parent.locatenew(name, pos)
        self.point.set_vel(ref, vel)


def gen_3d_coords(spindle, letter, ref):

    coords = [dynamicsymbols("{}_{}".format(u, letter))
              for u in 'xyz']
    spindle.q_ind.extend(coords)

    coords_vel = [dynamicsymbols("v{}_{}".format(u, letter))
                  for u in 'xyz']

    spindle.u_ind.extend(coords_vel)
    pos = (coords[0] * ref.x +
           coords[1] * ref.y +
           coords[2] * ref.z)
    vel = (coords_vel[0] * ref.x +
           coords_vel[1] * ref.y +
           coords_vel[2] * ref.z)
    return pos, vel


class Spb(Organite):
    ''' A spindle pole body aka Microtubule Organizing Center aka Centrosome
    '''
    def __init__(self, side, spindle):

        name = 'spbR' if side == 1 else 'spbL'
        pos = side * spindle.length * spindle.S.x / 2
        vel = side * spindle.elong_rate * spindle.S.x / 2
        self.side = side
        self.spindle = spindle
        Organite.__init__(self, name, spindle.center, spindle.S, pos, vel)
        self.spindle.points.append(self.point)
        mu_s = parameters['mu_s']
        viscous = ViscousDrag.new(self.spindle, self.spindle.S,
                                  self.point, mu_s)


class Chromosome(Organite):
    '''
    A pair of sister chromatids
    '''

    def __init__(self, spindle, ch_id):
        name = 'ch_{}'.format(ch_id)
        self.id = ch_id
        self.spindle = spindle
        self.add_centromeres()

    def add_centromeres(self):
        '''
        Creates the centromeres
        '''
        kappa_c = parameters['kappa_c']
        d_0 = parameters['d_0']
        mu_c = parameters['mu_c']
        self.cen_A = Centromere(self.spindle, self, 'A')
        self.cen_B = Centromere(self.spindle, self, 'B')
        cen_spring = DampedSpring.new(self.spindle,
                                      self.spindle.S,
                                      self.cen_A.point,
                                      self.cen_B.point,
                                      kappa_c, mu_c, d_0)


class Centromere(Organite):

    def __init__(self, spindle, chromosome, tag):
        self.spindle = spindle
        self.tag = tag
        self.id = chromosome.id
        self.chromosome = chromosome
        name = 'cen_{}{}'.format(chromosome.id, tag)
        letter = '{}{}'.format(chromosome.id, tag)
        # sign = -1 if tag == 'A' else 1
        # pos = sign * (chromosome.strech
        #          + parameters['d_0'])* chromosome.C.x / 2
        # vel = sign * chromosome.strech_vel * chromosome.C.x /2
        pos, vel = gen_3d_coords(spindle, letter, spindle.S)
        Organite.__init__(self, name, spindle.center, spindle.S, pos, vel)
        self.spindle.points.append(self.point)
        mu_ch = parameters['mu_ch']

        viscous = ViscousDrag.new(self.spindle,
                                  self.spindle.S,
                                  self.point, mu_ch)
        self.plugsites = []
        for m in range(self.spindle.Mk):
            self.add_attachsite(m)

    def add_attachsite(self, site_id):

        site = PlugSite(self, site_id)
        self.plugsites.append(site)
        kappa_k = parameters['kappa_k']
        mu_k = parameters['mu_k']
        spring = DampedSpring.new(self.spindle, self.spindle.S,
                                  self.point, site.point,
                                  kappa_k, mu_k)


class PlugSite(Organite):

    def __init__(self, centromere, site_id):

        self.spindle = centromere.spindle
        self.site_id = (centromere.id, site_id,
                        centromere.tag)
        name = 'as_{}{}{}'.format(*self.site_id)
        letter = '{}{}{}'.format(*self.site_id)
        self.centromere = centromere
        S = self.spindle.S

        pos, vel = gen_3d_coords(self.spindle, letter, S)

        Organite.__init__(self, name, self.spindle.center, S, pos, vel)
        self.lbda = symbols('lambda_{}{}{}'.format(*self.site_id))
        self.rho = symbols('rho_{}{}{}'.format(*self.site_id))
        self.pi = symbols('pi_{}{}{}'.format(*self.site_id))

        self.pi_to_lr = {self.pi: -self.lbda + self.rho}
        self.lr_to_pi = {self.lbda: self.pi * (self.pi - 1)/2,
                         self.rho: self.pi * (self.pi + 1)/2}

        self.spindle.points.append(self.point)
        self.kt_MT_forces()

    def kt_MT_forces(self):
        F_k = parameters['F_k']
        V_k = parameters['V_k']
        kt_pull_L = LinearFV.new(self.spindle, self.spindle.S,
                                 self.spindle.spbL.point, self.point,
                                 self.lbda * F_k, V_k, -1)
        kt_pull_R = LinearFV.new(self.spindle, self.spindle.S,
                                 self.point, self.spindle.spbR.point,
                                 self.rho * F_k, V_k, -1)
