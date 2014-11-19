# -*- coding: utf-8 -*-
"""
This module defines all the objects considered in the simulation,
It uses cython used for the computer intensive bits

* Chromosome, Spb, Spindle, PlugSite are the spindle components
* KinetoDynamics wraps all the simulation internals: forces, vectors,
  and the system of equations.
"""

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import random
import numpy as np
cimport cython
cimport numpy as np
from cpython cimport bool

from .components cimport Spindle, Spb, Chromosome, Centromere, PlugSite
from .components import Spindle, Spb, Chromosome, Centromere, PlugSite

__all__ = ["KinetoDynamics"]

RIGHT = 1
LEFT = -1
a = 0
b = 1

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef class KinetoDynamics(object):
    """
    This class wraps all the simulation internals.
    """
    cdef public Spindle spindle
    cdef public Spb spbR, spbL
    cdef public unicode initial_plug
    cdef public bool simulation_done
    cdef public object params
    cdef public int num_steps
    cdef public float duration
    cdef public float dt
    cdef public list chromosomes
    cdef public np.ndarray B_mat, At_mat, A0_mat, Bn_mat, A0n_mat
    cdef public bool anaphase
    cdef public list all_plugsites
    cdef public int time_point
    cdef public np.ndarray speeds, speeds_n
    cdef public object prng

    def __init__(self, parameters, initial_plug='null', prng=None):
        """
        KinetoDynamics instenciation method

        :param parameters: A dictionnary of parameters as obtained from a
            xml_handler.ParamTree instance
        :type parameters: ParamTree instance

        :param initial_plug: Defines globally the initial attachment states.
            This argument can have the following values:
                * 'null': all kinetochores are detached
                * 'amphitelic': all chromosmes are amphitelic
                * 'random': all attachement site can be bound to
                        either pole or deteched with equal prob.
                * 'monotelic': right kinetochores are attached to the same pole,
                           left ones are detached
                * 'syntelic' : all kinetochores are attached to the same pole
        :type initial_plug: string or None
        """

        if not prng:
            self.prng = np.random.RandomState()
        else:
            self.prng = prng

        self.params = parameters
        L0 = self.params['L0']
        N = int(self.params['N'])
        Mk = int(self.params['Mk'])
        self.duration = self.params['span']
        self.dt = self.params['dt']
        self.num_steps = int(self.duration / self.dt)
        self.spindle = Spindle(self)
        self.spbR = Spb(self.spindle, RIGHT, L0)  # right spb (RIGHT = 1)
        self.spbL = Spb(self.spindle, LEFT, L0)  # left one (LEFT = -1)
        self.initial_plug = initial_plug#.decode('UTF-8')
        cdef Chromosome ch
        self.chromosomes = []
        for n in range(N):
            ch = Chromosome(self.spindle, n)
            self.chromosomes.append(ch)
        cdef int dim = 1 + N * (1 + Mk) * 2
        self.B_mat = np.zeros((dim, dim), dtype=float)
        self.Bn_mat = np.zeros((dim, dim), dtype=float)
        self.calc_B()
        self.A0_mat, self.A0n_mat = self.time_invariantA()
        self.At_mat = np.zeros((dim, dim), dtype=float)
        self.simulation_done = False
        self.anaphase = False
        self.all_plugsites = self.spindle.get_all_plugsites()
        self.speeds = np.zeros(dim)

    cdef np.ndarray _concat_matrices(self, np.ndarray M0, np.ndarray M1):
        ''' Puts the two passed matrices in a bigger block matrix
        [[M0 0]
         [0 M1]]
        '''
        cdef int n0, n1
        n0 = M0.shape[0]
        n1 = M1.shape[0]
        cdef np.ndarray[DTYPE_t, ndim = 2] M
        M = np.zeros((n0+n1,n0+n1))
        try:
            M[:n0, :n0] = M0
            M[n0:, n0:] = M1
        except ValueError:
            raise('Only accepting square matrices')
        return M

    cdef _idx(self, int side, int n, int m=-1):
        """
        :return: The index dictionnary
        """
        cdef int Mk, N, idx
        Mk = int(self.params['Mk'])
        N = int(self.params['N'])
        ### If n == -1, returns the index of the centromere
        idx = (2 * n + side) * (Mk + 1) + 2 + m
        return idx

    def one_step(self, int time_point):
        """
        Elementary step :

            1. Solving the equation for this time point (solve()).
            2. Updating position according to new speeds (position_update())
        """
        self.time_point = time_point
        self._one_step(time_point)
        if time_point == (self.num_steps - 1):
            self.simulation_done = True
            self.reset_positions()

    cdef void _one_step(self, int time_point):
        if not self.anaphase:
            self.plug_unplug(time_point)
        self.solve()
        self.position_update(time_point)

    cdef solve(self):
        cdef np.ndarray[DTYPE_t, ndim = 1] X, normal, Xt, C, Cn, Ct, pos_dep
        cdef np.ndarray[DTYPE_t, ndim = 2] A, An, At, Bt
        X = self.get_pos_vector()
        normal = self.get_normal_vector()
        Xt = np.concatenate([X, normal[1:]])
        A = self.calc_A()
        An = self.A0n_mat
        At = self._concat_matrices(A, An[1:, 1:])
        Bt = self._concat_matrices(self.B_mat, self.Bn_mat[1:, 1:])
        C, Cn = self.calc_C()
        Ct = np.concatenate([C, Cn[1:]])
        pos_dep = np.dot(Bt, Xt) + Ct
        self.speeds = np.linalg.solve(At, -pos_dep)
        cdef int n0 = C.shape[0]
        self.speeds, self.speeds_n = self.speeds[:n0], self.speeds[n0:]

    cdef np.ndarray get_pos_vector(self):
        """
        :return: a vector of the positions of each components
        """
        cdef int N, M, n, m
        N = int(self.params['N'])
        Mk = int(self.params['Mk'])
        cdef np.ndarray[DTYPE_t, ndim=1] X
        X = np.zeros(1 + 2*N * ( Mk + 1 ))
        X[0] = self.spbR.pos
        cdef int an, bn, anm, bnm
        cdef Chromosome ch
        for n in range(N):
            an = self._idx(0, n)
            bn = self._idx(1, n)
            ch = self.chromosomes[n]
            X[an] = ch.cen_A.pos
            X[bn] = ch.cen_B.pos
            for m in range(Mk):
                anm = self._idx(0, n, m)
                bnm = self._idx(1, n, m)
                X[anm] = ch.cen_A.plugsites[m].pos
                X[bnm] = ch.cen_B.plugsites[m].pos
        return X

    cdef np.ndarray get_normal_vector(self):
        """
        :return: a vector of the positions of each components
        """
        cdef int N, M, n, m
        N = int(self.params['N'])
        Mk = int(self.params['Mk'])
        cdef np.ndarray[DTYPE_t, ndim=1] X
        X = np.zeros(1 + 2*N * ( Mk + 1 ))
        X[0] = self.spbR.normal
        cdef int an, bn, anm, bnm
        cdef Chromosome ch
        for n in range(N):
            an = self._idx(0, n)
            bn = self._idx(1, n)
            ch = self.chromosomes[n]
            X[an] = ch.cen_A.normal
            X[bn] = ch.cen_B.normal
            for m in range(Mk):
                anm = self._idx(0, n, m)
                bnm = self._idx(1, n, m)
                X[anm] = ch.cen_A.plugsites[m].normal
                X[bnm] = ch.cen_B.plugsites[m].normal
        return X

    def calc_A(self):
        """
        :return: the matrix containing the linear terms of the equation set
        :math:`\mathbf{A}\dot{X} + \mathbf{B}X + C = 0`
        """
        return self._calc_A()

    cdef _calc_A(self):
        cdef np.ndarray[DTYPE_t, ndim=2] A
        self.time_dependentA()
        A = self.A0_mat + self.At_mat
        return A

    cpdef time_invariantA(self):
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        cdef float muc = self.params['muc']
        cdef float muk = self.params['muk']
        cdef float mus = self.params['mus']
        cdef float mu_t = self.params['mu_t']
        cdef float Mu_t = self.params['Mu_t']
        cdef float Vmz = self.params['Vmz']
        cdef float Fmz = self.params['Fmz']
        cdef float muco = self.params['muco']
        cdef int dims = 1 + 2 * N * ( Mk + 1 )
        ### Position
        cdef np.ndarray[DTYPE_t, ndim=2] A0
        A0 = np.zeros((dims, dims))
        ### Normal
        cdef np.ndarray[DTYPE_t, ndim=2] A0n
        A0n = np.zeros((dims, dims))
        ### We'll strip the first col & row later [no interaction at the SPB]

        A0[0, 0] = - 2 * mus - 4 * Fmz / Vmz

        cdef int n, m
        cdef Chromosome ch
        for n in range(N):
            ch = self.chromosomes[n]
            an = self._idx(0, n)
            bn = self._idx(1, n)
            A0[an, an] = - Mk * muk - muc + muco
            A0[bn, bn] = - Mk * muk - muc - muco
            A0n[an, an] = - Mk * mu_t - Mu_t
            A0n[bn, bn] = - Mk * mu_t - Mu_t
            for m in range(Mk):
                anm = self._idx(0, n, m)
                bnm = self._idx(1, n, m)
                A0[anm, anm] = - muk
                A0[anm, an] = muk
                A0[an, anm] = muk
                A0[bnm, bnm] = - muk
                A0[bnm, bn] = muk
                A0[bn, bnm] = muk

                A0n[anm, anm] = - mu_t
                A0n[anm, an] = mu_t
                A0n[an, anm] = mu_t
                A0n[bnm, bnm] = - mu_t
                A0n[bnm, bn] = mu_t
                A0n[bn, bnm] = mu_t
        return A0, A0n

    cdef time_dependentA(self):
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        cdef float pi_nmA, pi_nmB
        cdef int pluggedA, pluggedB
        cdef int n, m, anm, bnm
        self.At_mat[0,0] = 0
        cdef Chromosome ch
        cdef PlugSite plugsite_A, plugsite_B
        for n in range(N):
            ch = self.chromosomes[n]
            for m in range(Mk):
                plugsite_A = ch.cen_A.plugsites[m]
                pi_nmA = plugsite_A.plug_state
                pluggedA = plugsite_A.plugged
                plugsite_B = ch.cen_B.plugsites[m]
                pi_nmB = plugsite_B.plug_state
                pluggedB = plugsite_B.plugged

                # Lenght dependance
                if not self.anaphase:
                    pi_nmA *= plugsite_A.calc_ldep()
                    pi_nmB *= plugsite_B.calc_ldep()

                    # pi_nmA *= plugsite_A.calc_attach_trans()
                    # pi_nmB *= plugsite_B.calc_attach_trans()

                #spbs diag terms:
                self.At_mat[0, 0] -= pluggedA + pluggedB
                anm = self._idx(0, n, m)
                bnm = self._idx(1, n, m)
                self.At_mat[anm, anm] = - pluggedA
                self.At_mat[0, anm] = pi_nmA
                self.At_mat[anm, 0] = pi_nmA
                #B side
                self.At_mat[bnm, bnm] = -  pluggedB
                self.At_mat[0, bnm] = pi_nmB
                self.At_mat[bnm, 0] = pi_nmB

    def calc_B(self):
        """
        :return: the matrix containing the linear terms of the equation set
        :math:`\mathbf{A}\dot{X} + \mathbf{B}X + C = 0`
        """
        return self._calc_B()

    cdef _calc_B(self):
        cdef float kappa_c, kappa_k, kappa_t, Kappa_t
        kappa_k = self.params['kappa_k']
        kappa_c = self.params['kappa_c']
        kappa_t = self.params['kappa_t']
        Kappa_t = self.params['Kappa_t']
        cdef np.ndarray[DTYPE_t, ndim=2] Bk, Bc
        Bk = self.kinetochore_B()
        Bc = self.cohesin_B()
        Bkn = self.kinetochore_B()
        Bcn = self.cohesin_B()
        ### maybe we should return this instead...
        self.B_mat = kappa_k * Bk + kappa_c * Bc
        self.Bn_mat = kappa_t * Bkn + Kappa_t * Bcn

    cdef kinetochore_B(self):
        cdef N, Mk, n, m
        Mk = int(self.params['Mk'])
        N = int(self.params['N'])
        cdef int dim = 1 + N * (1 + Mk) * 2
        cdef int an, bn, anm, bnm
        cdef np.ndarray[DTYPE_t, ndim = 2] Bk
        Bk = np.zeros((dim, dim), dtype=float)
        for n in range(N):
            an = self._idx(0, n)
            bn = self._idx(1, n)
            Bk[an, an] = - Mk
            Bk[bn, bn] = - Mk
            for m in range(Mk):
                anm = self._idx(0, n, m)
                bnm = self._idx(1, n, m)
                Bk[anm, anm] = - 1
                Bk[bnm, bnm] = - 1
                Bk[an, anm] = 1
                Bk[anm, an] = 1
                Bk[bn, bnm] = 1
                Bk[bnm, bn] = 1
        return Bk

    cdef cohesin_B(self):
        cdef int N, Mk
        Mk = int(self.params['Mk'])
        N = int(self.params['N'])
        cdef int dim = 1 + N * (1 + Mk) * 2
        cdef np.ndarray[DTYPE_t, ndim = 2] Bc
        Bc = np.zeros((dim, dim), dtype=float)
        cdef int an, bn
        for n in range(N):
            an = self._idx(0, n)
            bn = self._idx(1, n)
            Bc[an, an] = - 1.
            Bc[bn, bn] = - 1.
            Bc[bn, an] = 1.
            Bc[an, bn] = 1.
        return Bc

    def calc_C(self):
        """
        :return: the matrix containing the linear terms of the equation set
        :math:`\mathbf{A}\dot{X} + \mathbf{B}X + C = 0`
        """
        return self._calc_C()

    cdef _calc_C(self):
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        cdef float Fmz = self.params['Fmz']
        cdef float d0 = self.params['d0']
        cdef float kappa_c = self.params['kappa_c']
        cdef float Kappa_t = self.params['Kappa_t']
        cdef float pi_nmA, pi_nmB
        cdef PlugSite plugsite_A, plugsite_B

        cdef np.ndarray[DTYPE_t, ndim = 1] C
        cdef np.ndarray[DTYPE_t, ndim = 1] Cn

        C = np.zeros(1 + N * (1 + Mk) * 2, dtype="float")
        Cn = np.zeros(1 + N * (1 + Mk) * 2, dtype="float") #cut line 0 afterwards
        C[0] = 2 * Fmz
        cdef int n, m, delta1, an, bn
        cdef int pluggedA, pluggedB
        cdef Chromosome ch
        for n in range(N):
            ch = self.chromosomes[n]
            delta1 = ch.delta1()
            an = self._idx(0, n)
            bn = self._idx(1, n)
            C[an] = - delta1 * kappa_c * d0
            C[bn] = delta1 * kappa_c * d0
            Cn[an] = np.pi * Kappa_t
            Cn[bn] = -np.pi * Kappa_t
            for m in range(Mk):
                anm = self._idx(0, n, m)
                bnm = self._idx(1, n, m)
                plugsite_A = ch.cen_A.plugsites[m]
                pi_nmA = plugsite_A.plug_state
                pluggedA = plugsite_A.plugged
                plugsite_B = ch.cen_B.plugsites[m]
                pi_nmB = plugsite_B.plug_state
                pluggedB = plugsite_B.plugged

                if not self.anaphase:
                    # Lenght dependance
                    pi_nmA *= plugsite_A.calc_ldep()
                    pi_nmB *= plugsite_B.calc_ldep()

                    # pi_nmA *= plugsite_A.calc_attach_trans()
                    # pi_nmB *= plugsite_B.calc_attach_trans()

                C[0] -= pluggedA + pluggedB
                C[anm] = pi_nmA
                C[bnm] = pi_nmB
        return C, Cn

    cdef plug_unplug(self, int time_point):
        """
        Let's play dices ...
        """
        cdef PlugSite plugsite
        for plugsite in self.all_plugsites:
            plugsite.plug_unplug(time_point)

    cdef position_update(self, int time_point):
        """
        Given the speeds obtained by solving A.x = btot and calculated switch events
        """
        cdef int Mk = int(self.params['Mk'])
        cdef int N = int(self.params['N'])
        cdef double dt = self.params['dt']
        cdef double Vk = self.params['Vk']
        cdef np.ndarray[DTYPE_t] speeds, speeds_n
        speeds = self.speeds
        speeds_n = self.speeds_n

        speeds *= Vk * dt # Back to real space
        speeds_n *= dt

        self.spbR.set_pos(self.spbR.pos + speeds[0], time_point)
        self.spbL.set_pos(self.spbL.pos - speeds[0], time_point)

        cdef int n, m, an, anm, bn, bnm
        cdef float new_pos
        cdef Chromosome ch
        cdef PlugSite plugsite
        for n in range(N):
            an = self._idx(0, n)
            bn = self._idx(1, n)
            ch = self.chromosomes[n]
            ch.cen_A.set_pos(ch.cen_A.pos + speeds[an], time_point)
            ch.cen_B.set_pos(ch.cen_B.pos + speeds[bn], time_point)
            ch.cen_A.set_normal(ch.cen_A.normal + speeds_n[an-1], time_point)
            ch.cen_B.set_normal(ch.cen_B.normal + speeds_n[bn-1], time_point)

            for m in range(Mk):
                anm = self._idx(0, n, m)
                plugsite = ch.cen_A.plugsites[m]
                new_pos = plugsite.pos + speeds[anm]
                plugsite.set_pos(new_pos, time_point)
                new_norm = plugsite.normal + speeds_n[anm-1]
                plugsite.set_normal(new_norm, time_point)

                bnm = self._idx(1, n, m)
                plugsite = ch.cen_B.plugsites[m]
                new_pos = plugsite.pos + speeds[bnm]
                plugsite.set_pos(new_pos, time_point)
                new_norm = plugsite.normal + speeds_n[bnm-1]
                plugsite.set_normal(new_norm, time_point)



    cdef reset_positions(self):
        """
        When a simu is done, reset all positions to t = 0
        """

        self.spbR.pos = self.spbR.traj[0]
        self.spbL.pos = self.spbL.traj[0]

        cdef int Mk = int(self.params['Mk'])
        cdef int N = int(self.params['N'])
        cdef int n, m, an, anm, bn, bnm
        cdef float new_pos
        cdef Chromosome ch
        cdef PlugSite plugsite
        for n in range(N):
            an = self._idx(0, n)
            bn = self._idx(1, n)
            ch = self.chromosomes[n]
            ch.cen_A.pos = ch.cen_A.traj[0]
            ch.cen_B.pos = ch.cen_B.traj[0]
            for m in range(Mk):
                anm = self._idx(0, n, m)
                bnm = self._idx(1, n, m)
                plugsite = ch.cen_A.plugsites[m]
                plugsite.pos = plugsite.traj[0]
                plugsite = ch.cen_B.plugsites[m]
                plugsite.pos = plugsite.traj[0]
