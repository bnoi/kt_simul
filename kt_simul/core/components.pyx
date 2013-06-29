"""
Simulated objects declaration. 'Objects' are for example: centromeres,
kinetochores, spindle pole body, etc.
"""

import logging

import random
import numpy as np
cimport numpy as np
cimport cython
from cpython cimport bool

__all__ = ["Spb", "Chromosome",
           "Centromere", "PlugSite", "Spindle"]

cdef class Spindle(object):
    def __init__(self, KD):
        self.KD = KD
        self.all_plugsites = []

    def get_all_plugsites(self):
        plugsites = []
        cdef Chromosome ch
        cdef PlugSite plugsite
        for ch in self.KD.chromosomes:
            for plugsite in ch.cen_A.plugsites:
                plugsites.append(plugsite)
            for plugsite in ch.cen_B.plugsites:
                plugsites.append(plugsite)
        return plugsites


cdef class Organite(object):
    """
    Base class for all the physical elements of the spindle

    Parameters
    ----------
    parent : an other subclass of :class:`Organite`
        from which the parameters are inheritated.
    init_pos : float, initial position

    Attributes
    ----------
    KD : a :class:`~spindle_dynamics.KinetoDynamics` instance
    pos : float, the position
    traj : ndarrat, the trajectory

    Methods
    -------
    set_pos(pos, time_point) : sets the position and updates the trajectory
    get_pos(time_point): returns the position at `time_point`
    """

    def __init__(self, parent, init_pos):
        self.parent = parent
        self.KD = parent.KD
        self.num_steps = parent.KD.num_steps
        self.traj = np.zeros(self.num_steps)
        self.pos = init_pos
        self.traj[0] = init_pos

    cdef void set_pos(self, float pos, int time_point=-1):
        """
        Sets the position. If `time_point` is provided, sets
        the corresponding value in `self.traj[time_point]`
        """
        self.pos = pos
        if pos > self.KD.spbR.pos:
            self.pos = self.KD.spbR.pos
        elif self.pos < self.KD.spbL.pos:
            self.pos = self.KD.spbL.pos
        if time_point >= 0:
            self.traj[time_point] = self.pos

    cdef float get_pos(self, int time_point=-1):
        """Returns the position.

        If `time_point` is -1 (default), returns the current position
        If `time_point` is not null, returns the position at time_point
        """
        if time_point == 0:
            return self.pos
        return self.traj[time_point]


cdef class Spb(Organite):
    """
    A spindle pole object.

    Attributes:
    ===========

    side: 1 or -1
    1 corresponds to the right SPB, i.e. to the Spb with x > 0.
    The left SPB, with `side = -1` corresponds to the daggered
    variable)
    """

    def __init__(self, spindle, side, L0):
        """
        side = 1 : left
        side = -1 : right
        L0 : spindle length
        """
        self.side = side
        init_pos = side * L0 / 2.
        Organite.__init__(self, spindle, init_pos)

cdef class Chromosome(Organite):
    """
    The chromosome, containing two centromeres ('A' and 'B')

    Parameters
    ----------
    spindle : :class:`~Spindle` instance

    """
    def __init__(self, spindle, ch_id):

        self.ch_id = ch_id
        self.id = self.ch_id
        d0 = spindle.KD.params['d0']
        L0 = spindle.KD.params['L0']
        center_pos = random.gauss(0, 0.2 * (L0 - d0))
        Organite.__init__(self, spindle, center_pos)
        self.cen_A = Centromere(self, 'A')
        self.cen_B = Centromere(self, 'B')
        self.correct_history = np.zeros((self.KD.num_steps, 2))
        self.correct_history[0] = self.correct()
        self.erroneous_history = np.zeros((self.KD.num_steps, 2))
        self.erroneous_history[0] = self.erroneous()

    cdef bool is_right_A(self):
        """
        Returns True if the majority of the attachments corresponds
        to centromere A bound to the right pole and centromere B bound
        to the left.
        """
        cdef int right_A, right_B
        right_A = self.cen_A.right_plugged() + self.cen_B.left_plugged()
        left_A =  self.cen_A.left_plugged() + self.cen_B.right_plugged()
        if right_A >= left_A:
            return True
        return False

    cdef int delta(self):
        """
        In case the centromeres swap (exchange side), the direction
        of the cohesin restoring force needs to be changed
        """
        return 1 if self.cen_A.pos < self.cen_B.pos else -1

    def correct(self):
        """
        returns the number of *correctly* plugged MTs
        """
        if self.is_right_A():
            return self.cen_A.right_plugged(), self.cen_B.left_plugged()
        else:
            return self.cen_A.left_plugged(), self.cen_B.right_plugged()

    def calc_erroneous_history(self):

        leftA, rightA = self.cen_A.calc_plug_history()
        leftB, rightB = self.cen_B.calc_plug_history()
        erroneous_hist = []
        for lA, rA, lB, rB in zip(leftA, rightA, leftB, rightB):
            erroneous = (rA, lB) if lA + rB > rA + lB else (lA, rB)
            erroneous_hist.append(erroneous)
        self.erroneous_history =  np.array(erroneous_hist)

    def calc_correct_history(self):
        leftA, rightA = self.cen_A.calc_plug_history()
        leftB, rightB = self.cen_B.calc_plug_history()
        correct_hist = []
        for lA, rA, lB, rB in zip(leftA, rightA, leftB, rightB):
            correct = (lA, rB) if lA + rB > rA + lB else (rA, lB)
            correct_hist.append(correct)
        self.correct_history = np.array(correct_hist)

    def erroneous(self):
        """
        Returns the number of *erroneously* plugged MTs
        """
        if self.is_right_A():
            return self.cen_A.left_plugged(), self.cen_B.right_plugged()
        else:
            return self.cen_A.right_plugged(), self.cen_B.left_plugged()

    cdef float pair_dist(self):
        return abs(self.cen_A.pos - self.cen_B.pos)

    cdef float plug_dist(self, float plugpos):
        return abs(self.center() - plugpos)

    cdef float center(self):
        return (self.cen_A.pos + self.cen_B.pos) / 2

    cdef np.ndarray center_traj(self):
        return (self.cen_A.traj + self.cen_B.traj) / 2

    cdef bool at_rightpole(self, float tol):
        """
        tol : tolerance distance
        """
        if self.cen_A.at_rightpole(tol) or self.cen_B.at_rightpole(tol):
            return True
        return False

    cdef bool at_leftpole(self, float tol):
        if self.cen_A.at_leftpole(tol) or self.cen_B.at_leftpole(tol):
            return True
        return False

    cdef bool at_pole(self, int side=0, float tol=0.01):
        if side == 1 and self.at_rightpole(tol):
            return True
        elif side == -1 and self.at_leftpole(tol):
            return True
        elif self.at_rightpole(tol) and self.at_leftpole(tol):
                return True
        else:
            return False


cdef class Centromere(Organite):
    """
    The centromere is where the plugsites are bound to the
    chromosome and where the cohesin spring restoring force, as
    well as the friction coefficient, are applied.
    This is a subclass of :class:`Organite`.

    Parameters
    ----------
    chromosome: a :class:`~Chromosome` instance
        the parent chromosome
    tag: {'A', 'B'}
       Side of the centromere. Note that the centromere
       side and the SPB side are not necesseraly related

    """
    def __init__(self, chromosome, tag):

        self.tag = tag
        self.chromosome = chromosome
        d0 = self.chromosome.KD.params['d0']
        if tag == 'A':
            init_pos = chromosome.pos - d0 / 2.
        elif tag == 'B':
            init_pos = chromosome.pos + d0 / 2.
        else:
            raise ValueError("the `tag` attribute must be 'A' or 'B'.")
        Organite.__init__(self, chromosome, init_pos)
        Mk = int(self.KD.params['Mk'])
        self.toa = 0  # time of arrival at pole
        self.plug_vector = np.zeros(Mk, dtype=np.int)
        self.plugsites = []
        cdef PlugSite ps
        for m in range(Mk):
            ps = PlugSite(self, m)
            self.plugsites.append(ps)
        self.calc_plug_vector()

    def is_attached(self):
        """
        Returns True if at least one plugsite is attached
        to at least one SPB
        """
        cdef PlugSite plugsite
        for plugsite in self.plugsites:
            if plugsite.plug_state != 0:
                return True
        return False

    cdef void calc_plug_vector(self):
        cdef np.ndarray[ITYPE_t] state
        cdef PlugSite plugsite
        state = np.array([plugsite.plug_state for plugsite
                          in self.plugsites])
        self.plug_vector = state

    def calc_plug_history(self):
        cdef np.ndarray[ITYPE_t, ndim = 2] state_hist
        cdef np.ndarray[ITYPE_t] right_hist, left_hist
        cdef PlugSite plugsite
        state_hist = np.array([plugsite.state_hist for plugsite in self.plugsites])
        right_hist = np.array(state_hist > 0).sum(axis=0)
        left_hist = np.array(state_hist < 0).sum(axis=0)
        return left_hist, right_hist

    cdef float P_attachleft(self):
        cdef float orientation
        orientation = self.KD.params['orientation']
        if orientation == 0:
            return 0.5
        cdef int lp, rp
        self.calc_plug_vector()
        lp = self.left_plugged()
        rp = self.right_plugged()
        if lp + rp == 0:
            return 0.5
        cdef float P_left
        P_left = 0.5 + orientation * (lp - rp) / (2 * (lp + rp))
        return P_left

    cdef int left_plugged(self):
        cdef int lp
        cdef np.ndarray[ITYPE_t] left_plugged
        left_plugged = self.plug_vector * (self.plug_vector - 1) / 2
        lp = left_plugged.sum()
        return lp

    cdef int right_plugged(self):
        cdef int rp
        cdef np.ndarray[ITYPE_t] right_plugged
        right_plugged = self.plug_vector * (1 + self.plug_vector) / 2
        rp = right_plugged.sum()
        return rp

    cdef bool at_rightpole(self, float tol=0.01):
        cdef float rightpole_pos
        rightpole_pos = self.KD.spbR.pos
        if abs(self.pos - rightpole_pos) < tol:
            return True
        return False

    cdef bool at_leftpole(self, float tol=0.01):
        cdef float leftpole_pos = self.KD.spbR.pos
        if abs(self.pos - leftpole_pos) < tol:
            return True
        return False

    def calc_toa(self, float tol=0.01):
        """
        Calculate time of arrivals
        """
        cdef np.ndarray dist_to_pole
        if self.KD.spbR.traj[-1] - self.traj[-1] < tol:
            dist_to_pole = self.KD.spbR.traj - self.traj
        elif self.traj[-1] - self.KD.spbL.traj[-1] < tol:
            dist_to_pole = self.traj - self.KD.spbL.traj
        else:
            self.toa = np.nan
            return False

        cdef int toa
        cdef float d
        toa = self.num_steps
        for d in dist_to_pole[::-1]:
            toa -= 1
            if d > tol:
                self.toa = toa * self.KD.params['dt']
                return True


cdef class PlugSite(Organite):
    """
    An attachment site object.

    Parameters:
    -----------
    centromere: a :class:`~Centromere` instance
    """

    def __init__(self, centromere, site_id):
        init_pos = centromere.pos
        Organite.__init__(self, centromere, init_pos)
        initial_plug = self.KD.initial_plug
        self.centromere = centromere
        self.tag = self.centromere.tag
        self.site_id = site_id

        if initial_plug == None:
            self.plug_state = random.randint(-1, 1)
        elif initial_plug == 'null':
            self.plug_state = 0
        elif initial_plug == 'amphitelic':
            self.plug_state = - 1 if self.tag == 'A' else 1
        elif initial_plug == 'random':
            self.plug_state = random.randint(-1, 1)
        elif initial_plug == 'monotelic':
            self.plug_state = - 1 if self.tag == 'A' else 0
        elif initial_plug == 'syntelic':
            self.plug_state = 1
        elif initial_plug == 'merotelic':
            self.plug_state = random.choice([-1,1])
        else:
            self.plug_state = initial_plug

        self.set_pos(init_pos)
        self.state_hist = np.zeros(self.KD.num_steps, dtype=np.int)
        self.state_hist[:] = self.plug_state
        self.P_att = 1 - np.exp(- self.KD.params['k_a'])

    cdef void set_plug_state(self, int state, int time_point=-1):
        self.plug_state = state
        self.plugged = 0 if state == 0 else 1
        self.state_hist[time_point:] = state

    cdef float calc_ldep(self):
        """
        The force term will increase or decrease the force applied to Kt
        according to kt - spb distance.

        Force term is calculated with linear function: ax + b

        Because we want to keep a "correct mean behaviour" according to Fmz,
        we compute lbase to keep: ldep * distance + lbase = 1

        * ldep: the lenght dependance factor parameter
        * ldep_balance: mean distance btw spb and Kt (measured by microscopy)
        """
        cdef double ldep
        cdef double ldep_balance
        cdef double lbase
        cdef double mt_length
        cdef double pole_pos
        cdef double force_term

        ldep = self.KD.params['ldep']
        ldep_balance = self.KD.params['ldep_balance']

        # The mean
        lbase = (1 - ldep * ldep_balance)

        if ldep > (1 / ldep_balance):
            # Remove length dependant behaviour
            # We don't want lbase < 0
            return 1

        pole_pos = self.KD.spbR.pos * self.plug_state
        mt_length = abs(pole_pos - self.pos)

        force_term = ldep * mt_length + lbase

        return force_term

    cdef void plug_unplug(self, int time_point):
        cdef float dice, side_dice
        dice = random.random()
        # Attachment
        if self.plug_state == 0 and dice < self.P_att:
            side_dice = random.random()
            P_left = self.centromere.P_attachleft()
            if side_dice < P_left:
                self.set_plug_state(-1, time_point)
            else:
                self.set_plug_state(1, time_point)
        # Detachment
        elif dice < self.P_det():
            self.set_plug_state(0, time_point)

    def is_correct(self, int time_point=-1):
        """
        Returns True if the plugsite is plugged
        correctly, i.e. doesn't contribute to an
        attachment error
        """

        if time_point < 0:
            plug_state = self.plug_state
        else:
            plug_state = self.state_hist[time_point]

        if plug_state == 0:
            return False
        if self.tag == 'A':
            if self.centromere.chromosome.is_right_A():
                return True if plug_state == 1 else False
            else:
                return True if plug_state == -1 else False
        else:
            if self.centromere.chromosome.is_right_A():
                return True if plug_state == -1 else False
            else:
                return True if plug_state == 1 else False

    cdef float P_det(self):
        cdef float d_alpha, k_d0
        d_alpha = self.KD.params['d_alpha']
        k_d0 = self.KD.params['k_a']
        if d_alpha == 0: return k_d0
        cdef float dist
        dist = abs(self.pos -
                   (self.centromere.chromosome.cen_A.pos +
                    self.centromere.chromosome.cen_B.pos) / 2.)
        if dist == 0: return 1.
        k_dc = k_d0  *  d_alpha / dist
        if k_dc > 1e4: return 1.
        return 1 - np.exp(-k_dc)
