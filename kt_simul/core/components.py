import logging

import numpy as np
from scipy import stats

from vispy.color import get_colormap

from ..mecabio import Structure
from ..mecabio import Point

log = logging.getLogger(__name__)


class Spindle(Structure):

    def __init__(self, name, parameters,
                 initial_plug='null',
                 prng=None):

        Structure.__init__(self, name)

        self.params = parameters
        if not prng:
            self.prng = np.random.RandomState()
        else:
            self.prng = prng

        self.initial_plug = initial_plug

        L = self.params['L0']
        N = int(self.params['N'])

        self.duration = self.params['span']
        self.dt = self.params['dt']

        self.attributes_df['plug_state'] = np.nan
        self.attributes_df['plug_state'] = self.attributes_df['plug_state'].astype("float")
        self.plug_state_idx = 1

        self.spbL = self.add_point(idx=0, init_pos=[-L/2, 0, 0], color="gray")
        self.spbR = self.add_point(idx=1, init_pos=[L/2, 0, 0], color="gray")
        self.add_link(self.spbL, self.spbR)

        self.initial_plug = initial_plug
        self.chromosomes = []

        # Generate colors
        cmap = get_colormap("husl")
        colors = cmap.map(np.linspace(0, 1, 3))

        # Convert colors to html
        def _c(color):
            return "#{0:02x}{1:02x}{2:02x}".format(*np.round(color*255).astype('int'))
        colors = [_c(color) for color in colors]

        for color, n in zip(colors, range(N)):
            ch = Chromosome(n, self, color)
            self.chromosomes.append(ch)

        self.update_geometry()

    def all_plugsites(self):
        for ch in self.chromosomes:
            for cen in (ch.cen_A, ch.cen_B):
                for ps in cen.plugsites:
                    yield ps


class Chromosome():

    def __init__(self, idx, spindle, color=None):
        self.id = idx
        d0 = spindle.params['d0']
        L0 = spindle.params['L0']
        self.spindle = spindle

        self.center_pos = np.array([self.spindle.prng.normal(0, 0.2*(L0 - d0)),
                                    self.spindle.prng.normal(0, d0),
                                    self.spindle.prng.normal(0, d0)])
        # spindle.add_point(idx=idx, center_pos)

        self.cen_A = Centromere(self, 'A', color)
        self.cen_B = Centromere(self, 'B', color)
        self.spindle.add_link(self.cen_A, self.cen_B)

    def is_right_A(self):
        """
        Returns True if the majority of the attachments corresponds
        to centromere A bound to the right pole and centromere B bound
        to the left.
        """

        right_A = self.cen_A.right_plugged() + self.cen_B.left_plugged()
        left_A = self.cen_A.left_plugged() + self.cen_B.right_plugged()
        if right_A >= left_A:
            return True
        return False

    def correct(self):
        """
        returns the number of *correctly* plugged MTs
        """
        if self.is_right_A():
            return self.cen_A.right_plugged(), self.cen_B.left_plugged()
        else:
            return self.cen_A.left_plugged(), self.cen_B.right_plugged()

    def erroneous(self):
        """
        Returns the number of *erroneously* plugged MTs
        """
        self.cen_A.calc_plug_vector()
        self.cen_B.calc_plug_vector()
        if self.is_right_A():
            return self.cen_A.left_plugged(), self.cen_B.right_plugged()
        else:
            return self.cen_A.right_plugged(), self.cen_B.left_plugged()

    def calc_erroneous_history(self):

        leftA, rightA = self.cen_A.calc_plug_history()
        leftB, rightB = self.cen_B.calc_plug_history()
        erroneous_hist = []
        for lA, rA, lB, rB in zip(leftA, rightA, leftB, rightB):
            erroneous = (rA, lB) if lA + rB > rA + lB else (lA, rB)
            erroneous_hist.append(erroneous)
        self.erroneous_history = np.array(erroneous_hist)

    def calc_correct_history(self):
        leftA, rightA = self.cen_A.calc_plug_history()
        leftB, rightB = self.cen_B.calc_plug_history()
        correct_hist = []
        for lA, rA, lB, rB in zip(leftA, rightA, leftB, rightB):
            correct = (lA, rB) if lA + rB > rA + lB else (rA, lB)
            correct_hist.append(correct)
        self.correct_history = np.array(correct_hist)

    def pair_dist(self):
        return self.cenA.dist(self.cenB)

    def plug_dist(self, plugpos):
        return np.linalg.norm(self.center() - plugpos)

    def center(self):
        return (self.cen_A.pos + self.cen_B.pos) / 2

    def center_traj(self):
        return (self.cen_A.traj + self.cen_B.traj) / 2

    def at_rightpole(self, tol):
        """
        tol : tolerance distance
        """
        if self.cen_A.at_rightpole(tol) or self.cen_B.at_rightpole(tol):
            return True
        return False

    def at_leftpole(self, tol):
        if self.cen_A.at_leftpole(tol) or self.cen_B.at_leftpole(tol):
            return True
        return False

    def at_pole(self, side=0, tol=0.01):
        if side == 1 and self.at_rightpole(tol):
            return True
        elif side == -1 and self.at_leftpole(tol):
            return True
        elif self.at_rightpole(tol) and self.at_leftpole(tol):
                return True
        else:
            return False


class Centromere(Point):

    def __init__(self, chromosome, tag, color=None):
        self.tag = tag
        self.chromosome = chromosome
        self.spindle = chromosome.spindle
        Mk = int(self.spindle.params['Mk'])
        d0 = self.spindle.params['d0']
        n = self.chromosome.id
        if tag == 'A':
            init_pos = chromosome.center_pos - d0 / 2.
            idx = (2 * n) * (Mk + 1) + 2
        elif tag == 'B':
            init_pos = chromosome.center_pos + d0 / 2.
            idx = (2 * n + 1) * (Mk + 1) + 2
        else:
            raise ValueError("the `tag` attribute must be 'A' or 'B'.")

        Point.__init__(self, self.spindle, idx, init_pos=init_pos, color=color)

        self.toa = 0  # time of arrival at pole
        self.plug_vector = np.zeros(Mk, dtype=np.int)
        self.plugsites = []
        self.plugsites_idx = []
        for m in range(Mk):
            ps = PlugSite(self, m)
            self.plugsites.append(ps)
            self.plugsites_idx.append(ps.idx)
        self.calc_plug_vector()

    def is_attached(self):
        """
        Returns True if at least one plugsite is attached
        to at least one SPB
        """
        for plugsite in self.plugsites:
            if plugsite.plug_state != 0:
                return True
        return False

    def calc_plug_vector(self):
        state = self.spindle.attributes_df.values[self.plugsites_idx, self.spindle.plug_state_idx]
        self.plug_vector = state

    def calc_plug_history(self):
        state_hist = np.array([plugsite.state_hist
                               for plugsite in self.plugsites])
        right_hist = np.array(state_hist > 0).sum(axis=0)
        left_hist = np.array(state_hist < 0).sum(axis=0)
        return left_hist, right_hist

    def P_attachleft(self):
        orientation = self.spindle.params['orientation']
        if orientation == 0:
            return 0.5
        self.calc_plug_vector()
        lp = self.left_plugged()
        rp = self.right_plugged()
        if lp + rp == 0:
            return 0.5
        P_left = 0.5 + orientation * (lp - rp) / (2 * (lp + rp))
        return P_left

    def left_plugged(self):
        left_plugged = self.plug_vector * (self.plug_vector - 1) // 2
        lp = left_plugged.sum()
        return lp

    def right_plugged(self):
        right_plugged = self.plug_vector * (1 + self.plug_vector) // 2
        rp = right_plugged.sum()
        return rp

    def at_rightpole(self, tol=0.01):

        if self.dist(self.spindle.spbR) < tol:
            return True
        return False

    def at_leftpole(self, tol=0.01):

        if self.dist(self.spindle.spbL) < tol:
            return True
        return False

    def calc_toa(self, tol=0.01):
        """
        Calculate time of arrivals
        """
        if self.spindle.spbR.dist(self) < tol:
            dist_to_pole = np.linalg.norm(self.spindle.spbR.traj -
                                          self.traj, axis=1)
        elif self.spindle.spbL.dist(self) < tol:
            dist_to_pole = np.linalg.norm(self.spindle.spbL.traj -
                                          self.traj, axis=1)
        else:
            self.toa = np.nan
            return False

        toa = self.num_steps
        for d in dist_to_pole[::-1]:
            toa -= 1
            if d > tol:
                self.toa = toa * self.spindle.params['dt']
                return True


class PlugSite(Point):

    def __init__(self, centromere, site_id):
        self.spindle = centromere.spindle
        self.centromere = centromere
        idx = self.centromere.idx + site_id + 1
        d0 = self.spindle.params['d0']
        init_pos = centromere.pos + self.spindle.prng.normal(0, 0.1*d0, 3)

        Point.__init__(self, self.spindle, idx, init_pos=init_pos)

        self.spindle.add_link(self.centromere, self)
        self.spindle.add_link(self.spindle.spbL, self)
        self.spindle.add_link(self.spindle.spbR, self)
        self.tag = self.centromere.tag
        self.current_side = ""
        self.site_id = site_id
        self._initial_plug()

    def _initial_plug(self):
        initial_plug = self.spindle.initial_plug
        if initial_plug is None:
            plug_state = self.spindle.prng.choice([-1, 0, 1])
        elif initial_plug == 'null':
            plug_state = 0
        elif initial_plug == 'amphitelic':
            plug_state = - 1 if self.tag == 'A' else 1
        elif initial_plug == 'random':
            plug_state = self.spindle.prng.choice([-1, 0, 1])
        elif initial_plug == 'monotelic':
            plug_state = - 1 if self.tag == 'A' else 0
        elif initial_plug == 'syntelic':
            plug_state = 1
        elif initial_plug == 'merotelic':
            plug_state = self.spindle.prng.choice([-1, 1])
        else:
            plug_state = initial_plug
        self.plug_state = plug_state

    @property
    def plug_state(self):
        return self.structure.attributes_df.at[self.idx, "plug_state"]

    @plug_state.setter
    def plug_state(self, state):
        self.plugged = 0 if state == 0 else 1
        self.structure.attributes_df.loc[self.idx, 'plug_state'] = state

    @property
    def state_hist(self):
        return self.structure.attributes_hist.xs(
            self.idx, axis='major').T['plug_state']

    def calc_ldep(self):
        """
        The force term will increase or decrease the force applied to Kt
        according to kt - spb distance.

        Force term is calculated with linear function: ax + b

        Because we want to keep a "correct mean behaviour" according to Fmz,
        we compute lbase to keep: ldep * distance + lbase = 1

        * ldep: the lenght dependance factor parameter
        * ldep_balance: mean distance btw spb and Kt (measured by microscopy)
        """

        # Disable ldep if current time point lower than time_ldep
        if ('time_ldep' in self.spindle.params.keys() and
                'time_ldep_inject' in self.spindle.params.keys()):
            if (self.spindle.params['time_ldep_inject'] and
                    self.spindle.time_point <
                    self.spindle.params['time_ldep']):
                return 1
            if (not self.spindle.params['time_ldep_inject'] and
                    self.spindle.time_point >
                    self.spindle.params['time_ldep']):
                return 1

        ldep = self.spindle.params['ldep']
        ldep_balance = self.spindle.params['ldep_balance']

        if ldep == 0:
            return 1
        # The mean
        lbase = (1 - ldep * ldep_balance)
        if ldep > (1 / ldep_balance):
            # Remove length dependant behaviour
            # We don't want lbase < 0
            return 1

        if self.plug_state == 1:
            pole = self.spindle.spbR
        elif self.plug_state == -1:
            pole = self.spindle.spbL
        else:
            return 0

        mt_length = self.dist(pole)
        ldep_factor = ldep * mt_length + lbase
        return ldep_factor

    def calc_ldep_for_attachment(self):
        """Using Cauchy distribution
        """
        # Disable ldep if current time point lower than time_ldep
        if ('time_ldep' in self.spindle.params.keys() and
                'time_ldep_inject' in self.spindle.params.keys()):
            if (self.spindle.params['time_ldep_inject'] and
               self.spindle.time_point < self.spindle.params['time_ldep']):
                return 1
            if (not self.spindle.params['time_ldep_inject'] and
               self.spindle.time_point > self.spindle.params['time_ldep']):
                return 1

        gamma = float(self.spindle.params['ldep_for_attachment_gamma'])
        mu = float(self.spindle.params['ldep_for_attachment_mu'])
        # N_mt = int(self.spindle.params['ldep_for_attachment_N_mt'])

        if gamma == 0:
            return 1

        if self.current_side == 'right':
            spb = self.spindle.spbR
        elif self.current_side == 'left':
            spb = self.spindle.spbL
        else:
            raise Exception('Wrong side idiot !')

        mt_size = self.dist(spb)
        spb_size = self.spindle.spbR.dist(self.spindle.spbL)

        # Get vector of spindle asize
        x = np.linspace(0, spb_size, 20)

        # Get Cauchy/Lorentz distribution according to current spindle size
        cauchy_cdf = stats.cauchy.pdf(x, loc=mu, scale=gamma)

        # Now compute a k_a prefactor for a given MT size
        # Then we normalize this pre factor to 1
        factor = stats.cauchy.pdf(mt_size, loc=mu, scale=gamma) / cauchy_cdf.mean()
        return factor

    def plug_unplug(self):
        """
        """
        dice = self.spindle.prng.rand()

        side_dice = self.spindle.prng.rand()
        P_left = self.centromere.P_attachleft()

        if side_dice < P_left:
            # Attachment
            self.current_side = "left"
            k_a = self.spindle.params['k_a'] * self.calc_ldep_for_attachment()
            pa = 1 - np.exp(-k_a)

            if self.plug_state == 0 and dice < pa:
                self.plug_state = -1
            # Detachment
            elif dice < self.P_det():
                self.plug_state = 0
        else:
            # Attachment
            self.current_side = "right"
            k_a = self.spindle.params['k_a'] * self.calc_ldep_for_attachment()
            pa = 1 - np.exp(-k_a)

            if self.plug_state == 0 and dice < pa:
                self.plug_state = 1
            # Detachment
            elif dice < self.P_det():
                self.plug_state = 0

    def is_correct(self, time_point=-1):
        """
        Returns True if the plugsite is plugged
        correctly, i.e. doesn't contribute to an
        attachment error
        """

        if time_point < 0:
            plug_state = self.plug_state
        else:
            plug_state = self.state_hist.loc[time_point]

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

    def P_det(self):
        """
        """

        d_alpha = self.spindle.params['d_alpha']
        k_d0 = self.spindle.params['k_d0']

        if d_alpha == 0:
            return k_d0

        ch_center = self.centromere.chromosome.center()
        dist = np.linalg.norm(self.pos - ch_center)

        if dist == 0:
            return 1

        # Cauchy distribution
        # x0 = 0
        # k_dc = k_d0 * (2 / (np.pi * d_alpha)) / (1 + ((dist - x0)
        #      / (d_alpha / 2)) ** 2)

        # Inverse distribution
        k_dc = k_d0 * d_alpha / dist

        if k_dc > 1e4:
            return 1

        return 1 - np.exp(-k_dc)

    # def calc_attach_trans(self):
    #     """
    #     """

    #     tau_trans = 10

    #     state = self.state_hist
    #     try:
    #         index_t0 = np.where(state != self.plug_state)[0][-1]
    #     except IndexError:
    #         index_t0 = 0
    #     dt = (current_t_index - index_t0) * self.spindle.dt
    #     trans_modulation = 1 - np.exp(-dt * (1 / tau_trans))
    #     return trans_modulation
