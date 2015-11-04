"""
This module provides the core simulation functionalities.

See Gay et al. J. Cell Biol., 2012 http://dx.doi.org/10.1083/jcb.201107124
The original framework was adapted from:
Civelekoglu-Scholey et al. Biophys. J. 90(11), 2006
http://dx.doi.org/10.1529/biophysj.105.078691
"""

from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import logging
import numpy as np
import collections

import matplotlib
import matplotlib.pyplot as plt

import tqdm

from ..core.components import Spindle
from ..core.dynamics import SpindleModel
from ..core import parameters
from ..utils.progress import print_progress
from ..utils.format import pretty_dict


log = logging.getLogger(__name__)

__all__ = ["Metaphase"]


class SimulationAlreadyDone(Exception):
    pass


class Metaphase(object):

    """
    An instance of the Metaphase class is a wrapper around
    the whole simulation.

    Parameters
    ----------

    paramtree : :class:`~kt_simul.io.ParamTree` instance or None
        The paramtree contains the parameters for the simulation
        if paramtree is None, the parameters are read
        from the file paramfile. Defaults to None.

    measuretree : :class:`~kt_simul.io.MeasureTree` instance or None
        The measuretree contains the observed characteristics
        of the mitosis e.g. metaphase spindle elongation rate, etc.
        if measuretree is None, the measures are read from the file
        indicated by the measurefile argument. Defaults to None.

    initial_plug : string or None
        Defines globally the initial attachment states. This argument can have
        the following values:
            - 'null': all kinetochores are detached
            - 'amphitelic': all chromosmes are amphitelic
            - 'random': all attachement site can be bound to either pole or
            detected with equal prob.
            - 'monotelic': right kinetochores are attached to the same pole,
            left ones are detached
            - 'syntelic' : all kinetochores are attached to the same pole

    keep_same_random_seed : bool
        To launch simulations with same random state

    """

    RANDOM_STATE = None

    def __init__(self,  paramtree=None, measuretree=None,
                 initial_plug='random', verbose=False,
                 keep_same_random_seed=False,
                 force_parameters=[], name='sim'):

        # Enable or disable log console
        self.verbose = verbose
        log = logging.getLogger(__name__)
        if not self.verbose:
            log.disabled = True
        else:
            log.disabled = False

        if paramtree is None:
            self.paramtree = parameters.get_default_paramtree()
        else:
            self.paramtree = paramtree

        if measuretree is None:
            self.measuretree = parameters.get_default_measuretree()
        else:
            self.measuretree = measuretree

        parameters.reduce_params(self.paramtree, self.measuretree,
                                 force_parameters=force_parameters)

        log.info('Parameters loaded')

        if keep_same_random_seed:
            self.prng = self.__class__.get_random_state()
        else:
            self.prng = np.random.RandomState()

        params = self.paramtree.relative_dic
        # Reset explicitely the unit parameters to their
        # dimentionalized value
        params['Vk'] = self.paramtree.absolute_dic['Vk']
        params['Fk'] = self.paramtree.absolute_dic['Fk']
        params['dt'] = self.paramtree.absolute_dic['dt']

        self.spindle = Spindle(name, params,
                               initial_plug=initial_plug,
                               prng=self.prng)
        self.spindle.update_geometry()
        self.model = SpindleModel(self.spindle)

        dt = self.paramtree.absolute_dic['dt']
        duration = self.paramtree.absolute_dic['span']

        self.num_steps = int(duration / dt)
        self.model.anaphase = False
        self.timelapse = np.arange(0, duration, dt)
        self.report = []
        self.delay = -1
        self.observations = {}
        self.analysis = {}

        log.info('Simulation initialized')
        log.disabled = False

        self.chrom_colors = ["red",
                             "green",
                             "blue"]

    @classmethod
    def get_random_state(cls):
        """
        """
        prng = np.random.RandomState()
        if not cls.RANDOM_STATE:
            cls.RANDOM_STATE = prng.get_state()
        else:
            prng.set_state(cls.RANDOM_STATE)
        return prng

    def __str__(self):
        lines = []
        lines.append('Metaphase class')
        try:
            lines.append('Parameters:')
            for line in str(self.paramtree.relative_dic).split(','):
                lines.append(line)
        except AttributeError:
            pass
        try:
            lines.append('Measures:')
            for line in str(self.measuretree.absolute_dic).split(','):
                lines.append(line)
        except AttributeError:
            pass

        return '\n'.join(lines)

    @property
    def time(self):
        return np.arange(0, self.paramtree['span'], self.model.dt)

    @property
    def time_anaphase(self):
        return self.analysis['real_t_A']

    @property
    def index_anaphase(self):
        return np.argwhere(self.time == self.time_anaphase)[0][0]

    def index_anaphase_before(self, t_shift):
        return np.argwhere(self.time == (self.time_anaphase - t_shift))[0][0]

    def simul(self, ablat=None, ablat_pos=0., progress=False):
        """
        The simulation main loop.

        Parameters
        ----------

        ablat : float, optional
            Timepoint at which ablation takes place. If None (default)
            no ablation is performed.

        """

        # Check is simulation has already be done
        if self.model.simulation_done:
            raise SimulationAlreadyDone("""A simulation is already done on this
                instance. Please create another Metaphase instance
                to launch a new simulation.""")

        kappa_c = self.model.params['kappa_c']

        if self.verbose:
            log.info('Running simulation')
        log_anaphase_onset = False

        self.analysis['real_t_A'] = self.model.params['t_A']

        for time_point in tqdm.tqdm(range(1, self.num_steps), total=self.num_steps,
                                    disable=not progress):

            # Ablation test
            if ablat == time_point:
                if self.verbose:
                    log.info("Performing ablation")
                self._ablation(time_point, pos=ablat_pos)

            # Anaphase transition ?
            if self._anaphase_test(time_point):
                if not log_anaphase_onset:
                    print_progress(-1)
                    if self.verbose:
                        log.info("Anaphase onset at %i / %i" %
                                 (time_point, self.num_steps))
                    log_anaphase_onset = True
            self.model.one_step(time_point)

        if self.verbose:
            log.info('Simulation done')

        self.model.params['kappa_c'] = kappa_c
        delay_str = "delay = %2d seconds" % self.delay
        self.report.append(delay_str)

        for ch in self.spindle.chromosomes:
            ch.calc_correct_history()
            ch.calc_erroneous_history()
            ch.cen_A.calc_toa()
            ch.cen_B.calc_toa()

    def get_report(self, time=0):
        """Print simulation state about a specific time point
        """
        params = self.paramtree.relative_dic

        report = collections.OrderedDict()

        report["Total time (span)"] = params["span"]
        report["Time precision (dt)"] = params["dt"]
        report["Chromosomes (N)"] = params["N"]
        report["Kinetochores (Mk)"] = params["Mk"]
        report["separate"] = ""

        report["Current time"] = "%i\n" % ((time + 1) * float(params["dt"]))

        report["separate"] = ""
        report["spbR pos"] = np.round(self.spindle.spbR.traj[time], 3)
        report["spbL pos"] = np.round(self.spindle.spbL.traj[time], 3)
        report["separate"] = ""

        for ch in self.spindle.chromosomes:
            chdict = collections.OrderedDict()
            chdict["correct"] = ch.correct_history[time]
            chdict["erroneous"] = ch.erroneous_history[time]
            for cent in [ch.cen_A, ch.cen_B]:
                cen_dict = collections.OrderedDict()
                cen_dict["position"] = np.round(cent.traj[time], 3)
                for site in cent.plugsites:
                    site_dict = collections.OrderedDict()
                    site_dict['position'] = np.round(site.traj[time], 3)
                    site_dict['Plug state'] = site.state_hist[time]

                    cen_dict["PlugSite %i" % site.site_id] = site_dict

                chdict["Centromere %s" % cent.tag] = cen_dict

            report["Chromosome %i" % ch.ch_id] = chdict
            report["separate"] = ""

        return pretty_dict(report)

    def _anaphase_test(self, time_point):
        """
        At anaphase onset, set the cohesin spring constant to 0 and
        self.model.anaphase to True.

        Returns
        -------
        bool
            True if anaphase has been executed.

        """

        t_A = int(self.model.params['t_A'])
        dt = self.model.params['dt']
        t = time_point * dt

        if self.model.anaphase:
            return True

        if t >= t_A and self._plug_checkpoint():
            self.delay = t - t_A
            self.analysis['real_t_A'] = t

            # Then we just get rid of cohesin
            self.model.params['kappa_c'] = 0
            nb_mero = self._mero_checkpoint()
            if nb_mero:
                s = ("There were %d merotelic MT at anaphase onset"
                     % nb_mero)
                self.report.append(s)
            self.model.anaphase = True
            return True
        return False

    def _ablation(self, time_point, pos=None, radius=0.2):
        """
        Simulates a laser ablation: detaches all the kinetochores
        and sets the midzone stall force Fmz to 0

        Parameters
        ----------

        pos : float or None, optional
            Position of the laser beam within the spindle

        """
        if pos is None:
            pos = self.spindle.spbR.pos
        if not self.spindle.spbL.pos.x <= pos[0] <= self.spindle.spbR.pos.x:
            log.warning('Missed shot, same player play again!')
            return
        self.model.params['Fmz'] = 0.
        self.model.params['k_a'] = 0.
        self.model.params['k_d0'] = 0.
        self.model.A0_mat = self.model.time_invariantA()

        for plugsite in self.spindle.all_plugsites():
            if pos < plugsite.pos.x and plugsite.plug_state == - 1:
                plugsite.plug_state = 0
            elif pos > plugsite.pos and plugsite.plug_state == 1:
                plugsite.plug_state = 0

    def _plug_checkpoint(self):
        """If the spindle assembly checkpoint is active, returns True
        if all chromosomes are plugged by at least one kMT, False
        otherwise.
        """
        sac = self.model.params['sac']
        if sac == 0:
            return True

        for ch in self.spindle.chromosomes:
            if not ch.cen_A.is_attached() or not ch.cen_B.is_attached():
                return False
        return True

    def _mero_checkpoint(self):
        """
        Returns
        -------

        int
            The total number of merotellic kT
        """
        nb_mero = 0
        for ch in self.spindle.chromosomes:
            if np.any(ch.erroneous()):
                nb_mero += sum(ch.erroneous())
        return nb_mero

    def _mplate_checkpoint(self):
        """
        Returns
        -------
        bool
            Returns True if each kinetochore is in the proper half
            of the spindle
        """
        for ch in self.spindle.chromosomes:
            ktR = ch.cen_A.pos.x
            ktL = ch.cen_B.pos.x
            if min(ktR, ktL) <= 0 and max(ktR, ktL) >= 0:
                return True
        return True

    def show(self):
        """
        Quickly show kymograph of current simulation.
        Matplotlib is required.
        """

        fig = plt.figure(figsize=(20, 15))

        cm = matplotlib.cm.get_cmap('Set1')
        colors = [cm(1 * i / len(self.spindle.chromosomes)) for i in range(len(self.spindle.chromosomes))]

        max_pos = max(self.spindle.spbL.traj.values.flatten().max(),
                      self.spindle.spbR.traj.values.flatten().max())

        for i, coord in enumerate(['x', 'y', 'z']):

            ax = fig.add_subplot(2, 2, i+1)
            ax.plot(self.time, self.spindle.spbL.traj[coord], color="black", lw=2)
            ax.plot(self.time, self.spindle.spbR.traj[coord], color="black", lw=2)

            ax.set_xlabel("Time (s)", fontsize=18)
            ax.set_ylabel("{} (um)".format(coord), fontsize=18)

            for ch, color in zip(self.spindle.chromosomes, colors):
                ax.plot(self.time, ch.cen_A.traj[coord], color=color, lw=2)
                ax.plot(self.time, ch.cen_B.traj[coord], color=color, lw=2)

            ax.set_ylim(-max_pos, max_pos)
            ax.set_xlim(0, self.time[-1])

            ax.axvline(self.time_anaphase, lw=1, alpha=0.8, color="black")

        fig.tight_layout()

        return fig

    def get_attachment_vector(self):
        """Get attachment states for all chromosomes and all timepoints.
        """

        att = []

        for ch in self.spindle.chromosomes:
            ch.calc_correct_history()
            ch.calc_erroneous_history()

            # Return attachment history for all timepoints
            c_hist = ch.correct_history
            e_hist = ch.erroneous_history
            state = [self.get_attachment(np.concatenate((c, e)))
                     for c, e in zip(c_hist, e_hist)]

            att.append(state)

        return np.array(att).T

    def get_attachment(self, state):
        """
        Get attachment state name according to biological classification:
            - amphitelic returns 0
            - monotelic returns 1
            - syntelic returns 2
            - merotelic returns 3
            - unattached returns 4
            - error means no attachment were found (issue with code) returns 5
        """

        def amphitelic(state):
            return (state[0] > 0 and state[1] > 0 and
                    state[2] == 0 and state[3] == 0)

        def monotelic(state):
            return (state[0] > 0 and state[1] == 0 and
                    state[2] == 0 and state[3] == 0)

        def syntelic(state):
            return (state[0] > 0 and state[1] == 0 and
                    state[2] == 0 and state[3] > 0)

        def merotelic(state):
            return (state[0] > 0 and state[2] > 0)

        def unattached(state):
            return (state[0] == 0 and state[1] == 0 and
                    state[2] == 0 and state[3] == 0)

        permuted_state = [state[1], state[0], state[3], state[2]]

        if amphitelic(state) or amphitelic(permuted_state):
            return 0
        elif monotelic(state) or monotelic(permuted_state):
            return 1
        elif syntelic(state) or syntelic(permuted_state):
            return 2
        elif merotelic(state) or merotelic(permuted_state):
            return 3
        elif unattached(state) or unattached(permuted_state):
            return 4
        else:
            return 5

    def get_attachment_names(self):
        """
        """

        names = ['amphitelic', 'monotelic', 'syntelic',
                 'merotelic', 'unattached', 'error']
        return names

    def __del__(self):
        """
        """
        if hasattr(self, 'spindle'):
            del self.spindle
        if hasattr(self, 'model'):
            del self.model
        if hasattr(self, 'analysis'):
            del self.analysis
        if hasattr(self, 'measuretree'):
            del self.measuretree
        if hasattr(self, 'paramtree'):
            del self.paramtree
