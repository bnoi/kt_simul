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

import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()},
                  reload_support=True)

from ..core.spindle_dynamics import KinetoDynamics
from ..io import ParamTree
from ..core import parameters
from ..utils.progress import pprogress
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
                 force_parameters=[]):

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

        self.KD = KinetoDynamics(params, initial_plug=initial_plug,
                                 prng=self.prng)
        dt = self.paramtree.absolute_dic['dt']
        duration = self.paramtree.absolute_dic['span']
        self.num_steps = int(duration / dt)
        self.KD.anaphase = False
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
        return np.arange(0, self.KD.duration, self.KD.dt)

    @property
    def time_anaphase(self):
        return self.analysis['real_t_A']

    @property
    def index_anaphase(self):
        return np.argwhere(self.time == self.time_anaphase)[0][0]

    def index_anaphase_before(self, t_shift):
        return np.argwhere(self.time == (self.time_anaphase - t_shift))[0][0]

    def simul(self, ablat=None, ablat_pos=0.):
        """
        The simulation main loop.

        Parameters
        ----------

        ablat : float, optional
            Timepoint at which ablation takes place. If None (default)
            no ablation is performed.

        """

        # Check is simulation has already be done
        if self.KD.simulation_done:
            raise SimulationAlreadyDone("""A simulation is already done on this
                instance. Please create another Metaphase instance
                to launch a new simulation.""")

        kappa_c = self.KD.params['kappa_c']

        if self.verbose:
            log.info('Running simulation')
        bef = 0
        log_anaphase_onset = False

        self.analysis['real_t_A'] = self.KD.params['t_A']

        for time_point in range(1, self.num_steps):

            progress = int((time_point * 100.0) / self.num_steps)

            if self.verbose and progress != bef:
                pprogress(int(progress))
                bef = progress

            # Ablation test
            if ablat == time_point:
                if self.verbose:
                    log.info("Performing ablation")
                self._ablation(time_point, pos=ablat_pos)

            # Anaphase transition ?
            if self._anaphase_test(time_point):
                if not log_anaphase_onset:
                    pprogress(-1)
                    if self.verbose:
                        log.info("Anaphase onset at %i / %i" %
                                 (time_point, self.num_steps))
                    log_anaphase_onset = True

            self.KD.one_step(time_point)
            # if time_point % 100 == 0:
            #     print self.KD.At_mat

        if self.verbose:
            pprogress(-1)

        if self.verbose:
            log.info('Simulation done')
        self.KD.params['kappa_c'] = kappa_c
        delay_str = "delay = %2d seconds" % self.delay
        self.report.append(delay_str)

        for ch in self.KD.chromosomes:
            ch.calc_correct_history()
            ch.calc_erroneous_history()
            ch.cen_A.calc_toa()
            ch.cen_B.calc_toa()

    def get_report(self, time=0):
        """
        Print simulation state about a specific time point

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
        report["spbR pos"] = round(self.KD.spbR.traj[time], 3)
        report["spbL pos"] = round(self.KD.spbL.traj[time], 3)
        report["separate"] = ""

        for ch in self.KD.chromosomes:
            chdict = collections.OrderedDict()
            chdict["correct"] = ch.correct_history[time]
            chdict["erroneous"] = ch.erroneous_history[time]
            for cent in [ch.cen_A, ch.cen_B]:
                cen_dict = collections.OrderedDict()
                cen_dict["position"] = round(cent.traj[time], 3)
                for site in cent.plugsites:
                    site_dict = collections.OrderedDict()
                    site_dict['position'] = round(site.traj[time], 3)
                    site_dict['Plug state'] = site.state_hist[time]

                    cen_dict["PlugSite %i" % site.site_id] = site_dict

                chdict["Centromere %s" % cent.tag] = cen_dict

            report["Chromosome %i" % ch.ch_id] = chdict
            report["separate"] = ""

        return pretty_dict(report)

    def _anaphase_test(self, time_point):
        """
        At anaphase onset, set the cohesin spring constent to 0 and
        self.KD.anaphase to True.

        Returns
        -------
        bool
            True if anaphase has been executed.

        """

        t_A = int(self.KD.params['t_A'])
        dt = self.KD.params['dt']
        t = time_point * dt
        if self.KD.anaphase:
            return True
        if t >= t_A and self._plug_checkpoint():
            self.delay = t - t_A
            self.analysis['real_t_A'] = t
            # Then we just get rid of cohesin
            self.KD.params['kappa_c'] = 0.
            self.KD.calc_B()
            nb_mero = self._mero_checkpoint()
            if nb_mero:
                s = ("There were %d merotelic MT at anaphase onset"
                     % nb_mero)
                self.report.append(s)
            self.KD.anaphase = True
            return True
        return False

    def _ablation(self, time_point, pos=None):
        """
        Simulates a laser ablation: detaches all the kinetochores
        and sets the midzone stall force Fmz to 0

        Parameters
        ----------

        pos : float or None, optional
            Position of the laser beam within the spindle

        """
        if pos is None:
            pos = self.KD.spbR.pos[0]
        if not self.KD.spbL.pos <= pos <= self.KD.spbR.pos:
            log.warning('Missed shot, same player play again!')
            return
        self.KD.params['Fmz'] = 0.
        self.KD.params['k_a'] = 0.
        self.KD.params['k_d0'] = 0.
        self.KD.A0_mat = self.KD.time_invariantA()

        for plugsite in self.KD.spindle.all_plugsites:
            if pos < plugsite.pos[0] and plugsite.plug_state == - 1:
                plugsite.set_plug_state(0, time_point)
            elif pos > plugsite.pos[0] and plugsite.plug_state == 1:
                plugsite.set_plug_state(0, time_point)

    def _plug_checkpoint(self):
        """
        If the spindle assembly checkpoint is active, returns True
        if all chromosomes are plugged by at least one kMT, False
        otherwise.

        """
        sac = self.KD.params['sac']
        if sac == 0:
            return True
        for ch in self.KD.chromosomes:
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
        for ch in self.KD.chromosomes:
            if np.any(ch.erroneous()):
                nb_mero += sum(ch.erroneous())
                # print "active checkpoint"
        return nb_mero

    def _mplate_checkpoint(self):
        """
        Returns
        -------
        bool
            Returns True if each kinetochore is in the proper half
            of the spindle
        """
        for ch in self.KD.chromosomes:
            ktR = ch.cen_A.pos[0]
            ktL = ch.cen_B.pos[0]
            if min(ktR, ktL) <= 0 and max(ktR, ktL) >= 0:
                return True
        return True

    def show(self):
        """
        Quickly show kymograph of current simulation.
        Matplotlib is required.
        """

        import matplotlib.pyplot as plt

        duration = self.KD.duration
        dt = self.KD.dt
        # steps = self.KD.num_steps
        times = np.arange(0, duration, dt)
        kts = self.KD.chromosomes

        h = len(kts) * 2 + 6
        fig = plt.figure(figsize=(14, h))

        mainrowspan = int(len(kts) * 2)
        tot = len(kts) * 2 + mainrowspan

        main_ax = plt.subplot2grid((tot, 1), (3, 0), rowspan=mainrowspan)
        axs = []
        for i in range(len(kts)):
            ax1 = plt.subplot2grid((tot, 1), (i, 0), rowspan=1, sharex=main_ax)
            ax2 = plt.subplot2grid((tot, 1), (i + len(kts) + mainrowspan, 0),
                                   rowspan=1, sharex=main_ax)
            axs.append((ax1, ax2))

        main_ax = self.get_kymo_plot(main_ax)

        for (ax1, ax2), color, kt in zip(axs, self.chrom_colors, kts):
            correctA = kt.correct_history.T[1]
            correctB = kt.correct_history.T[0]
            errA = kt.erroneous_history.T[1]
            errB = kt.erroneous_history.T[0]

            ax1.plot(times, correctA, color=color, alpha=0.8)
            ax1.plot(times, errA, color=color, alpha=0.8, ls='--')

            ax2.plot(times, correctB, color=color, alpha=0.8)
            ax2.plot(times, errB, color=color, alpha=0.8, ls='--')

            ax1.set_yticks(list(range(0, 5, 2)))
            ax2.set_yticks(list(range(0, 5, 2)))

            ax1.grid()
            ax2.grid()

            for lab in ax1.xaxis.get_majorticklabels():
                lab.set_visible(False)

            for lab in ax2.xaxis.get_majorticklabels():
                lab.set_visible(False)

        for lab in ax1.xaxis.get_majorticklabels():
            lab.set_visible(True)

        for lab in ax2.xaxis.get_majorticklabels():
            lab.set_visible(True)

        axs[-1][1].set_xlabel('Time (s)')
        main_ax.set_ylabel('Distance from center (um)')

        # fig.suptitle("Kymograph")

        fig.tight_layout()
        plt.show()

        return fig

    def show_kymo(self, lims=[-3, 3]):
        """
        """
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(12, 7))
        main_ax = plt.subplot(111)

        main_ax = self.get_kymo_plot(main_ax, lims=lims)

        main_ax.set_xlabel('Time (s)')
        main_ax.set_ylabel('Distance from center (um)')

        # fig.suptitle("Kymograph")

        fig.tight_layout()
        plt.show()

        return fig

    def get_kymo_plot(self, ax, lims=[-3, 3]):
        """
        """

        import matplotlib

        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import matplotlib.cm as cmx

        dt = self.KD.dt
        duration = self.KD.duration
        anaphase = self.KD.params['t_A']
        times = np.arange(0, duration, dt)
        spbA = self.KD.spbL.traj
        spbB = self.KD.spbR.traj
        kts = self.KD.chromosomes

        ax.plot(times, spbA, color='black')
        ax.plot(times, spbB, color='black')

        if len(self.chrom_colors) == len(kts):
            def colors(x): return self.chrom_colors[x]
        else:
            cm = plt.get_cmap('gist_rainbow')
            norm = colors.Normalize(vmin=0, vmax=len(kts))
            sm = cmx.ScalarMappable(norm=norm, cmap=cm)
            colors = sm.to_rgba

        for i, kt in enumerate(kts):
            ax.plot(times, kt.cen_A.traj, color=colors(i), alpha=0.8)
            ax.plot(times, kt.cen_B.traj, color=colors(i), alpha=0.8)
        ax.grid()

        for lab in ax.xaxis.get_majorticklabels():
            lab.set_visible(False)

        ax.axvline(anaphase, color='black')
        if lims:
            ax.set_ylim(*lims)

        return ax

    def kymo_figure(self, axis_min=5, size=(2*5, 5)):
        """
        """

        colors = ["black",  # SPB
                  "#004bff",  # Kt 1
                  "#009d1c",  # Kt 2
                  "#ff2f00"]  # Kt 3

        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=size)
        ax = plt.subplot(111)

        drawer = ax.plot

        dt = self.KD.dt
        duration = self.KD.duration
        anaphase = self.KD.params['t_A']
        times = np.arange(0, duration, dt)
        spbA = self.KD.spbL.traj
        spbB = self.KD.spbR.traj
        kts = self.KD.chromosomes

        times = times / 60

        for i, kt in enumerate(kts):
            ax.plot(times, kt.cen_A.traj, color=colors[i+1], lw=3)
            ax.plot(times, kt.cen_B.traj, color=colors[i+1], lw=3)

        # Draw SPB
        drawer(times, spbA, label="SPB A", color=colors[0], lw=3)
        drawer(times, spbB, label="SPB B", color=colors[0], lw=3)

        # Set axis limit
        ax.set_xlim(min(times), max(times))
        ax.set_ylim(-4, 4)

        import matplotlib

        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

        majorLocator = matplotlib.ticker.MultipleLocator(axis_min)
        ax.xaxis.set_major_locator(majorLocator)
        ax.minorticks_off()

        majorFormatter = matplotlib.ticker.FuncFormatter(lambda x, y: "")
        ax.xaxis.set_major_formatter(majorFormatter)
        majorFormatter = matplotlib.ticker.FuncFormatter(lambda x, y: "")
        ax.yaxis.set_major_formatter(majorFormatter)

        majorLocator = matplotlib.ticker.MultipleLocator(2)
        ax.yaxis.set_major_locator(majorLocator)

        for i in ax.spines.values():
            i.set_linewidth(2)
            i.set_color('black')

        ax.grid(b=True, which='major', color='#555555', linestyle='-', alpha=1)

        return fig

    def get_attachment_vector(self):
        """Get attachment states for all chromosomes and all timepoints.
        """

        att = []

        time = np.arange(0, self.KD.duration, self.KD.dt)
        anaphase = self.KD.params['t_A']
        index_anaphase = np.argwhere(time == anaphase)[0][0]

        for ch in self.KD.chromosomes:
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
        if hasattr(self, 'KD'):
            del self.KD
        if hasattr(self, 'analysis'):
            del self.analysis
        if hasattr(self, 'measuretree'):
            del self.measuretree
        if hasattr(self, 'paramtree'):
            del self.paramtree
